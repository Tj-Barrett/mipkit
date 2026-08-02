import os, json, copy
import shlex
import subprocess
from functools import reduce

from rich.progress import Progress, TimeElapsedColumn, BarColumn, TextColumn
from scipy.spatial.distance import cdist

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, GetPeriodicTable, Descriptors
from rdkit.Geometry import Point3D
from rdkit.Chem.rdchem import BondType
from rdkit import rdBase
# our work
from MIPkit.constants.constants import (
    fm2smiles,
    decode_fms,
    sulfur_exception,
    atom2PT,
)
from MIPkit.utils.check_availability import check_availability
from MIPkit.utils.check_reaction import check_reaction
from MIPkit.utils.determine_cost import determine_recipe_cost
from MIPkit.utils.parser import sanitize_inputs
from MIPkit.utils.print_stats import (
    print_stats,
    print_mw_distributions,
    print_update,
    print_bondstats,
)
from MIPkit.utils.print_molecule import print_molecule
from MIPkit.utils.read_pdb import read_pdb
from MIPkit.run.run_acpype import acpype_do
from MIPkit.run.run_gmx import gmx_do
from MIPkit.run.run_obabel import (
    obabel_combine,
    obabel_do,
    obabel_initiators,
    obabel_solvent
)
from MIPkit.utils.utils import (
    isnumber,
    mol_to_dictionary,
    isRXN
)
from MIPkit.initiators.initiators_yaml import init_attackers


################################################################
#
# Reaction Helpers
#
# These are shared by the explicit initiation, explicit
# polymerization, and implicit polymerization pathways.
#
################################################################

def _count_bond_order(bonds):
    """Return the total bond order for a set of bonds."""
    total = 0
    for bond in bonds:
        if bond.GetBondType() == BondType.SINGLE:
            total += 1
        elif bond.GetBondType() == BondType.DOUBLE:
            total += 2
        elif bond.GetBondType() == BondType.TRIPLE:
            total += 3
    return total


def _fix_overvalent_bonds(mol, bond_sets):
    """
    For each set of bonds at a reaction site, if the running bond-order
    sum reaches 5, reset every bond in that set to SINGLE.

    Mirrors the original per-bond accumulation check so that the fix
    triggers as soon as the cumulative order hits exactly 5.
    """
    for bonds in bond_sets:
        running = 0
        for bond in bonds:
            if bond.GetBondType() == BondType.SINGLE:
                running += 1
            elif bond.GetBondType() == BondType.DOUBLE:
                running += 2
            elif bond.GetBondType() == BondType.TRIPLE:
                running += 3
            if running == 5:
                for b in bonds:
                    mol.GetBondWithIdx(b.GetIdx()).SetBondType(BondType.SINGLE)


def _assign_explicit_hs(mol, unreacted, unreacted_types, species_count, typed=False):
    """
    Assign explicit hydrogen counts to unreacted terminal atoms.

    When typed=True (vinyl-vinyl implicit pathway), 'search' atoms receive
    one extra hydrogen compared to 'receptor' atoms:
        search  : bond order 1 -> 3 H, 2 -> 2 H, 3 -> 1 H
        receptor: bond order 1 -> 2 H, 2 -> 1 H, 3 -> 0 H

    All other pathways (vinyl-radical, self-interaction, initiation) use
    receptor-style counts for every atom (typed=False).
    """
    for _t, idx in enumerate(unreacted):
        _type = (
            unreacted_types[_t]
            if unreacted_types and _t < len(unreacted_types)
            else "receptor"
        )
        _pb = _count_bond_order(mol.GetAtomWithIdx(idx).GetBonds())
        mol.GetAtomWithIdx(idx).SetNoImplicit(True)

        if typed and _type == "search":
            hs_map = {1: 3, 2: 2, 3: 1}
        else:
            hs_map = {1: 2, 2: 1, 3: 0}

        hs = hs_map.get(_pb)
        if hs is not None:
            mol.GetAtomWithIdx(idx).SetNumExplicitHs(hs)
        else:
            print(
                f" Bond Count too high. Valence Error likely on Atom {idx} in"
                f" X{species_count}. Please report this on Github."
            )


def _renumber_and_finalize_mol(mol, species_count):
    """
    Renumber atoms by atom-map number, stamp PDB residue info as UNL,
    and run a UFF geometry optimisation.  Returns the renumbered molecule.
    """
    new_order = tuple(
        zip(
            *sorted(
                [
                    (j, k)
                    for k, j in enumerate(
                        [a.GetAtomMapNum() for a in mol.GetAtoms()]
                    )
                ]
            )
        )
    )[1]
    mi = Chem.AtomPDBResidueInfo()
    mi.SetResidueName("UNL")
    [a.SetMonomerInfo(mi) for a in mol.GetAtoms()]
    mi.SetResidueNumber(species_count)
    [a.SetMonomerInfo(mi) for a in mol.GetAtoms()]
    mol = Chem.RenumberAtoms(mol, new_order)
    uff_status = AllChem.UFFOptimizeMolecule(mol, maxIters=1000)
    if uff_status == 1:
        print(f" Warning: UFF optimization did not converge for X{species_count}. Geometry may be distorted.")
    elif uff_status == -1:
        print(f" Warning: UFF has no parameters for one or more atoms in X{species_count}. Skipped geometry optimization.")
    return mol


def _update_polymer_indexes(
    mol,
    species_count,
    fmsmiles,
    fmterminals,
    polymers,
    search,
    targa,
    receptor,
    targb,
    buffer,
):
    """
    Register the newly formed polymer X{species_count} in fmsmiles,
    fmterminals, and polymers.

    search / targa   : terminal-atom pair for molecule A.  Pass None to skip.
    receptor / targb : terminal-atom pair for molecule B.  Pass None to skip.
    buffer           : atom-index offset of molecule B (0 for self-interactions).
    """
    smol = Chem.RemoveAllHs(mol)
    fmsmiles[f"X{species_count}"] = Chem.MolToSmiles(smol, canonical=True)
    # Stage 1: cache the heavy-atom-only Mol alongside the index remap below, so
    # later cycles can look this species up instead of re-reading X{n}.pdb from
    # disk. `smol` is exactly what "rdkitmol"/`match` indexes into.

    smatch = json.loads(smol.GetProp("_smilesAtomOutputOrder"))
    match = [int(s) for s in smatch]

    fmterminals[f"X{species_count}"] = {}

    if search is not None and targa is not None:
        msearch = match[search]
        mtarga  = match[targa]
        fmterminals[f"X{species_count}"][msearch] = mtarga
        fmterminals[f"X{species_count}"][mtarga]  = msearch

    if receptor is not None and targb is not None:
        mreceptor = match[receptor + buffer]
        mtargb    = match[targb    + buffer]
        fmterminals[f"X{species_count}"][mreceptor] = mtargb
        fmterminals[f"X{species_count}"][mtargb]    = mreceptor

    polymers[f"X{species_count}"] = {"rdkitmol": match, "mol": smol}


def _verify_atom_correspondence(mol, positions, atoms=None, context=""):
    """
    Sanity check that `positions` (and optionally `atoms`, a parsed PDB atom-name
    list) line up 1:1 with `mol`'s atom order before positions are applied via
    SetAtomPosition. The whole coordinate-transfer pattern in this file assumes
    the post-MD PDB atom order matches the RDKit-built order; this makes that
    assumption fail loudly instead of silently writing coordinates onto the
    wrong atoms.

    Atom-count mismatches are unambiguous and raise. Element-symbol mismatches
    are a best-effort check (the PDB-atom-name-to-element lookup isn't
    exhaustive for arbitrary functional-monomer atom names) and only warn.
    """
    n_mol = mol.GetNumAtoms()
    if len(positions) != n_mol:
        raise ValueError(
            f"Atom-order mismatch{f' ({context})' if context else ''}: "
            f"molecule has {n_mol} atoms but {len(positions)} coordinates were "
            f"provided. Refusing to apply positions positionally."
        )
    if atoms is not None:
        symbol_table = atom2PT()
        for _m, atm in enumerate(atoms):
            if _m >= n_mol:
                break
            expected = symbol_table.get(atm, atm[0] if atm else "?")
            actual = mol.GetAtomWithIdx(_m).GetSymbol()
            if expected != actual:
                print(
                    f"Warning: possible atom-order mismatch{f' ({context})' if context else ''} "
                    f"at index {_m}: pdb atom '{atm}' (~{expected}) vs molecule atom {actual}."
                )
                break


################################################################
#
# Main
#
################################################################

'''

We need to refactor with classes
-> molecules into class
-> 



'''

def complex(args):
    # Get arguments for complexation
    pargs = sanitize_inputs(args)

    # read functional monomer structure
    precomplex = read_pdb(pargs.cplx, hydrogens=False)
    totmon = len(precomplex)
    fms = []
    for mol in precomplex:
        fms.append(mol['name'])

    pargs.fms = fms

    check_availability(pargs)
    check_reaction(pargs)
    determine_recipe_cost(pargs)

    if pargs.gpu_id:
        gpustr = f"Accelerating with {pargs.gpu_id}"
    else:
        gpustr = ""

    print(
        "\nComplexation using {} cores.\n{}\nPrecomplex : {} \nTemplate : {}.".format(pargs.ncpu,gpustr,pargs.cplx,pargs.template_name),
    )

    ################################################################
    #
    # Main
    #
    ################################################################
    
    # initialize counter
    i = 0

    columns = [
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TimeElapsedColumn(),
    ]

    with Progress(*columns) as pb:
        # overall
        pb1 = pb.add_task("[steel_blue1]  Complexation", total=pargs.ncycles)

        ################################################################
        #
        # Main
        #
        ################################################################
        RDLogger.DisableLog("rdApp.*")

        bmbt_logs_dir = pargs.logdir
        name = f"bmbt-polymerization.log"
        statname = f"bmbt-poly-stats.txt"

        # --------------------------------------------------------------

        # embedding parameters
        # enforceChirality False is vital for macromolecules
        # https://github.com/rdkit/rdkit/issues/4765
        rdversion = rdBase.rdkitVersion.split('.')
        if len(rdversion) > 0 and int(rdversion[0]) == 2024 and int(rdversion[1]) == 9:
            params = AllChem.ETKDGv3()
            params.useRandomCoords = True
            params.maxAttempts = 10000
            params.enforceChirality = False
        elif len(rdversion) > 0 and int(rdversion[0]) == 2025:
            params = AllChem.ETKDGv3()
            params.useRandomCoords = True
            # params.maxAttempts = 10000
            params.enforceChirality = False
        else:
            print(
                f"Note: RDKit version {rdBase.rdkitVersion} has not been explicitly "
                "tested with MIPkit; falling back to generic embedding parameters."
            )
            params = AllChem.ETKDGv3()
            params.useRandomCoords = True
            params.maxAttempts = 10000
            params.enforceChirality = False

        # --------------------------------------------------------------

        defm = decode_fms()
        fmsmiles = fm2smiles()
        sulfurs = sulfur_exception()

        PT = GetPeriodicTable()

        if pargs.template_name:
            template = read_pdb(pargs.template_name, hydrogens=False)

        ################################################################
        #
        # Pre - Polymerization Minimization
        #
        ################################################################

        if args.min and not pargs.restart:
            # Process Precomplex
            min_complex = []

            for j, fm in enumerate(precomplex):
                """
                Determine polymers from pdb, determine fms from smiles
                --
                some of the polymer smiles were coming out weird which is why we are not using them

                """
                receptor_atoms = []
                hcount = []

                _fm = fm["name"]
                mol = Chem.MolFromSmiles(fmsmiles[defm[_fm]])
                mol = Chem.RemoveAllHs(mol)

                embed_success = -1
                while embed_success < 0:
                    embed_success = AllChem.EmbedMolecule(mol, params)
                    if embed_success == -1:
                        print("Embedding failure. Retrying...")
                Chem.SanitizeMol(mol)

                positions = precomplex[j]["coords"]

                _verify_atom_correspondence(
                    mol, positions, precomplex[j].get("atoms"), context=f"min pass, {_fm}"
                )

                conf = mol.GetConformer()

                """
                Set positions to those used in the pdb file

                Then add hydrogens, get locations, bonds, and types to write out

                """

                for _m in range(mol.GetNumAtoms()):
                    x, y, z = positions[_m]
                    conf.SetAtomPosition(_m, Point3D(x, y, z))

                #
                # Remove hydrogens from unreacted atoms.
                #

                for _i, idx in enumerate(receptor_atoms):
                    atom = mol.GetAtomWithIdx(idx)
                    atom.SetNoImplicit(True)
                    atom.SetNumExplicitHs(hcount[_i])

                mol = Chem.rdmolops.AddHs(mol, addCoords=True)

                Chem.SanitizeMol(mol)

                new_fm = mol_to_dictionary(mol, fm["name"], PT)
                min_complex.append(new_fm)

            # Obabel
            logname = f"obabel-minimization.log"
            out_complex = obabel_do("min", -1, min_complex, pargs, logname)

            # Acpype
            logname = f"acpype-minimization.log"
            acpype_do("min", -1, min_complex, out_complex, pargs, logname)

            # Combined
            logname = f"combination-minimization.log"
            obabel_combine("min",-1,out_complex, pargs, logname)

            # GMX
            logname = f"gmx-minimization.log"
            gmx_run = gmx_do()
            gmx_run.params('min',-1,-1, pargs, logname)
            incomplex, intemplate, ininitiator = gmx_run.run()

            pargs.intemplate = intemplate
            pargs.incomplex = incomplex
            pargs.ininitiator = ininitiator

            if pargs.clean:
                workdir_q = shlex.quote(pargs.workdir)
                subprocess.run(f"rm {workdir_q}rxn*", shell=True)
                subprocess.run(f"rm {workdir_q}RXN*", shell=True)
                subprocess.run(f"rm {workdir_q}\\#*", shell=True)

        ################################################################
        #
        # Polymerization
        #
        ################################################################

        print_stats(
            pargs.workdir, init=True, precomplex=pargs.cplx, template=pargs.template_name,
        )

        print_mw_distributions(
            pargs.workdir, init=True, precomplex=pargs.cplx, template=pargs.template_name,
        )

        print_bondstats(
            pargs.workdir,
            statname,
            init=True,
            precomplex=pargs.cplx,
            template=pargs.template_name,
        )

        species_count = 0
        polymers = {}
        fmterminals = {}

        #
        # does this need to be open?
        #
        # was nice for the beginning but probably not necessary
        #

        if pargs.restart:
            start_num = pargs.restart_number
            # write out recipe
            recipe = []
            for r in precomplex:
                recipe.append(r["name"])

        else:
            start_num = -1

            # write out recipe
            recipe = []
            for r in precomplex:
                recipe.append(defm[r["name"]])

        with open(os.path.join(bmbt_logs_dir, f"{name}"), "w+") as bmbt:
            bmbt.write(
                "Functional Monomer Polymerization\nEpitope : {}\nRecipe : {}\n".format(pargs.template_name,recipe)
            )
            print_update(0, 0, 0, 0, init=True)


            for i in range(start_num+1, pargs.ncycles):
                bmbt.write(
                    f"\n----------------------------------\nStep {i}\n----------------------------------\n\n"
                )

                if pargs.restart and i == start_num+1:
                    prefix = f"step{i-1}"
                    pb.update(pb1, advance=i)

                    workdir = pargs.workdir
                    
                    incomplex = f"{workdir}{prefix}-fullcomplex-prod.pdb"
                    pargs.incomplex = incomplex

                    if pargs.protein:
                        intemplate = f"{workdir}{prefix}-protein-prod.pdb"
                    elif pargs.template:
                        intemplate = f"{workdir}{prefix}-template-prod.pdb"
                    else:
                        intemplate = None
                    pargs.intemplate = intemplate
                        
                    if pargs.explicit:
                        ininitiator = f"{workdir}{prefix}-initiator-prod.pdb"

                        # -explicit on the restart command line gives the
                        # original starting counts, not what's left after
                        # however many cycles already reacted before the
                        # crash. Re-sync pargs.initiator (type -> remaining
                        # count) from what's actually still in the reloaded
                        # initiator file, using the same grouping the
                        # in-loop update below uses.
                        restart_initiators = read_pdb(
                            ininitiator, hydrogens=True, initiators=True
                        )
                        inits = {}
                        for init in restart_initiators:
                            inits[init['name']] = inits.get(init['name'], 0) + 1
                        pargs.initiator = inits
                        del restart_initiators
                    else:
                        ininitiator = None
                    pargs.ininitiator = ininitiator

                    fms = read_pdb(incomplex)
                    for fm in fms:
                        if fm['name'][0] == 'X':
                            mol = Chem.MolFromPDBFile(
                                pargs.compdir+f"{fm['name']}.pdb",
                                sanitize=True,
                                removeHs=True,
                                proximityBonding=False,
                            )
                            rstsmiles = (
                                Chem.MolToSmiles(mol, canonical=True)
                            )

                            smatch = json.loads(
                                mol.GetProp("_smilesAtomOutputOrder")
                            )
                            match = []
                            for s in smatch:
                                match.append(int(s))

                            # Stage 1: warm the in-memory Mol cache with this
                            # already-parsed, already-Hs-free mol (removeHs=True
                            # above), matching what "rdkitmol"/match indexes into,
                            # so the rest of the run doesn't need to re-read
                            # X{n}.pdb.
                            tmp = {"rdkitmol": match, "mol": mol}
                            polymers[fm['name']] = tmp

                            del rstsmiles, smatch, match, tmp

                    # resume numbering after the highest X-species already on
                    # disk, so a restarted run doesn't reuse (and overwrite)
                    # an existing polymer's id
                    if polymers:
                        species_count = max(
                            int(name[1:]) for name in polymers if name[0] == "X"
                        ) + 1

                # read functional monomer structure
                if args.min and i == 0:
                    if pargs.template_name:
                        precomplex = read_pdb(incomplex, hydrogens=False)
                        template = read_pdb(intemplate, hydrogens=False)
                    else:
                        precomplex = read_pdb(incomplex, hydrogens=False)

                elif i > 0:
                    if pargs.template_name:
                        precomplex = read_pdb(incomplex, hydrogens=False)
                        template = read_pdb(intemplate, hydrogens=False)
                    else:
                        precomplex = read_pdb(incomplex, hydrogens=False)

                fms = []
                for mol in precomplex:
                    fms.append(mol['name'])

                # If no minimization, add initiators at runtime in gromacs
                #
                # only supporting aps, azo, abn 
                initiator_atoms = {}
                initiator_coords = []
                initiator_molecules = []

                if (args.explicit and args.min) or (args.explicit and i>0):
                    initiator_molecules = read_pdb(ininitiator, hydrogens=True, initiators=True)
                    init_attack = init_attackers()
                    for j, inits in enumerate(initiator_molecules):
                        if inits['name'] in ['ABN', 'AZO', 'APS']:
                            atoms, num = init_attack.load(inits['name'])
                            init_name = inits['name']

                            for _j, atm in enumerate(inits['atoms']):
                                if atm in atoms and _j in num:
                                    cc = tuple(inits['coords'][_j])
                                    initiator_coords.append(inits['coords'][_j])
                                    initiator_atoms[tuple(cc)] = [j, init_name, _j, -1]

                else:
                    initiator_atoms = None
                    initiator_coords = None

                # get array of all heavy protein atoms for statistics
                #
                # goal is to find bonding distances -- < 3.0 A for on target, > 3.0 for off target polymerization
                if pargs.template_name:
                    template_coords = []
                    template_res = []
                    template_atom = []
                    for res in template:
                        for j, pt in enumerate(res["coords"]):
                            template_coords.append(pt)
                            template_res.append(res["name"])
                            template_atom.append(res["atoms"][j])

                # Deterimining available terminals and getting statsitics of previous step
                search_atoms = {}
                search_coords = []
                receptor_atoms = {}
                receptor_coords = []

                nmono = 0
                npoly = 0
                mws = []

                for j, fms in enumerate(precomplex):
                    if fms["name"][0:1] == "X":
                        fm = fms["name"]

                        unreacted = []
                        mol = Chem.MolFromPDBFile(
                            pargs.compdir+f"{fm}.pdb",
                            sanitize=True,
                            removeHs=False,
                            proximityBonding=False,
                        )
                        for atom in mol.GetAtoms():
                            if (
                                atom.GetAtomicNum() == 6
                                and atom.GetNumImplicitHs() == 1
                            ):
                                unreacted.append(atom.GetIdx())
                        del mol

                        mol = Chem.MolFromPDBFile(
                            pargs.compdir+f"{fm}.pdb",
                            sanitize=True,
                            removeHs=True,
                            proximityBonding=False,
                        )
                        npoly += 1

                    else:
                        fm = defm[fms["name"]]
                        mol = Chem.MolFromSmiles(fmsmiles[fm])
                        nmono += 1

                    AllChem.EmbedMolecule(mol, params)
                    Chem.SanitizeMol(mol)

                    # Blanking on actual terms
                    #
                    # CH2 === CH-R
                    # ^       ^
                    # search  receptor
                    search = []
                    receptor = []
                    bns = []

                    # if "vinyl" in rxntype:
                    for bn, bond in enumerate(mol.GetBonds()):
                        if bond.GetBondType() == BondType.DOUBLE:
                            a = bond.GetBeginAtomIdx()
                            b = bond.GetEndAtomIdx()

                            _a = mol.GetAtomWithIdx(a)
                            _b = mol.GetAtomWithIdx(b)

                            if (
                                _a.GetAtomicNum() == 6
                                and _b.GetAtomicNum() == 6
                                and len(_a.GetNeighbors()) == 1
                            ):
                                search.append(a)
                                receptor.append(b)

                                bns.append(bn)

                            elif (
                                _b.GetAtomicNum() == 6
                                and _a.GetAtomicNum() == 6
                                and len(_b.GetNeighbors()) == 1
                            ):
                                search.append(b)
                                receptor.append(a)
                                bns.append(bn)

                    receptor = list(receptor)
                    search = list(search)

                    # get unreacted polymer terminals from polymer[]

                    fmterminals[fm] = {}

                    for _t, t in enumerate(search):
                        fmterminals[fm][t] = receptor[_t]
                        fmterminals[fm][receptor[_t]] = t

                    # searchs search for receptor match

                    for _t, t in enumerate(search):
                        cc = fms["coords"][t]

                        search_atoms[tuple(cc)] = [j, fm, t, bns[_t]]
                        search_coords.append(tuple(cc))

                    # receptor and unreacted

                    for _t, t in enumerate(receptor):
                        cc = fms["coords"][t]

                        receptor_atoms[tuple(cc)] = [j, fm, t, bns[_t]]
                        receptor_coords.append(tuple(cc))

                    if fms["name"][0:1] == "X":
                        for _t, t in enumerate(unreacted):
                            cc = fms["coords"][t]

                            receptor_atoms[tuple(cc)] = [j, fm, t, "RADICAL"]
                            receptor_coords.append(tuple(cc))

                    # get weights
                    mws.append(Chem.Descriptors.ExactMolWt(mol))
                    del mol
                totmon = len(mws)
                
                print_update(i, nmono, npoly, mws)
                print_stats(pargs.workdir, i, nmono, npoly, mws)
                print_mw_distributions(pargs.workdir, i, totmon, mws)

                reacted_mols = []
                step_complex = []

                n_fm = []

                '''
                Polymerization
                
                If distance between available atoms is less than the cutoff
                find minimum, polymerize, fix bond order, then make a new molecule
                
                For explicit polymerization, check against locations of APS/AZO/etc + Polymers
                For implicit, polymerization check against probability, other polymers, and initiated polymers
                
                '''
                if pargs.explicit:
                    pop_initiators = []

                    for j, seed in enumerate(receptor_atoms):
                        '''
                        Explicit Reaction against initiators
                        '''
                        if (pargs.explicit and pargs.min) or (pargs.explicit and i>0):
                            init_dists = cdist([seed], initiator_coords) < pargs.cutoff_CO
                            init_dists = init_dists.ravel()

                            for d, _d in enumerate(init_dists):
                                # Initiating a monomer
                                pta = receptor_coords[j]
                                ptb = initiator_coords[d]

                                na, fma, targa, bna = receptor_atoms[pta]
                                nb, fmb, targb, bnb = initiator_atoms[ptb]

                                if _d and fma[0:1] != "X":

                                    ta = f"{fma} {na}"

                                    # not reacted and different molecules
                                    n_fm_set = set(n_fm)
                                    if  ta not in n_fm_set and nb not in pop_initiators:
                                        unreacted = []
                                        unreacted_types = []
                                        hcount = []

                                        # monomer
                                        if fma in sulfurs:
                                            mola = Chem.MolFromPDBFile(
                                                pargs.sulfurdir+f"/{fma}.pdb",
                                                sanitize=True,
                                                removeHs=True,
                                                proximityBonding=False,
                                            )
                                        else:
                                            mola = Chem.MolFromSmiles(fmsmiles[fma])

                                        embed_success = -1
                                        while embed_success < 0:
                                            embed_success = AllChem.EmbedMolecule(mola, params)
                                            if embed_success == -1:
                                                print("Embedding failure. Retrying...")
                                        mola = Chem.RemoveAllHs(mola)

                                        buffer = mola.GetNumAtoms()

                                        # initiator
                                        molb = Chem.MolFromPDBFile(
                                            pargs.initiatordir+f"/{init_name.lower()}.pdb",
                                            sanitize=True,
                                            removeHs=True,
                                            proximityBonding=False,
                                        )

                                        embed_success = -1
                                        while embed_success < 0:
                                            embed_success = AllChem.EmbedMolecule(molb, params)
                                            if embed_success == -1:
                                                print("Embedding failure. Retrying...")
                                        molb = Chem.RemoveAllHs(molb)

                                        # ---- Explicit initiation: monomer + initiator ----
                                        # change bond order (sniff test)
                                        bt1 = mola.GetBondWithIdx(bna).GetBondType()

                                        if bt1 == BondType.DOUBLE:
                                            mola.GetBondWithIdx(bna).SetBondType(BondType.SINGLE)

                                            n_fm.append(ta)
                                            pop_initiators.append(nb)

                                            bmbt.write(f"Initiating {fma} {targa} as X{species_count}\n")

                                            # set 3D coordinates from the pdb
                                            positionsa = precomplex[na]["coords"]
                                            atomsa     = precomplex[na]["atoms"]
                                            positionsb = initiator_molecules[nb]["coords"]
                                            atomsb     = initiator_molecules[nb]["atoms"]

                                            _verify_atom_correspondence(mola, positionsa, atomsa, context=f"initiation, {fma} (A)")
                                            _verify_atom_correspondence(molb, positionsb, atomsb, context=f"initiation, {fmb} (B)")

                                            confa = mola.GetConformer()
                                            confb = molb.GetConformer()

                                            for _m in range(mola.GetNumAtoms()):
                                                x, y, z = positionsa[_m]
                                                confa.SetAtomPosition(_m, Point3D(x, y, z))

                                            for _m in range(molb.GetNumAtoms()):
                                                x, y, z = positionsb[_m]
                                                confb.SetAtomPosition(_m, Point3D(x, y, z))

                                            # combine molecules and form the new bond
                                            # https://gist.github.com/dominiquesydow/399899aed83c5583cfabdfd18b3cd7dd
                                            merged = reduce(Chem.CombineMols, [mola, molb])
                                            merged_editable = Chem.EditableMol(merged)
                                            merged_editable.AddBond(
                                                targa, targb + buffer,
                                                order=Chem.rdchem.BondType.SINGLE,
                                            )
                                            combined_mol = merged_editable.GetMol()

                                            pba = combined_mol.GetAtomWithIdx(targa).GetBonds()
                                            _fix_overvalent_bonds(combined_mol, [pba])
                                            Chem.SanitizeMol(combined_mol)

                                            # make unreacted search terminal
                                            # (flipped: search neighbour acts as receptor)
                                            search = fmterminals[fma][targa]
                                            unreacted.append(search)

                                            editable = Chem.EditableMol(combined_mol)
                                            combined_mol = editable.GetMol()

                                            _assign_explicit_hs(combined_mol, unreacted, None, species_count)

                                            # set all other places in initiator to no hydrogens
                                            # init_attack[APS] gives (['O', 'O'], [3, 4])
                                            no_hs_here = []
                                            for site in init_attack.load(fmb)[1]:
                                                if site == targb:
                                                    continue
                                                else:
                                                    atom = combined_mol.GetAtomWithIdx(site + buffer)
                                                    atom.SetNoImplicit(True)
                                                    atom.SetNumExplicitHs(0)
                                                    atom.SetFormalCharge(-1)

                                            combined_mol = Chem.rdmolops.AddHs(combined_mol, addCoords=True)
                                            Chem.SanitizeMol(combined_mol)

                                            combined_mol = _renumber_and_finalize_mol(combined_mol, species_count)

                                            reacted = mol_to_dictionary(combined_mol, f"X{species_count}", PT)
                                            reacted_mols.append(na)
                                            step_complex.append(reacted)

                                            w = Chem.PDBWriter(pargs.compdir + f"X{species_count}.pdb")
                                            w.write(combined_mol)
                                            w.close()

                                            # update indexes (mol A only; initiator has no tracked terminal)
                                            _update_polymer_indexes(
                                                combined_mol, species_count,
                                                fmsmiles, fmterminals, polymers,
                                                search=search, targa=targa,
                                                receptor=None, targb=None,
                                                buffer=buffer,
                                            )

                                            if pargs.template_name:
                                                print_bondstats(
                                                    pargs.workdir, statname, i, species_count,
                                                    fma, positionsa, atomsa,
                                                    fmb, positionsb, atomsb,
                                                    template_coords, template_res, template_atom,
                                                )

                                            species_count += 1

                                            del (
                                                mola, molb, confa, confb,
                                                combined_mol, merged, merged_editable,
                                            )


                    # Update initiator file to not include any reacted ones
                    if len(pop_initiators) > 0:
                        step_initiator_molecules = []
                        for j, init in enumerate(initiator_molecules):
                            if j in pop_initiators:
                                continue
                            else:
                                step_initiator_molecules.append(init)

                        inits = {}
                        init_names = []
                        for init in step_initiator_molecules:
                            if init['name'] in init_names:
                                inits[init['name']] += 1
                            else:
                                inits[init['name']] = 1
                                init_names.append(init['name'])

                        initiator_molecules = step_initiator_molecules    

                        del step_initiator_molecules

                        pargs.initiator = inits  
                
                for j, seed in enumerate(search_atoms):
                    # Distance matrix between searcher and all recepetors
                    dists = cdist([seed], receptor_coords) < pargs.cutoff_CC
                    dists = dists.ravel()


                    '''
                    Implicit and Initiated monomers 

                    Implicit :
                        Mol A
                        - Search atom reacts with receptor atom of B
                        - Receptor atom is capped with hydrogen
                        - if pargs.nocap is true, recepetor is NOT capped, and can react
                        Mol B
                        - Receptor atom reacts with search atom of A
                        - Search atom is made a 'RADICAL'
                    
                    Explicit :
                        Mol A
                        - Search atom is a 'RADICAL' and reacts with receptor of B
                        Mol B
                        - Receptor atom reacts with search atom of A
                        - Search atom is made a 'RADICAL' 
                    '''
                    for d, _d in enumerate(dists):

                        pta = search_coords[j]
                        ptb = receptor_coords[d]

                        na, fma, targa, bna = search_atoms[pta]
                        nb, fmb, targb, bnb = receptor_atoms[ptb]

                        # Allow Reaction
                        allow_rxn = False
                        if pargs.implicit:
                            # Implicit
                            if isRXN(args, pargs.rxn_prob):
                                allow_rxn = True
                            else:
                                allow_rxn = False
                        
                        elif pargs.explicit and bnb == "RADICAL":
                            allow_rxn = True


                        # Proceed if allowed to react
                        if _d and allow_rxn:
                            ta = f"{fma} {na}"
                            tb = f"{fmb} {nb}"

                            # A is always search, B is always receptor
                            self_interaction = False
                            if ta == tb and fma[0:1] == "X" and fmb[0:1] == "X":
                                self_interaction = True

                            # not reacted and different molecules
                            n_fm_set = set(n_fm)
                            if (ta not in n_fm_set
                                and tb not in n_fm_set
                                and ta != tb):
                                unreacted = []
                                unreacted_types = []
                                hcount = []

                                # Initialize polymer atom mappings
                                matcha = None
                                matchb = None

                                # molecule A
                                if fma[0:1] == "X":
                                    mola = Chem.MolFromPDBFile(
                                        pargs.compdir+f"{fma}.pdb",
                                        sanitize=True,
                                        removeHs=False,
                                        proximityBonding=False,
                                    )
                                    embed_success = -1
                                    while embed_success < 0:
                                        embed_success = AllChem.EmbedMolecule(
                                            mola, params
                                        )
                                        if embed_success == -1:
                                            print("Embedding failure. Retrying...")

                                    # look for valence mismatch  -- polymerized locations will be missing one hydrogen
                                    for atom in mola.GetAtoms():
                                        if (
                                            atom.GetAtomicNum() == 6
                                            and atom.GetNumImplicitHs() == 1
                                            and atom.GetIdx() != targa
                                        ):
                                            unreacted.append(atom.GetIdx())
                                            unreacted_types.append("RADICAL")
                                            hc = 0
                                            for nbr in atom.GetNeighbors():
                                                if nbr.GetAtomicNum() == 1:
                                                    hc += 1
                                            hcount.append(hc)

                                    del mola
                                    mola = Chem.MolFromPDBFile(
                                        pargs.compdir+f"{fma}.pdb",
                                        sanitize=True,
                                        removeHs=True,
                                        proximityBonding=False,
                                    )
                                    if fma not in polymers:
                                        raise ValueError(
                                            f"Polymer {fma} not found in polymer registry. "
                                            f"Check if it was properly created in previous step."
                                        )
                                    matcha = polymers[fma]["rdkitmol"]
                                elif fma in sulfurs:
                                    mola = Chem.MolFromPDBFile(
                                        pargs.sulfurdir+f"/{fma}.pdb",
                                        sanitize=True,
                                        removeHs=True,
                                        proximityBonding=False,
                                    )
                                else:
                                    mola = Chem.MolFromSmiles(fmsmiles[fma])

                                embed_success = -1
                                while embed_success < 0:
                                    embed_success = AllChem.EmbedMolecule(mola, params)
                                    if embed_success == -1:
                                        print("Embedding failure. Retrying...")
                                mola = Chem.RemoveAllHs(mola)

                                buffer = mola.GetNumAtoms()

                                # molecule B
                                if fmb[0:1] == "X":
                                    molb = Chem.MolFromPDBFile(
                                        pargs.compdir+f"{fmb}.pdb",
                                        sanitize=True,
                                        removeHs=False,
                                        proximityBonding=False,
                                    )
                                    embed_success = -1
                                    while embed_success < 0:
                                        embed_success = AllChem.EmbedMolecule(
                                            molb, params
                                        )
                                        if embed_success == -1:
                                            print("Embedding failure. Retrying...")

                                    # look for valence mismatch  -- polymerized locations will be missing one hydrogen
                                    for atom in molb.GetAtoms():
                                        if (
                                            atom.GetAtomicNum() == 6
                                            and atom.GetNumImplicitHs() == 1
                                            and atom.GetIdx() != targb
                                        ):
                                            unreacted.append(atom.GetIdx() + buffer)
                                            unreacted_types.append("RADICAL")
                                            hc = 0
                                            for nbr in atom.GetNeighbors():
                                                if nbr.GetAtomicNum() == 1:
                                                    hc += 1
                                            hcount.append(hc)

                                    del molb

                                    molb = Chem.MolFromPDBFile(
                                        pargs.compdir+f"{fmb}.pdb",
                                        sanitize=True,
                                        removeHs=True,
                                        proximityBonding=False,
                                    )
                                    if fmb not in polymers:
                                        raise ValueError(
                                            f"Polymer {fmb} not found in polymer registry. "
                                            f"Check if it was properly created in previous step."
                                        )
                                    matchb = polymers[fmb]["rdkitmol"]
                                elif fmb in sulfurs:
                                    molb = Chem.MolFromPDBFile(
                                        pargs.sulfurdir+f"/{fmb}.pdb",
                                        sanitize=True,
                                        removeHs=True,
                                        proximityBonding=False,
                                    )
                                else:
                                    molb = Chem.MolFromSmiles(fmsmiles[fmb])

                                embed_success = -1
                                while embed_success < 0:
                                    embed_success = AllChem.EmbedMolecule(molb, params)
                                    if embed_success == -1:
                                        print("Embedding failure. Retrying...")
                                molb = Chem.RemoveAllHs(molb)

                                #
                                # Vinyl - vinyl bonding of implicit system
                                #

                                # ---- Vinyl-vinyl bonding (implicit or explicit radical) ----
                                if isnumber(bna) and isnumber(bnb):
                                    bt1 = mola.GetBondWithIdx(bna).GetBondType()
                                    bt2 = molb.GetBondWithIdx(bnb).GetBondType()

                                    if bt1 == BondType.DOUBLE and bt2 == BondType.DOUBLE:
                                        mola.GetBondWithIdx(bna).SetBondType(BondType.SINGLE)
                                        molb.GetBondWithIdx(bnb).SetBondType(BondType.SINGLE)

                                        # ta != tb is already guaranteed by the enclosing
                                        # "not reacted and different molecules" check above
                                        n_fm.append(ta)
                                        n_fm.append(tb)

                                        bmbt.write(
                                            f"Bonding {fma} {targa} and {fmb} {targb}"
                                            f" ({targb + buffer}) as X{species_count}\n"
                                        )

                                        # set 3D coordinates from the pdb
                                        positionsa = precomplex[na]["coords"]
                                        atomsa     = precomplex[na]["atoms"]
                                        positionsb = precomplex[nb]["coords"]
                                        atomsb     = precomplex[nb]["atoms"]

                                        _verify_atom_correspondence(mola, positionsa, atomsa, context=f"vinyl-vinyl, {fma} (A)")
                                        _verify_atom_correspondence(molb, positionsb, atomsb, context=f"vinyl-vinyl, {fmb} (B)")

                                        confa = mola.GetConformer()
                                        confb = molb.GetConformer()

                                        if fma[0:1] == "X":
                                            for _m in matcha:
                                                x, y, z = positionsa[_m]
                                                confa.SetAtomPosition(_m, Point3D(x, y, z))
                                        else:
                                            for _m in range(mola.GetNumAtoms()):
                                                x, y, z = positionsa[_m]
                                                confa.SetAtomPosition(_m, Point3D(x, y, z))

                                        if fmb[0:1] == "X":
                                            for _m in matchb:
                                                x, y, z = positionsb[_m]
                                                confb.SetAtomPosition(_m, Point3D(x, y, z))
                                        else:
                                            for _m in range(molb.GetNumAtoms()):
                                                x, y, z = positionsb[_m]
                                                confb.SetAtomPosition(_m, Point3D(x, y, z))

                                        # combine molecules and form the new bond
                                        # https://gist.github.com/dominiquesydow/399899aed83c5583cfabdfd18b3cd7dd
                                        merged = reduce(Chem.CombineMols, [mola, molb])
                                        merged_editable = Chem.EditableMol(merged)
                                        merged_editable.AddBond(
                                            targa, targb + buffer,
                                            order=Chem.rdchem.BondType.SINGLE,
                                        )
                                        combined_mol = merged_editable.GetMol()

                                        pba = combined_mol.GetAtomWithIdx(targa).GetBonds()
                                        pbb = combined_mol.GetAtomWithIdx(targb + buffer).GetBonds()
                                        _fix_overvalent_bonds(combined_mol, [pba, pbb])
                                        Chem.SanitizeMol(combined_mol)

                                        # assign terminal hydrogen counts
                                        # search neighbour gets +1 H vs receptor neighbour
                                        search   = fmterminals[fma][targa]
                                        receptor = fmterminals[fmb][targb]
                                        unreacted.append(search)
                                        unreacted_types.append("search")
                                        unreacted.append(receptor + buffer)
                                        unreacted_types.append("receptor")

                                        editable = Chem.EditableMol(combined_mol)
                                        combined_mol = editable.GetMol()

                                        _assign_explicit_hs(
                                            combined_mol, unreacted, unreacted_types,
                                            species_count, typed=True,
                                        )

                                        combined_mol = Chem.rdmolops.AddHs(combined_mol, addCoords=True)
                                        Chem.SanitizeMol(combined_mol)

                                        combined_mol = _renumber_and_finalize_mol(combined_mol, species_count)

                                        reacted = mol_to_dictionary(combined_mol, f"X{species_count}", PT)
                                        reacted_mols.append(na)
                                        reacted_mols.append(nb)
                                        step_complex.append(reacted)

                                        w = Chem.PDBWriter(pargs.compdir + f"X{species_count}.pdb")
                                        w.write(combined_mol)
                                        w.close()

                                        _update_polymer_indexes(
                                            combined_mol, species_count,
                                            fmsmiles, fmterminals, polymers,
                                            search=search,    targa=targa,
                                            receptor=receptor, targb=targb,
                                            buffer=buffer,
                                        )

                                        if pargs.template_name:
                                            print_bondstats(
                                                pargs.workdir, statname, i, species_count,
                                                fma, positionsa, atomsa,
                                                fmb, positionsb, atomsb,
                                                template_coords, template_res, template_atom,
                                            )

                                        species_count += 1

                                        del (
                                            mola, molb, confa, confb,
                                            combined_mol, merged, merged_editable,
                                        )

                                    else:
                                        continue

                                #
                                # Bonding a vinyl functional group to a free radical carbon
                                #

                                # ---- Vinyl-radical bonding (explicit propagation) ----
                                elif (isnumber(bna) and bnb == "RADICAL") or (
                                    isnumber(bnb) and bna == "RADICAL"
                                ):
                                    if isnumber(bna):
                                        bt = mola.GetBondWithIdx(bna).GetBondType()
                                    else:
                                        bt = molb.GetBondWithIdx(bnb).GetBondType()

                                    if bt == BondType.DOUBLE:
                                        if isnumber(bna):
                                            mola.GetBondWithIdx(bna).SetBondType(BondType.SINGLE)
                                        if isnumber(bnb):
                                            molb.GetBondWithIdx(bnb).SetBondType(BondType.SINGLE)

                                        n_fm.append(ta)
                                        n_fm.append(tb)

                                        bmbt.write(
                                            f"Bonding {fma} {targa} and {fmb} {targb}"
                                            f" ({targb + buffer}) as X{species_count}\n"
                                        )

                                        # set 3D coordinates from the pdb
                                        positionsa = precomplex[na]["coords"]
                                        atomsa     = precomplex[na]["atoms"]
                                        positionsb = precomplex[nb]["coords"]
                                        atomsb     = precomplex[nb]["atoms"]

                                        _verify_atom_correspondence(mola, positionsa, atomsa, context=f"vinyl-radical, {fma} (A)")
                                        _verify_atom_correspondence(molb, positionsb, atomsb, context=f"vinyl-radical, {fmb} (B)")

                                        confa = mola.GetConformer()
                                        confb = molb.GetConformer()

                                        if fma[0:1] == "X":
                                            for _m in matcha:
                                                x, y, z = positionsa[_m]
                                                confa.SetAtomPosition(_m, Point3D(x, y, z))
                                        else:
                                            for _m in range(mola.GetNumAtoms()):
                                                x, y, z = positionsa[_m]
                                                confa.SetAtomPosition(_m, Point3D(x, y, z))

                                        if fmb[0:1] == "X":
                                            for _m in matchb:
                                                x, y, z = positionsb[_m]
                                                confb.SetAtomPosition(_m, Point3D(x, y, z))
                                        else:
                                            for _m in range(molb.GetNumAtoms()):
                                                x, y, z = positionsb[_m]
                                                confb.SetAtomPosition(_m, Point3D(x, y, z))

                                        # combine molecules and form the new bond
                                        # https://gist.github.com/dominiquesydow/399899aed83c5583cfabdfd18b3cd7dd
                                        merged = reduce(Chem.CombineMols, [mola, molb])
                                        merged_editable = Chem.EditableMol(merged)
                                        merged_editable.AddBond(
                                            targa, targb + buffer,
                                            order=Chem.rdchem.BondType.SINGLE,
                                        )
                                        combined_mol = merged_editable.GetMol()

                                        pba = combined_mol.GetAtomWithIdx(targa).GetBonds()
                                        pbb = combined_mol.GetAtomWithIdx(targb + buffer).GetBonds()
                                        _fix_overvalent_bonds(combined_mol, [pba, pbb])
                                        Chem.SanitizeMol(combined_mol)

                                        # collect terminal atoms for H assignment
                                        # (both vinyl and radical sides need terminal H assignment)
                                        search   = None
                                        receptor = None
                                        if isnumber(bna) or bna == "RADICAL":
                                            search = fmterminals[fma][targa]
                                            unreacted.append(search)
                                            unreacted_types.append("search")
                                        if isnumber(bnb) or bnb == "RADICAL":
                                            receptor = fmterminals[fmb][targb]
                                            unreacted.append(receptor + buffer)
                                            unreacted_types.append("receptor")

                                        editable = Chem.EditableMol(combined_mol)
                                        combined_mol = editable.GetMol()

                                        # search closes to a normal terminus (+1 H); receptor
                                        # stays one short, becoming the new radical site
                                        _assign_explicit_hs(
                                            combined_mol, unreacted, unreacted_types,
                                            species_count, typed=True,
                                        )

                                        combined_mol = Chem.rdmolops.AddHs(combined_mol, addCoords=True)
                                        Chem.SanitizeMol(combined_mol)

                                        combined_mol = _renumber_and_finalize_mol(combined_mol, species_count)

                                        reacted = mol_to_dictionary(combined_mol, f"X{species_count}", PT)
                                        reacted_mols.append(na)
                                        reacted_mols.append(nb)
                                        step_complex.append(reacted)

                                        w = Chem.PDBWriter(pargs.compdir + f"X{species_count}.pdb")
                                        w.write(combined_mol)
                                        w.close()

                                        _update_polymer_indexes(
                                            combined_mol, species_count,
                                            fmsmiles, fmterminals, polymers,
                                            search=search,    targa=targa,
                                            receptor=receptor, targb=targb,
                                            buffer=buffer,
                                        )

                                        if pargs.template_name:
                                            print_bondstats(
                                                pargs.workdir, statname, i, species_count,
                                                fma, positionsa, atomsa,
                                                fmb, positionsb, atomsb,
                                                template_coords, template_res, template_atom,
                                            )

                                        species_count += 1

                                        del (
                                            mola, molb, confa, confb,
                                            combined_mol, merged, merged_editable,
                                        )

                            # not reacted and same molecule
                            elif (
                                ta not in n_fm_set
                                and tb not in n_fm_set
                                and self_interaction
                            ):
                                unreacted = []
                                unreacted_types = []

                                hcount = []

                                # molecule A
                                if fma[0:1] == "X":
                                    mola = Chem.MolFromPDBFile(
                                        pargs.compdir+f"{fma}.pdb",
                                        sanitize=True,
                                        removeHs=False,
                                        proximityBonding=False,
                                    )
                                    embed_success = -1
                                    while embed_success < 0:
                                        embed_success = AllChem.EmbedMolecule(
                                            mola, params
                                        )
                                        if embed_success == -1:
                                            print("Embedding failure. Retrying...")

                                    # look for valence mismatch  -- polymerized locations will be missing one hydrogen
                                    for atom in mola.GetAtoms():
                                        if (
                                            atom.GetAtomicNum() == 6
                                            and atom.GetNumImplicitHs() == 1
                                            and atom.GetIdx() != targa
                                        ):
                                            unreacted.append(atom.GetIdx())
                                            unreacted_types.append("RADICAL")
                                            hc = 0
                                            for nbr in atom.GetNeighbors():
                                                if nbr.GetAtomicNum() == 1:
                                                    hc += 1
                                            hcount.append(hc)

                                    del mola
                                    mola = Chem.MolFromPDBFile(
                                        pargs.compdir+f"{fma}.pdb",
                                        sanitize=True,
                                        removeHs=True,
                                        proximityBonding=False,
                                    )
                                    if fma not in polymers:
                                        raise ValueError(
                                            f"Polymer {fma} not found in polymer registry. "
                                            f"Check if it was properly created in previous step."
                                        )
                                    matcha = polymers[fma]["rdkitmol"]
                                elif fma in sulfurs:
                                    mola = Chem.MolFromPDBFile(
                                        pargs.sulfurdir+f"/{fma}.pdb",
                                        sanitize=True,
                                        removeHs=True,
                                        proximityBonding=False,
                                    )
                                else:
                                    mola = Chem.MolFromSmiles(fmsmiles[fma])
                                embed_success = -1
                                while embed_success < 0:
                                    embed_success = AllChem.EmbedMolecule(mola, params)
                                    if embed_success == -1:
                                        print("Embedding failure. Retrying...")
                                mola = Chem.RemoveAllHs(mola)

                                no_bond = True
                                atoma = mola.GetAtomWithIdx(targa)
                                for bn, bond in enumerate(atoma.GetBonds()):
                                    a = bond.GetBeginAtomIdx()
                                    b = bond.GetEndAtomIdx()
                                    if a == targb or b == targb:
                                        no_bond = False

                                #
                                # Vinyl - vinyl bonding
                                #
                                #
                                #

                                # ---- Self-interaction: vinyl-vinyl bonding ----
                                if isnumber(bna) and isnumber(bnb) and no_bond:
                                    bt1 = mola.GetBondWithIdx(bna).GetBondType()
                                    bt2 = mola.GetBondWithIdx(bnb).GetBondType()

                                    if bt1 == BondType.DOUBLE and bt2 == BondType.DOUBLE:
                                        mola.GetBondWithIdx(bna).SetBondType(BondType.SINGLE)
                                        mola.GetBondWithIdx(bnb).SetBondType(BondType.SINGLE)

                                        n_fm.append(ta)

                                        bmbt.write(
                                            f"Bonding {fma} {targa} and {fmb} {targb}"
                                            f" ({targb}) as X{species_count}\n"
                                        )

                                        # set 3D coordinates from the pdb
                                        positionsa = precomplex[na]["coords"]
                                        atomsa     = precomplex[na]["atoms"]
                                        positionsb = precomplex[nb]["coords"]
                                        atomsb     = precomplex[nb]["atoms"]

                                        _verify_atom_correspondence(mola, positionsa, atomsa, context=f"self-interaction vinyl-vinyl, {fma}")

                                        confa = mola.GetConformer()

                                        if fma[0:1] == "X":
                                            for _m in matcha:
                                                x, y, z = positionsa[_m]
                                                confa.SetAtomPosition(_m, Point3D(x, y, z))
                                        else:
                                            for _m in range(mola.GetNumAtoms()):
                                                x, y, z = positionsa[_m]
                                                confa.SetAtomPosition(_m, Point3D(x, y, z))

                                        # single molecule – just add the intramolecular bond
                                        # https://gist.github.com/dominiquesydow/399899aed83c5583cfabdfd18b3cd7dd
                                        merged_editable = Chem.EditableMol(mola)
                                        merged_editable.AddBond(
                                            targa, targb,
                                            order=Chem.rdchem.BondType.SINGLE,
                                        )
                                        combined_mol = merged_editable.GetMol()
                                        Chem.SanitizeMol(combined_mol)

                                        # assign terminal atoms and hydrogen counts
                                        search   = fmterminals[fma][targa]
                                        receptor = fmterminals[fma][targb]
                                        unreacted.append(search)
                                        unreacted_types.append("search")
                                        unreacted.append(receptor)
                                        unreacted_types.append("receptor")

                                        editable = Chem.EditableMol(combined_mol)
                                        combined_mol = editable.GetMol()

                                        _assign_explicit_hs(
                                            combined_mol, unreacted, unreacted_types, species_count,
                                        )

                                        combined_mol = Chem.rdmolops.AddHs(combined_mol, addCoords=True)
                                        Chem.SanitizeMol(combined_mol)

                                        combined_mol = _renumber_and_finalize_mol(combined_mol, species_count)

                                        reacted = mol_to_dictionary(combined_mol, f"X{species_count}", PT)
                                        reacted_mols.append(na)
                                        reacted_mols.append(nb)
                                        step_complex.append(reacted)

                                        w = Chem.PDBWriter(pargs.compdir + f"X{species_count}.pdb")
                                        w.write(combined_mol)
                                        w.close()

                                        # buffer=0 for self-interactions (single molecule)
                                        _update_polymer_indexes(
                                            combined_mol, species_count,
                                            fmsmiles, fmterminals, polymers,
                                            search=search,    targa=targa,
                                            receptor=receptor, targb=targb,
                                            buffer=0,
                                        )

                                        if pargs.template_name:
                                            print_bondstats(
                                                pargs.workdir, statname, i, species_count,
                                                fma, positionsa, atomsa,
                                                fmb, positionsb, atomsb,
                                                template_coords, template_res, template_atom,
                                            )

                                        species_count += 1

                                        del mola, confa, combined_mol

                                    else:
                                        continue

                                #
                                # Bonding a vinyl functional group to a free radical carbon
                                #

                                # ---- Self-interaction: vinyl-radical bonding ----
                                elif (
                                    isnumber(bna) and bnb == "RADICAL" and no_bond
                                ) or (isnumber(bnb) and bna == "RADICAL" and no_bond):
                                    if isnumber(bna):
                                        bt = mola.GetBondWithIdx(bna).GetBondType()
                                    else:
                                        bt = mola.GetBondWithIdx(bnb).GetBondType()

                                    if bt == BondType.DOUBLE:
                                        if isnumber(bna):
                                            mola.GetBondWithIdx(bna).SetBondType(BondType.SINGLE)
                                        if isnumber(bnb):
                                            mola.GetBondWithIdx(bnb).SetBondType(BondType.SINGLE)

                                        n_fm.append(ta)

                                        bmbt.write(
                                            f"Bonding {fma} {targa} and {fmb} {targb}"
                                            f" ({targb}) as X{species_count}\n"
                                        )

                                        # set 3D coordinates from the pdb
                                        positionsa = precomplex[na]["coords"]
                                        atomsa     = precomplex[na]["atoms"]
                                        positionsb = precomplex[nb]["coords"]
                                        atomsb     = precomplex[nb]["atoms"]

                                        _verify_atom_correspondence(mola, positionsa, atomsa, context=f"self-interaction vinyl-radical, {fma}")

                                        confa = mola.GetConformer()

                                        if fma[0:1] == "X":
                                            for _m in matcha:
                                                x, y, z = positionsa[_m]
                                                confa.SetAtomPosition(_m, Point3D(x, y, z))
                                        else:
                                            for _m in range(mola.GetNumAtoms()):
                                                x, y, z = positionsa[_m]
                                                confa.SetAtomPosition(_m, Point3D(x, y, z))

                                        # single molecule – intramolecular bond
                                        # https://gist.github.com/dominiquesydow/399899aed83c5583cfabdfd18b3cd7dd
                                        merged_editable = Chem.EditableMol(mola)
                                        merged_editable.AddBond(
                                            targa, targb,
                                            order=Chem.rdchem.BondType.SINGLE,
                                        )
                                        combined_mol = merged_editable.GetMol()
                                        Chem.SanitizeMol(combined_mol)

                                        # collect terminal atoms (only those with double bonds)
                                        search   = None
                                        receptor = None
                                        if isnumber(bna):
                                            search = fmterminals[fma][targa]
                                            unreacted.append(search)
                                            unreacted_types.append("search")
                                        if isnumber(bnb):
                                            receptor = fmterminals[fma][targb]
                                            unreacted.append(receptor)
                                            unreacted_types.append("receptor")

                                        editable = Chem.EditableMol(combined_mol)
                                        combined_mol = editable.GetMol()

                                        _assign_explicit_hs(
                                            combined_mol, unreacted, unreacted_types, species_count,
                                        )

                                        combined_mol = Chem.rdmolops.AddHs(combined_mol, addCoords=True)
                                        Chem.SanitizeMol(combined_mol)

                                        combined_mol = _renumber_and_finalize_mol(combined_mol, species_count)

                                        reacted = mol_to_dictionary(combined_mol, f"X{species_count}", PT)
                                        reacted_mols.append(na)
                                        reacted_mols.append(nb)
                                        step_complex.append(reacted)

                                        w = Chem.PDBWriter(pargs.compdir + f"X{species_count}.pdb")
                                        w.write(combined_mol)
                                        w.close()

                                        # buffer=0 for self-interactions (single molecule)
                                        _update_polymer_indexes(
                                            combined_mol, species_count,
                                            fmsmiles, fmterminals, polymers,
                                            search=search,    targa=targa,
                                            receptor=receptor, targb=targb,
                                            buffer=0,
                                        )

                                        if pargs.template_name:
                                            print_bondstats(
                                                pargs.workdir, statname, i, species_count,
                                                fma, positionsa, atomsa,
                                                fmb, positionsb, atomsb,
                                                template_coords, template_res, template_atom,
                                            )

                                        species_count += 1

                                        del mola, confa, combined_mol

                ################################################################
                #
                # Process non-reacted species
                #
                ################################################################

                for j, fm in enumerate(precomplex):
                    if j in set(reacted_mols):
                        continue

                    else:
                        """
                        Determine polymers from pdb, determine fms from smiles
                        --
                        some of the polymer smiles were coming out weird which is why we are not using them

                        """
                        receptor_atoms = []
                        hcount = []

                        if fm["name"][0:1] == "X":
                            _fm = fm["name"]
                            cached_mol = polymers[_fm]["mol"]

                            # look for valence mismatch  -- polymerized locations will be missing one hydrogen
                            for atom in cached_mol.GetAtoms():
                                if (
                                    atom.GetAtomicNum() == 6
                                    and atom.GetNumImplicitHs() == 1
                                ):
                                    receptor_atoms.append(atom.GetIdx())
                                    hc = 0
                                    for nbr in atom.GetNeighbors():
                                        if nbr.GetAtomicNum() == 1:
                                            hc += 1
                                    hcount.append(hc)

                            # Stage 1: work from a copy of the cached in-memory Mol
                            # instead of re-reading X{n}.pdb from disk (same atom
                            # order as the disk file, since the cache holds the
                            # exact heavy-atom-only mol that was derived from it).
                            # A copy, not the cached object itself, since
                            # AllChem.EmbedMolecule below mutates in place and we
                            # don't want conformers accumulating on the cache.
                            mol = Chem.Mol(cached_mol)
                            match = polymers[_fm]["rdkitmol"]

                        else:
                            _fm = fm["name"]
                            mol = Chem.MolFromSmiles(fmsmiles[defm[_fm]])
                            mol = Chem.RemoveAllHs(mol)

                        embed_success = -1
                        while embed_success < 0:
                            embed_success = AllChem.EmbedMolecule(mol, params)
                            if embed_success == -1:
                                print("Embedding failure. Retrying...")
                        Chem.SanitizeMol(mol)

                        positions = precomplex[j]["coords"]

                        _verify_atom_correspondence(
                            mol, positions, precomplex[j].get("atoms"), context=f"non-reacted species, {fm['name']}"
                        )

                        conf = mol.GetConformer()

                        """
                        Set positions to those used in the pdb file

                        Then add hydrogens, get locations, bonds, and types to write out

                        """

                        if fm["name"][0:1] == "X":
                            for _m in match:
                                x, y, z = positions[_m]
                                conf.SetAtomPosition(_m, Point3D(x, y, z))

                        else:
                            for _m in range(mol.GetNumAtoms()):
                                x, y, z = positions[_m]
                                conf.SetAtomPosition(_m, Point3D(x, y, z))

                        #
                        # Remove hydrogens from unreacted atoms.
                        #

                        for _i, idx in enumerate(receptor_atoms):
                            atom = mol.GetAtomWithIdx(idx)
                            atom.SetNoImplicit(True)
                            atom.SetNumExplicitHs(hcount[_i])

                        mol = Chem.rdmolops.AddHs(mol, addCoords=True)

                        Chem.SanitizeMol(mol)

                        new_fm = mol_to_dictionary(mol, fm["name"], PT)
                        step_complex.append(new_fm)

                """
                Write out labels used at this stage
                """

                labels = []
                for p in step_complex:
                    labels.append(p["name"])

                pargs.fms = labels

                bmbt.write(f"\n{labels}\n")

                ################################################################
                #
                # Openbabel Run
                #
                ################################################################
                logname = f"obabel-step{i}.log"
                out_complex = obabel_do("rxn",i, step_complex, pargs, logname)

                ################################################################
                #
                # acpype / antechamber Run
                #
                ################################################################
                logname = f"acpype-step{i}.log"
                acpype_do("rxn", i, step_complex, out_complex, pargs, logname)

                ################################################################
                #
                # Combine Template and Complex
                #
                ################################################################
                logname = f"combination-step{i}.log"
                obabel_combine("rxn",i, out_complex, pargs, logname)

                if args.explicit and initiator_molecules is not None and len(initiator_molecules) > 0:
                    obabel_initiators(i, initiator_molecules, pargs, logname)

                # obabel_solvent("rxn",i, pargs, logname)

                ################################################################
                #
                # GROMACS Run
                #
                ################################################################
                logname = f"gmx-step{i}.log"

                gmx_run = gmx_do()
                gmx_run.params('rxn', i, pargs.ncycles, pargs, logname)
                incomplex, intemplate, ininitiator = gmx_run.run()

                pargs.intemplate = intemplate
                pargs.incomplex = incomplex
                pargs.ininitiator = ininitiator

                pb.update(pb1, advance=1)

        if pargs.clean and pargs.template_name:
            # keep template_posres.itp: it's only (re)generated on a run's
            # very first cycle (see run_acpype.py's "just do once" gate), so
            # deleting it here would break any future -restart of this run
            workdir_q = shlex.quote(pargs.workdir)
            subprocess.run(
                f"find {workdir_q} -maxdepth 1 -name 'template*' ! -name 'template_posres.itp' -delete",
                shell=True,
            )

        ################################################################
        #
        # Stats for final step only
        #
        ################################################################
        nmono = 0
        npoly = 0
        mws = []
        for j, fms in enumerate(precomplex):
            if fms["name"][0:1] == "X":
                fm = fms["name"]
                mol = Chem.MolFromPDBFile(
                    pargs.compdir+f"{fm}.pdb",
                    sanitize=True,
                    removeHs=True,
                    proximityBonding=False,
                )
                npoly += 1

            else:
                fm = defm[fms["name"]]
                mol = Chem.MolFromSmiles(fmsmiles[fm])
                nmono += 1

            AllChem.EmbedMolecule(mol, params)
            Chem.SanitizeMol(mol)

            # get weights
            mws.append(Chem.Descriptors.ExactMolWt(mol))
            del mol

        print_stats(pargs.workdir, i + 1, nmono, npoly, mws)
        print_mw_distributions(pargs.workdir, i + 1, totmon, mws)
