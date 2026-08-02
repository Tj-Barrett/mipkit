import re, subprocess
from rdkit import rdBase
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, RemoveHs
from MIPkit.constants.constants import (
    fm2smiles,
    encode_fms,
)
from MIPkit.utils.read_pdb import read_pdb
from MIPkit.utils.print_molecule import print_molecule

def str_find(fm, mol):
    if fm in mol:
        return True
    else:
        return False

def id_mono_det(mol, ref_smiles):
    len_ref = []
    for ref in ref_smiles:
        len_ref.append(len(ref))

    molstr = ""
    mol_atoms = []
    for atm in mol['atoms']:
        if atm != "H":
            molstr+=atm
            mol_atoms.append(atm)

    if molstr in ref_smiles:
        return [ref_smiles.index(molstr)], [[0, len(molstr)]], mol_atoms
    else:
        ms = 0
        outstr = ""

        reacted_fms = []
        indexes = []

        while molstr != outstr:
            matched_this_pass = False
            for i, ref in enumerate(ref_smiles):
                p_str = molstr[ms:ms+len(ref)]

                if str_find( ref , p_str):
                    _ms = ms
                    ms += len(ref)
                    outstr+=ref

                    indexes.append([_ms, ms])

                    reacted_fms.append(i)
                    matched_this_pass = True

            if not matched_this_pass:
                print(f"Warning: could not fully decompose molecule into known FM fragments (stopped at offset {ms}/{len(molstr)}).")
                break

        return reacted_fms, indexes, mol_atoms

def identify(pargs, recipe, part, outname=None):
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

    workdir = pargs.workdir

    detected_fms={}
    for fm in recipe:
        detected_fms[fm]=0

    #
    # unreacted = monomers
    #
    # polymers = polymers
    #

    unreacted = []
    polymers = []

    partmol = read_pdb(part, hydrogens=True)

    for mol in partmol:
        if mol['name'][0]=='X':
            polymers.append(mol)
        else:
            unreacted.append(mol)

    fm_smiles = fm2smiles()
    smiles = {}
    for fm in recipe:
        smiles[fm]=fm_smiles[fm]

    rdkit_mols = []
    ref_smiles = []

    for fm in smiles:
        _mol = 0

        mol = Chem.MolFromSmiles(smiles[fm])

        for at in mol.GetAtoms():
            if at.GetAtomicNum() != 1:
                _mol += 1

        rdkit_mols.append(_mol)

        a = re.sub(r'[^a-zA-Z]', '', smiles[fm])
        ref_smiles.append(a)


    all_ufm = []
    for mol in unreacted:

        ufm = id_mono_det(mol, ref_smiles)
        all_ufm.append(ufm)

    all_fms = []
    all_indexes = []
    all_atoms = []
    for mol in polymers:

        fms, indexes, mol_atoms = id_mono_det(mol, ref_smiles)
        all_fms.append(fms)
        all_indexes.append(indexes)
        all_atoms.append(mol_atoms)

    names = []
    for n in detected_fms:
        names.append(n)

    enfm = encode_fms()
    for i, mol in enumerate(polymers):

        mol = Chem.MolFromPDBFile(
            part,
            sanitize=True,
            removeHs=False,
            proximityBonding=True,
        )

        # all_indexes[i]'s start/end are sequential positions among HEAVY atoms
        # only (id_mono_det skips "H" entries) -- map those back to this mol's
        # actual atom indices, which include interspersed hydrogens.
        heavy_idx_map = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1]

        for j, se in enumerate(all_indexes[i]):
            nm = enfm[names[all_fms[i][j]]]

            start = se[0]
            end = se[1]

            mi = Chem.AtomPDBResidueInfo()
            mi.SetResidueName(nm)

            target_idxs = set(heavy_idx_map[start:end])
            for a in mol.GetAtoms():
                if a.GetIdx() in target_idxs:
                    a.SetMonomerInfo(mi)

        # make hydrogens inherit heavy atom residue name
        for a in mol.GetAtoms():
            if a.GetAtomicNum() == 1:
                neigh = a.GetNeighbors()
                heavy = neigh[0]
                hm = heavy.GetMonomerInfo()
                mi.SetResidueName(hm.GetResidueName())
                a.SetMonomerInfo(mi)

        # mol = Chem.AddHs(mol)
        w = Chem.PDBWriter(f"{part[0:-4]}-id.pdb")
        w.write(mol)
        w.close()

    # out = []
    # out.extend(new_mols)
    # out.extend(unreacted)
    #
    # if outname is not None:
    #     print_molecule(out,
    #         f'{workdir}{outname}',
    #         conect=True,
    #         encode=False,)
    #
    # return out