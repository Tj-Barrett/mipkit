import os
import shlex
import subprocess
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

from networkx.algorithms import isomorphism
import networkx as nx

from MIPkit.utils.read_pdb import read_pdb
from MIPkit.utils.print_molecule import print_molecule
from MIPkit.constants.constants import fm2smiles, sulfur_exception

def autodock_vina(log, step, pargs, template, fm, vina_args):
    logdir = pargs.logdir
    docked = f"ligand-poses{step}.pdbqt"
    dockedpdb = f"ligand-poses{step}.pdb"
    docked_q = shlex.quote(docked)
    dockedpdb_q = shlex.quote(dockedpdb)

    mode = "w+" if step == 0 else "a"
    with open(os.path.join(logdir, log), mode) as f:
        vina_cms = ''
        vina_cms+= f"--receptor {shlex.quote(template)} "
        vina_cms+= f"--ligand {shlex.quote(fm)} "

        ptx, pty, ptz = vina_args['center']
        vina_cms+= f"--center_x {ptx} --center_y {pty} --center_z {ptz} "

        ptx, pty, ptz = vina_args['size']
        vina_cms+= f" --size_x {ptx} --size_y {pty} --size_z {ptz}"
        vina_cms+= f" --exhaustiveness {pargs.exhaustiveness} --energy_range {pargs.energy_range} --cpu {pargs.ncpu} "

        vina_cms+= f" --out tmpvina/{docked_q} "

        output = subprocess.run(
                f'{shlex.quote(pargs.vina)} {vina_cms}',
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
        if output.returncode != 0:
            raise RuntimeError(f"vina failed (exit code {output.returncode}): {output.stderr}")
        output = output.stdout.split('\n')

        get_first = False
        stats = []

        for line in output:
            if len(line) > 0:
                if line == "Writing output ... done.":
                    get_first = False

                if get_first:
                    s = line.split()
                    conformer = {   'vina_rmsd_lb':float(s[2]),
                                    'vina_rmsd_ub':float(s[3]),
                                    'vina_affinity':float(s[1]),}
                    stats.append(conformer)
                if line == "-----+------------+----------+----------":
                    get_first = True
                f.write(line+'\n')

    convert = subprocess.run(
        f"{shlex.quote(pargs.obabel)} -ipdbqt tmpvina/{docked_q} -opdb -O tmpvina/{dockedpdb_q}",
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
    )
    if convert.returncode == 0:
        subprocess.run(
            f"rm tmpvina/{docked_q}",
            shell=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

    mols = read_pdb(f"tmpvina/{dockedpdb}", conect=True, decode=False, hydrogens=True)

    if not mols:
        raise ValueError("vina failed to produce valid molecular poses")

    return mols, stats

def preprocessing(pargs, ligand=False, receptor=False):

    obabel = pargs.obabel
    obabel_q = shlex.quote(obabel)

    if not os.path.isdir("tmpvina"):
        os.mkdir("tmpvina")

    if ligand:
        _ligand = os.path.splitext(os.path.basename(ligand))[0] + ".pdbqt"

        print("Converting ligand")
        if os.path.splitext(ligand)[-1] == ".sdf":
            subprocess.run(
                f"{obabel_q} -isdf {shlex.quote(ligand)} -opdbqt -O tmpvina/{shlex.quote(_ligand)} -xr",
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )
        elif os.path.splitext(ligand)[-1] == ".pdb":
            subprocess.run(
                f"{obabel_q} -ipdb {shlex.quote(ligand)} -opdbqt -O tmpvina/{shlex.quote(_ligand)} -xr",
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )
        else:
            raise ValueError(f"Unsupported ligand file type: {ligand}. Expected .sdf or .pdb.")

    if receptor:
        _receptor = os.path.splitext(os.path.basename(receptor))[0] + ".pdbqt"
        receptor_q = shlex.quote(receptor)
        _receptor_q = shlex.quote(_receptor)

        print("Converting receptor")
        if os.path.splitext(receptor)[-1] == ".sdf":
            subprocess.run(
                f"{obabel_q}  -isdf {receptor_q} -opdbqt -O tmpvina/{_receptor_q} -xr",
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )
        elif os.path.splitext(receptor)[-1] == ".pdb":
            subprocess.run(
                f"{obabel_q}  -ipdb {receptor_q} -opdbqt -O tmpvina/{_receptor_q} -xr",
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )
        else:
            raise ValueError(f"Unsupported receptor file type: {receptor}. Expected .sdf or .pdb.")

    if ligand and receptor:
        return f"tmpvina/{_ligand}", f"tmpvina/{_receptor}"
    elif ligand:
        return f"tmpvina/{_ligand}"
    elif receptor:
        return f"tmpvina/{_receptor}"


def update_vina(pargs, iter, fm, ligandmol, basename):
    RDLogger.DisableLog("rdApp.*")
    fm2sm = fm2smiles()
    sulfurs = sulfur_exception()

    if iter == 0:
        targetpdbqt = basename

        targetpdb = os.path.splitext(os.path.basename(basename))[0] + f".pdb"
        targetpdb = "tmpvina/" + targetpdb

    else:
        targetpdbqt = os.path.splitext(os.path.basename(basename))[0] + f"{iter-1}.pdbqt"
        targetpdbqt = "tmpvina/" + targetpdbqt

        targetpdb = os.path.splitext(os.path.basename(basename))[0] + f"{iter-1}.pdb"
        targetpdb = "tmpvina/" + targetpdb

    subprocess.run(
        f"{shlex.quote(pargs.obabel)} -ipdbqt {shlex.quote(targetpdbqt)} -opdb -O {shlex.quote(targetpdb)}",
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
    )

    '''
    We need to reorder the molecule from gnina.
    gnina uses tree function to determine stiff and flexible parts which reorders them

    https://ljmartin.github.io/sideprojects/pdbqt_to_mol.html

    See LJ Martin's post for the networkx graph comparison to map locations efficiently

    '''
    #
    # We will take the ligands position and map it to our known molecule
    #
    ligandmol['name'] = fm
    print_molecule([ligandmol],
                    "tmpvina/tmp.pdb",
                    conect=True)

    ligand = Chem.MolFromPDBFile(
        "tmpvina/tmp.pdb",
        sanitize=True,
        removeHs=False,
        proximityBonding=False,
    )
    ligand = Chem.RemoveAllHs(ligand)
    positions = ligand.GetConformer().GetPositions()
    n_atoms   = ligand.GetNumAtoms()

    #
    # This will be our target for mapping
    #
    target = Chem.MolFromSmiles(fm2sm[fm])
    target = Chem.RemoveAllHs(target)
    AllChem.EmbedMolecule(target)
    conf = target.GetConformer()

    # Use a graph of node points with bond edges
    g_matching = isomorphism.GraphMatcher(
        mol_to_nx(target),
        mol_to_nx(ligand),
        node_match=isomorphism.categorical_node_match('atomic_num', -1),
    )

    if not g_matching.is_isomorphic():
        raise ValueError('VINA and Reference graphs do not match')

    for i in range(n_atoms):
        x,y,z = positions[g_matching.mapping[i]]
        conf.SetAtomPosition(i,Point3D(x,y,z))

    w = Chem.PDBWriter(
        "tmpvina/tmp.pdb"
    )
    w.write(target)
    w.close()

    mols = read_pdb(targetpdb, hydrogens=True)
    ligandmol = read_pdb("tmpvina/tmp.pdb", hydrogens=True)

    ligandmol = ligandmol[0]
    ligandmol['name']=fm
    mols.append(ligandmol)

    # Print template for next step
    templatepdb = os.path.splitext(os.path.basename(basename))[0] + f"{iter}.pdb"
    templatepdb = "tmpvina/" + templatepdb

    templatepdbqt = os.path.splitext(os.path.basename(basename))[0] + f"{iter}.pdbqt"
    templatepdbqt = "tmpvina/" + templatepdbqt

    print_molecule(mols, templatepdb,  amino = 'pdb')
    
    subprocess.run(
        f"{shlex.quote(pargs.obabel)} -ipdb {shlex.quote(templatepdb)} -opdbqt -O {shlex.quote(templatepdbqt)} -xr",
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
    )

    return templatepdb, templatepdbqt

def mol_to_nx(mol):
    #
    # See the implementation here:
    #
    # https://ljmartin.github.io/sideprojects/pdbqt_to_mol.html
    #
    G = nx.Graph()

    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(),
                   atomic_num=atom.GetAtomicNum())
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx())
    return G