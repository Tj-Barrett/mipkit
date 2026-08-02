# Wrapper for gnina
import shlex
import subprocess, os
from MIPkit.utils.read_pdb import read_pdb
from MIPkit.utils.print_molecule import print_molecule
from MIPkit.constants.constants import fm2smiles, sulfur_exception
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

from networkx.algorithms import isomorphism
import networkx as nx

def gnina(log, step, pargs, template, fm, vina_args):
    logdir = pargs.logdir
    docked = f"ligand-poses{step}.pdb"
    docked_q = shlex.quote(docked)

    mode = "w+" if step == 0 else "a"
    with open(os.path.join(logdir, f"{log}"), mode) as f:
        gnina_cms = ''
        gnina_cms+= f"--receptor {shlex.quote(template)} "
        gnina_cms+= f"--ligand {shlex.quote(fm)} "

        ptx, pty, ptz = vina_args['center']
        gnina_cms+= f"--center_x {ptx} --center_y {pty} --center_z {ptz} "

        ptx, pty, ptz = vina_args['size']
        gnina_cms+= f" --size_x {ptx} --size_y {pty} --size_z {ptz}"

        gnina_cms+= f" --exhaustiveness {pargs.exhaustiveness} --cpu {pargs.ncpu} "


        gnina_cms+= f" -o tmpgnina/{docked_q} " #--minimize_single_full, --scoring = vina

        output = subprocess.run(
                f'{shlex.quote(pargs.gnina)} {gnina_cms}',
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
        if output.returncode != 0:
            raise RuntimeError(f"gnina failed (exit code {output.returncode}): {output.stderr}")
        output = output.stdout.split('\n')
        get_first = False
        stats = {'affinity':0,
                 'intramol':0,
                 'posescore':0,
                 'cnnaffinity':0,}

        for line in output:
            if get_first:
                s = line.split()
                stats['affinity'] = float(s[1])
                stats['intramol'] = float(s[2])
                stats['posescore'] = float(s[3])
                stats['cnnaffinity'] = float(s[4])
                get_first = False
            if line == "-----+------------+------------+------------+----------":
                get_first = True
            f.write(line+'\n')

    mols = read_pdb(f"tmpgnina/{docked}", conect=True, decode=False, hydrogens=True)

    if not mols:
        raise ValueError("gnina failed to produce valid molecular poses")
    return mols[-1], stats
    
def update_gnina(pargs, iter, fm, ligandmol, basename):
    RDLogger.DisableLog("rdApp.*")
    fm2sm = fm2smiles()
    sulfurs = sulfur_exception()

    if iter == 0:
        targetpdb = basename
    else:
        targetpdb = os.path.splitext(os.path.basename(basename))[0] + f"{iter-1}.pdb"
        targetpdb = "tmpgnina/" + targetpdb

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
                    "tmpgnina/tmp.pdb",
                    conect=True)

    ligand = Chem.MolFromPDBFile(
        "tmpgnina/tmp.pdb",
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
        raise ValueError('GNINA and Reference graphs do not match')

    for i in range(n_atoms):
        x,y,z = positions[g_matching.mapping[i]]
        conf.SetAtomPosition(i,Point3D(x,y,z))

    w = Chem.PDBWriter(
        "tmpgnina/tmp.pdb"
    )
    w.write(target)
    w.close()

    mols = read_pdb(targetpdb, hydrogens=True)
    ligandmol = read_pdb("tmpgnina/tmp.pdb", hydrogens=True)

    if not ligandmol:
        raise ValueError("Failed to read ligand molecule from tmpgnina/tmp.pdb")
    ligandmol = ligandmol[0]
    ligandmol['name']=fm
    mols.append(ligandmol)

    # New template 
    templatepdb = os.path.splitext(os.path.basename(basename))[0] + f"{iter}.pdb"
    templatepdb = "tmpgnina/" + templatepdb

    print_molecule(mols, templatepdb,  amino = 'pdb')

    return templatepdb

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