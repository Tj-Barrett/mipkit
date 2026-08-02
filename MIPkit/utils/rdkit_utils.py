import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem, RemoveHs

from MIPkit.aminos.handle_aminos import import_amino_lib
from MIPkit.constants.constants import encode_fms, fm2smiles

def generate_amino_sdfs(hydrogens=True, writesdf=True):
    sdfdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../sdfs/"))
    # https://gist.github.com/leelasd/43219a222bf57d3e01c2c83f2ad9b031
    rdprep = {}
    aminos = import_amino_lib()

    for am in aminos:
        mol = Chem.MolFromSmiles(aminos[am])
        mol = Chem.AddHs(mol)

        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol, 1000)

        if not hydrogens:
            mol = RemoveHs(mol)

        if writesdf:
            w = Chem.SDWriter(f"{sdfdir}/{am}.sdf")
            w.write(mol)
            w.close()

    return  # rdprep

def generate_sdfs(hydrogens=True, writesdf=True):
    sdfdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../sdfs/"))
    # https://gist.github.com/leelasd/43219a222bf57d3e01c2c83f2ad9b031
    rdprep = {}
    fms = fm2smiles()
    for fm in fms:
        mol = Chem.MolFromSmiles(fms[fm])
        mol = Chem.AddHs(mol)

        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol, 1000)

        if not hydrogens:
            mol = RemoveHs(mol)

        # rdprep[fm]=mol

        if writesdf:
            w = Chem.SDWriter(f"{sdfdir}/{fm}.sdf")
            w.write(mol)
            w.close()

    return  # rdprep

def write_fm_pdb(arg):
    fms = fm2smiles()
    enfm = encode_fms()
    try:
        fm = arg.write_pdb[0]
        mol = Chem.MolFromSmiles(fms[fm])
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol, 1000)

        # https://sourceforge.net/p/rdkit/mailman/message/36404394/
        mi = Chem.AtomPDBResidueInfo()
        mi.SetResidueName(enfm[fm])

        [a.SetMonomerInfo(mi) for a in mol.GetAtoms()]

        w = Chem.PDBWriter(f"{fm}.pdb")
        w.write(mol)
        w.close()

        print(f"Wrote {fm}.pdb. {fm} is encoded as {enfm[fm]}.")

    except (KeyError, ValueError, IOError, RuntimeError) as e:
        print(f"Not a recognized FM from the FM library or RDKit operation failed: {e}")
        sys.exit(1)