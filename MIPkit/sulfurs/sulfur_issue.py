from MIPkit.utils.print_molecule import print_molecule

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, GetPeriodicTable, Draw
from rdkit.Geometry import Point3D
from rdkit.Chem.rdchem import BondType

PT = GetPeriodicTable()


def mol_with_atom_and_molecule_index(mol):
    """
    Taken from: https://iwatobipen.wordpress.com/2017/02/25/draw-molecule-with-atom-index-in-rdkit/

    Add label to each molecule's atom: "atom_name:atom_index"
    """

    for idx in range(mol.GetNumAtoms()):
        mol.GetAtomWithIdx(idx).SetProp(
            "molAtomMapNumber", str(mol.GetAtomWithIdx(idx).GetIdx())
        )

    return mol


smiles = "CC(C)(CS(=O)(=O)O)NC(=O)C=C"

mol = Chem.MolFromSmiles(smiles)

mol_with_atom_and_molecule_index(mol)

Draw.MolsToGridImage([mol])
# Draw.ShowMol(mol)

# mol = Chem.RemoveAllHs(mol)
# mol = Chem.rdmolops.AddHs(mol, addCoords=True)
# Chem.SanitizeMol(mol)
# mol_with_atom_and_molecule_index(mol)

# Draw.MolsToGridImage([mol])
# Draw.ShowMol(mol)


mol = Chem.rdmolops.AddHs(mol, addCoords=True)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
Chem.SanitizeMol(mol)

w = Chem.PDBWriter(f"AMPSA.pdb")
w.write(mol)
w.close()

# atoms = []
# coords = []
# bonds = {}
# for _m, atom in enumerate(mol.GetAtoms()):
#     ipos = mol.GetConformer().GetAtomPosition(_m)
#     numa = atom.GetAtomicNum()
#     coords.append( (ipos.x, ipos.y, ipos.z) )
#     atoms.append(PT.GetElementSymbol(numa))
#     bonds[_m]=[]


# for _b, bond in enumerate(mol.GetBonds()):
#     a = bond.GetBeginAtomIdx()
#     b = bond.GetEndAtomIdx()
#     bonds[a].append(b)

# # for b in bonds:
# #     if len(bonds[b])>0:
# #         print(f"{b} {bonds[b]}")

# new_fm = {'name':'',
#         'number':'',
#         'atoms':[],
#         'coords':[],
#         'bonds':[]}

# new_fm['name']='tt'
# new_fm['number']=1
# new_fm['atoms']=atoms
# new_fm['coords']=coords
# new_fm['bonds']=bonds

# all_fms = []
# for i in range(0,3):
#     fm = {'name':'',
#         'number':'',
#         'atoms':[],
#         'coords':[],
#         'bonds':[]}
#     crds = []
#     for j in coords:
#         a = j[0]
#         b = j[1]+10.0*float(i)
#         c = j[2]
#         crds.append((a,b,c))

#     fm['name']='tt'
#     fm['number']=i+1
#     fm['atoms']=atoms
#     fm['coords']=crds
#     fm['bonds']=bonds
#     all_fms.append(fm)

# print_molecule(all_fms, "work/test.pdb" ,conect=True, encode=False)
