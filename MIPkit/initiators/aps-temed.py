
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, GetPeriodicTable, Descriptors, Draw
from rdkit.Geometry import Point3D
from rdkit.Chem.rdchem import BondType

aps = "[NH4+].O=S(=O)([O-])OOS([O-])(=O)=O.[NH4+]"
temed = "CN(C)CCN(C)C"

maps = Chem.MolFromSmiles(aps)
mtemed = Chem.MolFromSmiles(temed)

params = AllChem.ETKDGv3()
params.useRandomCoords = True
params.maxAttempts = 10000
params.enforceChirality = False

nn = ["aps", "temed"]

for i, mmol in enumerate([maps, mtemed]):
	mol = Chem.rdmolops.AddHs(mmol, addCoords=True)
	Chem.SanitizeMol(mol)
	AllChem.EmbedMolecule(mol, params)
	
	w = Chem.PDBWriter(
		f"{nn[i]}.pdb"
	)
	w.write(mol)
	w.close()


aps_fragment = "O=S(=O)([O])[O]"
maps = Chem.MolFromSmiles(aps_fragment)
# mol = Chem.rdmolops.AddHs(maps, addCoords=True)

d2d = Draw.MolDraw2DSVG(300,300)
d2d.DrawMolecule(maps)
d2d.FinishDrawing()
svg = d2d.GetDrawingText()
with open('APS.svg','w+') as outf:
    outf.write(svg)