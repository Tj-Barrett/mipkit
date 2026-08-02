
from functools import reduce

from MIPkit.utils.read_pdb import read_pdb
from MIPkit.initiators.initiators_yaml import init_attackers

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

params = AllChem.ETKDGv3()
params.useRandomCoords = True
params.enforceChirality = False


fm = 'MAA'
initiators = ['abn', 'azo', 'aps', ]

#
#
# Implicit
#
#
mola = Chem.MolFromPDBFile(
        f"pdbs/{fm}.pdb",
        sanitize=True,
        removeHs=False,
        proximityBonding=False,
    )

molb = Chem.MolFromPDBFile(
        f"pdbs/{fm}.pdb",
        sanitize=True,
        removeHs=False,
        proximityBonding=False,
    )

AllChem.EmbedMolecule(mola, params)
AllChem.Compute2DCoords(mola)

d = rdMolDraw2D.MolDraw2DCairo(250, 200) # or MolDraw2DSVG to get SVGs
d.drawOptions().addStereoAnnotation = True
d.drawOptions().addAtomIndices = True
d.drawOptions().addBondIndices = True
d.DrawMolecule(mola)
d.FinishDrawing()
d.WriteDrawingText('outputa.png')   


mola.GetBondWithIdx(4).SetBondType(
    Chem.rdchem.BondType.SINGLE
)
molb.GetBondWithIdx(4).SetBondType(
    Chem.rdchem.BondType.SINGLE
)

AllChem.EmbedMolecule(mola, params)
AllChem.Compute2DCoords(mola)

d = rdMolDraw2D.MolDraw2DCairo(250, 200) # or MolDraw2DSVG to get SVGs
d.drawOptions().addStereoAnnotation = True
d.drawOptions().addAtomIndices = True
d.drawOptions().addBondIndices = True
d.DrawMolecule(mola)
d.FinishDrawing()
d.WriteDrawingText('outputa_1.png')   

mola.GetAtomWithIdx(1).SetNoImplicit(True)
mola.GetAtomWithIdx(1).SetNumExplicitHs(1)
mola = Chem.AddHs(mola)
AllChem.EmbedMolecule(mola, params)
AllChem.Compute2DCoords(mola)

d = rdMolDraw2D.MolDraw2DCairo(250, 200) # or MolDraw2DSVG to get SVGs
d.drawOptions().addStereoAnnotation = True
d.drawOptions().addAtomIndices = True
d.drawOptions().addBondIndices = True
d.DrawMolecule(mola)
d.FinishDrawing()
d.WriteDrawingText('outputa_2.png')   

#
#
# Explicit
#
#

buffer = molb.GetNumAtoms()

for init in initiators:

    molc = Chem.MolFromPDBFile(
            f"pdbs/{init}.pdb",
            sanitize=True,
            removeHs=True,
            proximityBonding=False,
        )
    
    aps_temed = []
    init_oxy  = []

    indexes = []

    aps_temed = read_pdb(f"pdbs/{init}.pdb", hydrogens=True, initiators=True)

    print(init)

    init_attack = init_attackers()
    for j, inits in enumerate(aps_temed):
        atoms, count = init_attack.load(inits['name'])
        if inits['name'] == 'APS':
            indexes.extend(count)
        else:
            indexes.append(count[0])
        for _j, atm in enumerate(inits['atoms']):
            if atm in atoms and _j in count:
                init_oxy.append(inits['coords'][_j])
    # print(init_oxy)

    if init == 'aps':
        molc.GetAtomWithIdx(4).SetNoImplicit(True)
        molc.GetAtomWithIdx(4).SetNumExplicitHs(0)

    AllChem.EmbedMolecule(molc, params)
    AllChem.Compute2DCoords(molc)

    d = rdMolDraw2D.MolDraw2DCairo(250, 200) # or MolDraw2DSVG to get SVGs
    d.drawOptions().addStereoAnnotation = True
    d.drawOptions().addAtomIndices = True
    d.drawOptions().addBondIndices = True
    d.DrawMolecule(molc)
    d.FinishDrawing()
    d.WriteDrawingText(f'{init}.png')   

    for site in indexes:
        merged = reduce(Chem.CombineMols, [molb, molc])

        merged_editable = Chem.EditableMol(merged)

        merged_editable.AddBond(
            1,
            site + buffer,
            order=Chem.rdchem.BondType.SINGLE,
        )

        combined_mol = merged_editable.GetMol()
        Chem.SanitizeMol(combined_mol)

        combined_mol.GetAtomWithIdx(2).SetNoImplicit(True)
        combined_mol.GetAtomWithIdx(2).SetNumExplicitHs(0)
        combined_mol = Chem.AddHs(combined_mol)

        AllChem.EmbedMolecule(combined_mol, params)
        AllChem.Compute2DCoords(combined_mol)

        d = rdMolDraw2D.MolDraw2DCairo(250, 200) # or MolDraw2DSVG to get SVGs
        d.drawOptions().addStereoAnnotation = True
        d.drawOptions().addAtomIndices = True
        d.drawOptions().addBondIndices = True
        d.DrawMolecule(combined_mol)
        d.FinishDrawing()
        d.WriteDrawingText(f'outputb_{init}_{site}.png')