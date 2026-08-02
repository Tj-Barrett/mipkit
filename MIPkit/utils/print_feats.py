import os, subprocess
from MIPkit.constants.constants import (
	fm2smiles,)

import yaml

from copy import deepcopy

from rdkit import RDConfig	
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.FeatMaps import FeatMaps
from rdkit.Chem.Lipinski import NumHAcceptors, NumHDonors
from rdkit.Chem import AllChem, Draw, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Crippen import MolLogP
from rdkit.Geometry import Point2D

def _append_featfile(text):
    with open('fm-aa-features.csv', 'a') as f:
        f.write(text)

def print_feats(init=False):

    if not init:
        if not os.path.isdir("AA-Feats"):
            os.mkdir("AA-Feats")
        else:
            subprocess.run(
                f"rm AA-Feats/*",
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )
        
        if not os.path.isdir("FM-Feats"):
            os.mkdir("FM-Feats")
        else:
            subprocess.run(
                f"rm FM-Feats/*",
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )

    constantsdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../constants/"))
    fmnamefile = f'{constantsdir}/fm2names.yaml'
    if os.path.isfile(fmnamefile):
        with open(fmnamefile) as stream:
            fmnamedic = yaml.safe_load(stream)
            

    # PHARMACOPHORE colors
    # aromatic : orange
    # hydrogen donor : blue
    # hydrogen acceptor : red
    # pos ion : yellow
    # neg ion : grey
    # hydrophobic : green

    if not init:
        with open('fm-aa-features.csv', 'w+'):
            pass

        colours = { 'Donor':[62/255,109/255,201/255],
                    'Acceptor':[201/255,62/255,62/255],
                    'NegIonizable':[119/255,116/255,120/255],
                    'PosIonizable':[153/255,62/255,201/255],
                    'Aromatic':[201/255,113/255,62/255],
                    'Hydrophobe':[65/255,201/255,79/255]}

    fdef = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef'))

    fm2sm = fm2smiles()

    aa_smiles = {'ALA': 'C[C@H](N)C=O', 
                'CYS': 'N[C@H](C=O)CS', 
                'ASP': 'N[C@H](C=O)CC(=O)O', 
                'GLU': 'N[C@H](C=O)CCC(=O)O', 
                'PHE': 'N[C@H](C=O)Cc1ccccc1', 
                'GLY': 'NCC=O', 
                'HIS': 'N[C@H](C=O)Cc1c[nH]cn1', 
                'ILE': 'CC[C@H](C)[C@H](N)C=O', 
                'LYS': 'NCCCC[C@H](N)C=O', 
                'LEU': 'CC(C)C[C@H](N)C=O', 
                'MET': 'CSCC[C@H](N)C=O', 
                'ASN': 'NC(=O)C[C@H](N)C=O', 
                'PRO': 'O=C[C@@H]1CCCN1', 
                'GLN': 'NC(=O)CC[C@H](N)C=O', 
                'ARG': 'N=C(N)NCCC[C@H](N)C=O', 
                'SER': 'N[C@H](C=O)CO', 
                'THR': 'C[C@@H](O)[C@H](N)C=O', 
                'VAL': 'CC(C)[C@H](N)C=O', 
                'TRP': 'N[C@H](C=O)Cc1c[nH]c2ccccc12',
                'TYR': 'N[C@H](C=O)Cc1ccc(O)cc1'}
                
    aminonames = {'ALA': 'Alanine', 
                'CYS': 'Cysteine', 
                'ASP': 'Aspartic Acid', 
                'GLU': 'Glutamic Acid', 
                'PHE': 'Phenylalanine', 
                'GLY': 'Glycine', 
                'HIS': 'Histidine', 
                'ILE': 'Isoleucine', 
                'LYS': 'Lysine', 
                'LEU': 'Leucine', 
                'MET': 'Methionine', 
                'ASN': 'Asparagine', 
                'PRO': 'Proline', 
                'GLN': 'Glutamine', 
                'ARG': 'Arginine', 
                'SER': 'Serine', 
                'THR': 'Threonine', 
                'VAL': 'Valine', 
                'TRP': 'Tryptophan',
                'TYR': 'Tyrosine'}

    params = AllChem.ETKDGv3()
    params.useRandomCoords = True
    # params.maxAttempts = 10000
    params.enforceChirality = False

    ms = {}

    if not init:
        print('FM       D   A   NI  PI  AR  HP  LHP  LogP')
        _append_featfile('Functional Monomer, Acronym, D, A, NI, PI, AR, HP, LHP, LogP, SMILES\n')

    for fm in fm2sm:
        smile = fm2sm[fm]

        mol = Chem.MolFromSmiles(smile)
        mol = Chem.rdmolops.AddHs(mol, addCoords=True)
        AllChem.EmbedMolecule(mol, params)
        conf = mol.GetConformer()
        Chem.SanitizeMol(mol)
        # ms.append(mol)

        fmParams = {}
        for k in fdef.GetFeatureFamilies():
            fparams = FeatMaps.FeatMapParams()
            fmParams[k] = fparams

        keep = ('Donor','Acceptor','NegIonizable','PosIonizable','Aromatic', 'Hydrophobe', 'LumpedHydrophobe')
        featLists = []
        rawFeats = fdef.GetFeaturesForMol(mol)
        # filter that list down to only include the ones we're intereted in 
        featLists.append([f for f in rawFeats if f.GetFamily() in keep])

        nums = [0,0,0,0,0,0,0]
        for feats in featLists:
            for feat in feats:
                # print(feat.GetFamily())
                if feat.GetFamily() == 'Donor':
                    nums[0]+=1
                elif feat.GetFamily() == 'Acceptor':
                    nums[1]+=1
                elif feat.GetFamily() == 'NegIonizable':
                    nums[2]+=1
                elif feat.GetFamily() == 'PosIonizable':
                    nums[3]+=1
                elif feat.GetFamily() == 'Aromatic':
                    nums[4]+=1
                elif feat.GetFamily() == 'Hydrophobe':
                    nums[5]+=1
                elif feat.GetFamily() == 'LumpedHydrophobe':
                    nums[6]+=1

        # NumDonors = NumHDonors(mol)
        # NumAcceptors = NumHAcceptors(mol)
        # NumNegIon = NumNegIonizable(mol)

        logp = round(MolLogP(mol),3)

        ms[fm] = []
        for i in nums:
            ms[fm].append(i)
        ms[fm].append(logp)
        
        if not init:
            print(f'{fm}'.ljust(9)+f'{nums[0]}   {nums[1]}   {nums[2]}   {nums[3]}   {nums[4]}   {nums[5]}   {nums[6]}   {logp}')
            _append_featfile(f'{fmnamedic[fm]}, {fm}, {nums[0]}, {nums[1]}, {nums[2]}, {nums[3]}, {nums[4]}, {nums[5]}, {nums[6]}, {logp}, {fm2sm[fm]}\n')


    if init:
        constantsdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../constants/"))
        with open(f'{constantsdir}/fm2pharmfeats.yaml', 'w') as outfile:
            yaml.dump(ms, outfile)

    if not init:
        mols = []
        high = []
        for fm in fm2sm:
            smile = fm2sm[fm]

            mol = Chem.MolFromSmiles(smile)
            mol = Chem.rdmolops.AddHs(mol, addCoords=True)
            AllChem.EmbedMolecule(mol, params)
            conf = mol.GetConformer()
            Chem.SanitizeMol(mol)
            mol = Chem.RemoveAllHs(mol)
            mol3d = deepcopy(mol)

            rdDepictor.Compute2DCoords(mol)
            
            fmParams = {}
            for k in fdef.GetFeatureFamilies():
                fparams = FeatMaps.FeatMapParams()
                fmParams[k] = fparams


            # keep = ('Donor','Acceptor','NegIonizable','PosIonizable','Aromatic', 'Hydrophobe', 'LumpedHydrophobe')
            keep = ('Donor','Acceptor','NegIonizable','PosIonizable', 'Hydrophobe','Aromatic')
            featLists = []
            rawFeats = fdef.GetFeaturesForMol(mol)
            # filter that list down to only include the ones we're intereted in 
            featLists.append([f for f in rawFeats if f.GetFamily() in keep])
            mols.append(mol)

            #
            # https://stackoverflow.com/questions/75735344/highlighting-smiles-features-with-rdkit
            #

            drawer = Draw.MolDraw2DCairo(400,400)
            drawer.DrawMolecule(mol,
                highlightAtoms=[],
                highlightAtomColors=[],
                highlightBonds=[],
                legend=fm)

            for fff in keep:
                for feats in featLists:
                    hit = []
                    for feat in feats:
                        if feat.GetFamily() == fff:
                            atom_ids = feat.GetAtomIds()
                            for atom in atom_ids:
                                atom_point = mol.GetConformer().GetAtomPosition(atom)
                                pt = Point2D(atom_point.x, atom_point.y)
                                # print(atom_point.x, atom_point.y, atom_point.z)
                                if fff == 'Donor':
                                    r,g,b = colours[fff]
                                    drawer.SetColour((r, g, b, 0.5))
                                    drawer.DrawArc( pt, 0.5, 0, 60)
                                elif fff == 'Acceptor':
                                    r,g,b = colours[fff]
                                    drawer.SetColour((r, g, b, 0.5))
                                    drawer.DrawArc( pt, 0.5,60,120)
                                elif fff == 'Hydrophobe':
                                    r,g,b = colours[fff]
                                    drawer.SetColour((r, g, b, 0.5))
                                    drawer.DrawArc( pt, 0.5,120,180)
                                elif fff == 'Aromatic':
                                    r,g,b = colours[fff]
                                    drawer.SetColour((r, g, b, 0.5))
                                    drawer.DrawArc( pt, 0.5,180,240)
                                elif fff == 'NegIonizable':
                                    r,g,b = colours[fff]
                                    drawer.SetColour((r, g, b, 0.5))
                                    drawer.DrawArc( pt, 0.5,240,300)
                                elif fff == 'PosIonizable':
                                    r,g,b = colours[fff]
                                    drawer.SetColour((r, g, b, 0.5))
                                    drawer.DrawArc( pt, 0.5,300,360)
            
            drawer.DrawMolecule(mol,
                highlightAtoms=[],
                highlightAtomColors=[],
                highlightBonds=[],
                legend=fm)
            # drawer.DrawMolecule(mol,
            # 	highlightAtoms=highlight_atoms,
            # 	highlightAtomColors=highlight_colors_atom,
            # 	highlightBonds=highlight_bonds,
            # 	legend=fm)
            drawer.FinishDrawing()
            drawer.WriteDrawingText(f'FM-Feats/{fm}.png')

            # high.append(highlight_colors_atom)

        # print(high)

        img = Draw.MolsToGridImage(mols,molsPerRow=10,legends=fm2sm)#, highlightAtomListsMatrix=high )#
        img.save(f'FMS.png')

        print('-------------------------------------------')
        _append_featfile('Amino Acid, Acronym, D, A, NI, PI, AR, HP, LHP, LogP, SMILES\n')
        for fm in aa_smiles:
            smile = aa_smiles[fm]

            mol = Chem.MolFromSmiles(smile)
            mol = Chem.rdmolops.AddHs(mol, addCoords=True)
            AllChem.EmbedMolecule(mol, params)
            conf = mol.GetConformer()
            Chem.SanitizeMol(mol)
            # ms.append(mol)

            fmParams = {}
            for k in fdef.GetFeatureFamilies():
                fparams = FeatMaps.FeatMapParams()
                fmParams[k] = fparams

            keep = ('Donor','Acceptor','NegIonizable','PosIonizable','Aromatic', 'Hydrophobe', 'LumpedHydrophobe')
            featLists = []
            rawFeats = fdef.GetFeaturesForMol(mol)
            # filter that list down to only include the ones we're intereted in 
            featLists.append([f for f in rawFeats if f.GetFamily() in keep])

            nums = [0,0,0,0,0,0,0]
            for feats in featLists:
                for feat in feats:
                    # print(feat.GetFamily())
                    if feat.GetFamily() == 'Donor':
                        nums[0]+=1
                    elif feat.GetFamily() == 'Acceptor':
                        nums[1]+=1
                    elif feat.GetFamily() == 'NegIonizable':
                        nums[2]+=1
                    elif feat.GetFamily() == 'PosIonizable':
                        nums[3]+=1
                    elif feat.GetFamily() == 'Aromatic':
                        nums[4]+=1
                    elif feat.GetFamily() == 'Hydrophobe':
                        nums[5]+=1
                    elif feat.GetFamily() == 'LumpedHydrophobe':
                        nums[6]+=1

            NumDonors = NumHDonors(mol)
            NumAcceptors = NumHAcceptors(mol)
            # NumNegIon = NumNegIonizable(mol)

            logp = round(MolLogP(mol),3)

            print(f'{fm}'.ljust(9)+f'{nums[0]}   {nums[1]}   {nums[2]}   {nums[3]}   {nums[4]}   {nums[5]}   {nums[6]}   {logp}')
            _append_featfile(f'{aminonames[fm]}, {fm}, {nums[0]}, {nums[1]}, {nums[2]}, {nums[3]}, {nums[4]}, {nums[5]}, {nums[6]}, {logp}, {aa_smiles[fm]}\n')


        mols = []
        high = []
        for fm in aa_smiles:
            smile = aa_smiles[fm]

            mol = Chem.MolFromSmiles(aa_smiles[fm])
            mol = Chem.rdmolops.AddHs(mol, addCoords=True)
            AllChem.EmbedMolecule(mol, params)
            conf = mol.GetConformer()
            Chem.SanitizeMol(mol)
            mol = Chem.RemoveAllHs(mol)
            # mol2d = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
            mol3d = deepcopy(mol)

            rdDepictor.Compute2DCoords(mol)


            fmParams = {}
            for k in fdef.GetFeatureFamilies():
                fparams = FeatMaps.FeatMapParams()
                fmParams[k] = fparams


            # keep = ('Donor','Acceptor','NegIonizable','PosIonizable','Aromatic', 'Hydrophobe', 'LumpedHydrophobe')
            keep = ('Donor','Acceptor','NegIonizable','PosIonizable', 'Hydrophobe','Aromatic')
            featLists = []
            rawFeats = fdef.GetFeaturesForMol(mol)
            # filter that list down to only include the ones we're intereted in 
            featLists.append([f for f in rawFeats if f.GetFamily() in keep])
            mols.append(mol)
            
            drawer = Draw.MolDraw2DCairo(400,400)
            drawer.DrawMolecule(mol,
                highlightAtoms=[],
                highlightAtomColors=[],
                highlightBonds=[],
                legend=fm)

            for fff in keep:
                for feats in featLists:
                    hit = []
                    for feat in feats:
                        if feat.GetFamily() == fff:
                            atom_ids = feat.GetAtomIds()
                            for atom in atom_ids:
                                atom_point = mol.GetConformer().GetAtomPosition(atom)
                                pt = Point2D(atom_point.x, atom_point.y)
                                # print(atom_point.x, atom_point.y, atom_point.z)
                                if fff == 'Donor':
                                    r,g,b = colours[fff]
                                    drawer.SetColour((r, g, b, 0.5))
                                    drawer.DrawArc( pt, 0.5, 0, 60)
                                elif fff == 'Acceptor':
                                    r,g,b = colours[fff]
                                    drawer.SetColour((r, g, b, 0.5))
                                    drawer.DrawArc( pt, 0.5,60,120)
                                elif fff == 'Hydrophobe':
                                    r,g,b = colours[fff]
                                    drawer.SetColour((r, g, b, 0.5))
                                    drawer.DrawArc( pt, 0.5,120,180)
                                elif fff == 'Aromatic':
                                    r,g,b = colours[fff]
                                    drawer.SetColour((r, g, b, 0.5))
                                    drawer.DrawArc( pt, 0.5,180,240)
                                elif fff == 'NegIonizable':
                                    r,g,b = colours[fff]
                                    drawer.SetColour((r, g, b, 0.5))
                                    drawer.DrawArc( pt, 0.5,240,300)
                                elif fff == 'PosIonizable':
                                    r,g,b = colours[fff]
                                    drawer.SetColour((r, g, b, 0.5))
                                    drawer.DrawArc( pt, 0.5,300,360)
            
            drawer.DrawMolecule(mol,
                highlightAtoms=[],
                highlightAtomColors=[],
                highlightBonds=[],
                legend=fm)
            drawer.FinishDrawing()
            drawer.WriteDrawingText(f'AA-Feats/{fm}.png')