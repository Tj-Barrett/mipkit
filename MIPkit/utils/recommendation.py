import os, subprocess
from MIPkit.constants.constants import (
	fm2smiles,)

import yaml

import numpy as np

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
from rdkit import DataStructs

def initialized_feature_matrix():

    constantsdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../constants/"))
    fmnamefile = f'{constantsdir}/fm2names.yaml'
    if os.path.isfile(fmnamefile):
        with open(fmnamefile) as stream:
            fmnamedic = yaml.safe_load(stream)

    fdef = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef'))

    fm2sm = fm2smiles()

    params = AllChem.ETKDGv3()
    params.useRandomCoords = True
    # params.maxAttempts = 10000
    params.enforceChirality = False

    ms = []

    feat_vectors = {}
    ms = []
    for fm in fm2sm:
        smile = fm2sm[fm]

        mol = Chem.MolFromSmiles(smile)
        mol = Chem.rdmolops.AddHs(mol, addCoords=True)
        AllChem.EmbedMolecule(mol, params)
        conf = mol.GetConformer()
        Chem.SanitizeMol(mol)
        ms.append(mol)

    lenms = len(ms)

    fpgen = AllChem.GetMorganGenerator(radius=2)
    # fpgen = AllChem.GetRDKitFPGenerator()
    fps = [fpgen.GetFingerprint(x) for x in ms]

    feat_matrix = np.zeros((len(fm2sm),len(fm2sm)),)

    for i in range(lenms):
        for j in range(lenms):
            feat_matrix[i,j] = DataStructs.TanimotoSimilarity(fps[i], fps[j])
    
    similarity = {}
    for i, afm in enumerate(fm2sm):
        similarity[afm] = feat_matrix[i,:].tolist()

    with open(f'{constantsdir}/fm2similarity.yaml', 'w') as outfile:
        yaml.dump(similarity, outfile)

def find_closest(fm, similarity):
    
    fmvec = similarity[fm]
    # print(fmvec)
    fmvev = list(fmvec)
    keys=list(similarity.keys() )

    fmvec, keys = zip(*sorted(zip(fmvec, keys), reverse=True))

    # print(f"Similarity to : {fm}")
    # print( "--------------------")

    recs = {}
    for i in range(1, min(6, len(keys))):
        recs[keys[i]] = round(fmvec[i],3)
        # print(f"{keys[i]}  :  {round(fmvec[i],3)}")
    return recs

def suggest_new_fms(fms, costs):
    constantsdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../constants/"))
    fmnamefile = f'{constantsdir}/fm2similarity.yaml'
    if os.path.isfile(fmnamefile):
        with open(fmnamefile) as stream:
            similarity = yaml.safe_load(stream)
    else:
        print(f'No similarity file {fmnamefile} available. Exiting...')
        exit()

    for i, f in enumerate(fms):
        print(f"{f}".ljust(10)+"  Similarity  "+"Cost (eur/Mol)".rjust(18) ) #+"Adjusted Cost".rjust(18)
        if costs[f]:
            print(f"                        "+f"{round(costs[f],2)}".rjust(18)) #+f"{round(costs[f]*fms[f],2)}".rjust(18)
        print( "------------------------------------------") #--------------------
        recs = find_closest(f, similarity)
        for rec in recs:
            if recs[rec] > 0.4:
                val = str(recs[rec]).rjust(12)
                if costs[rec] is not None:
                    cost = str(round(costs[rec],2))
                    adjcost = str(round(costs[rec]*fms[f],2))
                else:
                    cost = str(costs[rec])
                    adjcost = str(costs[rec])
                print(rec.ljust(10)+val+cost.rjust(20)) #+adjcost.rjust(18)
        print()
        