import os, sys
import random
from datetime import datetime
import yaml

import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondType

from MIPkit.constants.constants import fm_equivalents, fm2smiles
from MIPkit.utils.read_pdb import read_pdb

################################################################
#
# Is Functions
#
################################################################

def isFloat(n):
    try:
        float(n)
        return True
    except (ValueError, TypeError):
        return False
    # return isinstance(n, float)

def isInt(n):
    try:
        int(n)
        return True
    except (ValueError, TypeError):
        return False
    # return isinstance(n, int)

def isString(s):
    try:
        return isinstance(s[0], str)
    except (TypeError, IndexError):
        return isinstance(s, str)

def isBool(s):
    return isinstance(s, bool)

def isnumber(n):
    try:
        float(n)
        return True
    except ValueError:
        return False

def isRXN(args, n):
    # if default seed
    seed = datetime.now().timestamp()
    implicit_seed = args.implicit_seed[0] if isinstance(args.implicit_seed, list) else args.implicit_seed
    if isFloat(n) and implicit_seed == "@TIME":
        random.seed(seed)
        return n >= random.random()
    elif n == 1:
        return True
    elif n is None:
        return True
    else:
        random.seed(implicit_seed)
        return n >= random.random()


################################################################
#
# Helper Functions
#
################################################################

def generate_config(templatepdb, pargs):

    if templatepdb:
        mols = read_pdb(templatepdb)

        coords = []
        for mol in mols:
            coords.extend(mol["coords"])

        xs = []
        ys = []
        zs = []
        for c in coords:
            xs.append(c[0])
            ys.append(c[1])
            zs.append(c[2])

        avgx = float(round(np.mean(xs), 3))
        avgy = float(round(np.mean(ys), 3))
        avgz = float(round(np.mean(zs), 3))

        maxx = np.max(xs)
        maxy = np.max(ys)
        maxz = np.max(zs)

        minx = np.min(xs)
        miny = np.min(ys)
        minz = np.min(zs)

        # Make docking dimensions a minimum of 10. Particularly important for small molecules that might be planar
        rangex = int(np.ceil( np.max([int(np.ceil(maxx - minx)), 100])))
        rangey = int(np.ceil( np.max([int(np.ceil(maxy - miny)), 100])))
        rangez = int(np.ceil( np.max([int(np.ceil(maxz - minz)), 100])))

        out = {
            "protein": templatepdb,
            "fms": pargs.fms,
            "scale": 1,
            "center": [avgx, avgy, avgz],
            "size": [rangex, rangey, rangez],
            "energy_range": pargs.energy_range,
            "exhaustiveness": pargs.exhaustiveness,
        }
        with open("generated_config.yaml", "w+") as outfile:
            yaml.safe_dump(out, outfile, default_flow_style=False)

        return out
    else:
        print("No template molecule provided. Skipping config generation.")
        sys.exit()

def generate_recipe(afms):
    fms = {}
    recipe = []
    fmsmiles = fm2smiles()
    for i, f in enumerate(afms):
        # if the number is even
        if i % 2 == 0 and f in fmsmiles:
            try:
                fms[f] = int(afms[i + 1])
                recipe.append(f)
            except (IndexError, ValueError, TypeError) as e:
                print(f"Error in FMs. Are you trying to dock without a specific quantity?")
                sys.exit()
        elif not isInt(f) and f not in fmsmiles:
            print(f"Functional Monomer : {f} is not in the library. Check for misspellings, or add it to MIPkit/constants/fm-list.yaml.")
            sys.exit()
        elif isInt(f):
            continue
            

    sanfm = sanitize_fms(recipe)
    tmp = fms
    fms = {}

    for i, sfm in enumerate(sanfm):
        fms[sfm] = tmp[recipe[i]]

    return fms

def generic_top(pargs, pro_itp, protein, itps, mols, posres=False, initiator=False, solvents=None):
    amberdir = pargs.amberdir
    initiatordir = pargs.initiatordir
    solventdir = pargs.solventdir

    initiator_str = ""
    if initiator:
        inits = ['aps','nh4','tmd','nnn','abn','azo']
        for i in inits:
            initiator_str += f'\n#include "{initiatordir}/{i}.itp"'

    solvents_str = ""
    solvents_cnt = ""
    if solvents:
        # for sol in solvents:
        solvents_str += f'#include "{solventdir}/{solvents}.itp"'

    if posres:
        gtop = f"""
; [ defaults ] specified by acpype already included via forcefield.itp
;
# include "{amberdir}/forcefield.itp"

; #include "amber03.ff/forcefield.itp"

{itps}
{pro_itp}
#ifdef POSRES
    #include "{posres}"
#endif

{initiator_str}
{solvents_str}

#include "{amberdir}/spc.itp"
#include "{amberdir}/ions.itp"

[ system ]
RXN in water

[ molecules ]
; Compound        nmols
{protein}
{mols}
"""
    else:
        gtop = f"""
; [ defaults ] specified by acpype already included via forcefield.itp
;
# include "{amberdir}/forcefield.itp"

; #include "amber03.ff/forcefield.itp"

{itps}
{pro_itp}

{initiator_str}
{solvents_str}

#include "{amberdir}/spc.itp"
#include "{amberdir}/ions.itp"

[ system ]
RXN in water

[ molecules ]
; Compound        nmols
{protein}
{mols}
"""
    return gtop

def moving_average(x, N):
    return np.convolve(x, np.ones(N)/N, mode='valid')

def mol_to_dictionary(mol, fm, PT):
    new_atoms = []
    new_coords = []
    new_bonds = {}
    for _m, atom in enumerate(mol.GetAtoms()):
        ipos = mol.GetConformer().GetAtomPosition(_m)
        numa = atom.GetAtomicNum()
        x = round(ipos.x, 3)
        y = round(ipos.y, 3)
        z = round(ipos.z, 3)
        new_coords.append((x, y, z))
        new_atoms.append(PT.GetElementSymbol(numa))
        new_bonds[_m] = []

    for _b, bond in enumerate(mol.GetBonds()):
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        new_bonds[a].append(b)
        if bond.GetBondType() == BondType.DOUBLE:
            new_bonds[a].append(b)
        elif bond.GetBondType() == BondType.TRIPLE:
            new_bonds[a].append(b)
            new_bonds[a].append(b)

    new_fm = {"name": "", "number": "", "atoms": [], "coords": [], "bonds": []}

    new_fm["name"] = fm  # fm['name']
    new_fm["number"] = 1
    new_fm["atoms"] = new_atoms
    new_fm["coords"] = new_coords
    new_fm["bonds"] = new_bonds
    return new_fm


def sanitize_fms(recipe):
    fm_equals = fm_equivalents()
    fme = []
    for fm in fm_equals:
        fme.append(fm)
    fme = set(fme)

    new_recipe = []
    for fm in recipe:
        if fm in fme:
            new_recipe.append(fm_equals[fm])
        else:
            new_recipe.append(fm)

    return new_recipe

def sanitize_protein(protein):
    # print(protein)
    mol = Chem.MolFromPDBFile(
        protein, sanitize=False, removeHs=False, proximityBonding=True
    )
    # mol = Chem.AddHs(mol)

    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol, 1000)

    w = Chem.SDWriter(protein + ".sdf")
    try:
        w.write(mol)
    finally:
        w.close()