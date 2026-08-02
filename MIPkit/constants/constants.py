"""
Functional monomers

fm_list() -- list of all available functional monomers
fm2smiles -- convert fm acronym to smiles
smiles2fm -- convert smiles to fm acronym

fmequivalent -- convert any common acronyms to the ones we use here.
                mostly just housekeeping in case I forget what I was
                calling them / we havent settled yet on a final acr.


"""
import yaml
import os

def fm_list():
    constantsdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../constants/"))
    fmlibfile = f'{constantsdir}/fm2smiles.yaml'

    if os.path.isfile(fmlibfile):
        with open(fmlibfile) as stream:
            fm2smiles = yaml.safe_load(stream)
    else:
        print("No fm2smiles yaml file found. Make sure to -init first.")
        exit()

    return fm2smiles

def fm_equivalents():
    constantsdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../constants/"))
    fmlibfile = f'{constantsdir}/fm2equal.yaml'

    if os.path.isfile(fmlibfile):
        with open(fmlibfile) as stream:
            fm2equal = yaml.safe_load(stream)
    else:
        print("No fm2equal yaml file found. Make sure to -init first.")
        exit()

    return fm2equal

def fm2smiles():
    # smiles from
    # https://drugs.ncats.io/drug/
    # https://pubchem.ncbi.nlm.nih.gov
    # suppliers, sigma, fisher, merck, etc
    constantsdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../constants/"))
    fmlibfile = f'{constantsdir}/fm2smiles.yaml'

    if os.path.isfile(fmlibfile):
        with open(fmlibfile) as stream:
            fm2smiles = yaml.safe_load(stream)
    else:
        print("No fm2smiles yaml file found. Make sure to -init first.")
        exit()

    return fm2smiles

def exclusions():
    standards = {
        #Initiators
        "APS":"APS",
        "TMD":"TMD",
        "AZO":"AZO",
        "ABN":"ABN",
        "NNN":"NNN",
        "NH4":"NH4",
    }
    return standards

def encode_fms():
    # PDB only writes first three letters, so encode to a known arbitrary string
    standards = {
        #Constants
        "RXN": "RXN",
        #Initiators
        "APS":"APS",
        "TMD":"TMD",
        "AZO":"AZO",
        "ABN":"ABN",
        "NNN":"NNN",
        "NH4":"NH4",
        #Solvents
        "DMS":"DMS",
        "CHL":"CHL",
        "MCN":"MCN",
        "DMF":"DMF",
        # Template
        "UNK":"UNK",
    }
    constantsdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../constants/"))
    encode_fms = f'{constantsdir}/encodefms.yaml'

    if os.path.isfile(encode_fms):
        with open(encode_fms) as stream:
            encodefms = yaml.safe_load(stream)
    
        fullencodefms = {**encodefms, **standards}
    else:
        print("No encode yaml file found. Make sure to -init first.")
        exit()

    return fullencodefms

def decode_fms():
    standards = {
        #Constants
        "RXN": "RXN",
        #Initiators
        "APS":"APS",
        "TMD":"TMD",
        "AZO":"AZO",
        "ABN":"ABN",
        "NNN":"NNN",
        "NH4":"NH4",
        #Solvents
        "DMS":"DMS",
        "CHL":"CHL",
        "MCN":"MCN",
        "DMF":"DMF",
        # Template
        "UNK":"UNK",
    }
    constantsdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../constants/"))
    decode_fms = f'{constantsdir}/decodefms.yaml'

    if os.path.isfile(decode_fms):
        with open(decode_fms) as stream:
            decodefms = yaml.safe_load(stream)
    
        fulldecodefms = {**decodefms, **standards}
    else:
        print("No decode yaml file found. Make sure to -init first.")
        exit()

    return fulldecodefms

def sulfur_exception():
    enfm = encode_fms()

    ampsa = enfm["AMPSA"]

    sulfurs = {
        ampsa: "AMPSA",
        "AMPSA": ampsa,
    }
    return sulfurs

def solvent_equivalents():
    solvents = {'dmso':'dms'}
    return solvents

"""
Protein residues

res_list --  a full residue list to sanity check against
res2oneshot -- convert resname to one letter representation
oneshot2res -- convert one letter representation back to resname

"""

def res_list():
    res_list = [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "MAS",
        # Templates
        "UNK",
    ]
    return res_list

def res2oneshot():
    res2oneshot = {
        "ALA": "A",
        "ARG": "R",
        "ASN": "N",
        "ASP": "D",
        "CYS": "C",
        "GLN": "Q",
        "GLU": "E",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LEU": "L",
        "LYS": "K",
        "MET": "M",
        "PHE": "F",
        "PRO": "P",
        "SER": "S",
        "THR": "T",
        "TRP": "W",
        "TYR": "Y",
        "VAL": "V",
    }
    return res2oneshot

def oneshot2res():
    oneshot2res = {
        "A": "ALA",
        "R": "ARG",
        "N": "ASN",
        "D": "ASP",
        "C": "CYS",
        "Q": "GLN",
        "E": "GLU",
        "G": "GLY",
        "H": "HIS",
        "I": "ILE",
        "L": "LEU",
        "K": "LYS",
        "M": "MET",
        "F": "PHE",
        "P": "PRO",
        "S": "SER",
        "T": "THR",
        "W": "TRP",
        "Y": "TYR",
        "V": "VAL",
    }
    return oneshot2res

def res_atoms_pdb():
    cartoon_res = {'ALA':['N', 'CA', 'C', 'O', 'CB'],
                   'ARG':['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
                   'ASN':['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'ND2'],
                   'ASP':['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'OD2'],
                   'CYS':['N', 'CA', 'C', 'O', 'CB', 'SG'],
                   'GLN':['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'NE2'],
                   'GLU':['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2'],
                   'GLY':['N', 'CA', 'C', 'O'],
                   'HIS':['N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'],
                   'ILE':['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1'],
                   'LEU':['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2'],
                   'LYS':['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ'],
                   'MET':['N', 'CA', 'C', 'O', 'CB', 'CG', 'SD', 'CE'],
                   'PHE':['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
                   'PRO':['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD'],
                   'SER':['N', 'CA', 'C', 'O', 'CB', 'OG'],
                   'THR':['N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2'],
                   'TRP':['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
                   'TYR':['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'],
                   'VAL':['N', 'CA', 'C', 'O', 'CB', 'CG1','CG2']}
    return cartoon_res

def res_atoms_gmx():
    # pdb2gmx sorts to a different order!
    cartoon_res = {'ALA':['N', 'CA', 'CB', 'C', 'O'],
                   'ARG':['N', 'CA', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2', 'C', 'O'],
                   'ASN':['N', 'CA', 'CB', 'CG', 'OD1', 'ND2', 'C', 'O'],
                   'ASP':['N', 'CA', 'CB', 'CG', 'OD1', 'OD2', 'C', 'O'],
                   'CYS':['N', 'CA', 'CB', 'SG', 'C', 'O'],
                   'GLN':['N', 'CA', 'CB', 'CG', 'CD', 'OE1', 'NE2', 'C', 'O'],
                   'GLU':['N', 'CA', 'CB', 'CG', 'CD', 'OE1', 'OE2', 'C', 'O'],
                   'GLY':['N', 'CA', 'C', 'O'],
                   'HIS':['N', 'CA', 'CB', 'CG', 'ND1', 'CE1', 'NE2', 'CD2', 'C', 'O'],
                   'ILE':['N', 'CA', 'CB', 'CG2', 'CG1', 'CD', 'C', 'O'],
                   'LEU':['N', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'C', 'O'],
                   'LYS':['N', 'CA', 'CB', 'CG', 'CD', 'CE', 'NZ', 'C', 'O'],
                   'MET':['N', 'CA', 'CB', 'CG', 'SD', 'CE', 'C', 'O'],
                   'PHE':['N', 'CA', 'CB', 'CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2', 'C', 'O'],
                   'PRO':['N', 'CD', 'CG', 'CB', 'CA', 'C', 'O'],
                   'SER':['N', 'CA', 'CB', 'OG', 'C', 'O'],
                   'THR':['N', 'CA', 'CB', 'CG2', 'OG1', 'C', 'O'],
                   'TRP':['N', 'CA', 'CB', 'CG', 'CD1', 'NE1', 'CE2', 'CZ2', 'CH2', 'CZ3', 'CE3', 'CD2', 'C', 'O'],
                   'TYR':['N', 'CA', 'CB', 'CG', 'CD1', 'CE1', 'CZ', 'OH', 'CE2', 'CD2', 'C', 'O'],
                   'VAL':['N', 'CA', 'CB', 'CG1','CG2', 'C', 'O']}
    return cartoon_res

def gmx_backbone():
    backbone = {'ALA':[0, 1, 3, 4],
                'ARG':[0, 1, 9, 10],
                'ASN':[0, 1, 6, 7],
                'ASP':[0, 1, 6, 7],
                'CYS':[0, 1, 4, 5],
                'GLN':[0, 1, 7, 8],
                'GLU':[0, 1, 7, 8],
                'GLY':[0, 1, 2, 3],
                'HIS':[0, 1, 7, 8],
                'ILE':[0, 1, 5, 6],
                'LEU':[0, 1, 5, 6],
                'LYS':[0, 1, 6, 7],
                'MET':[0, 1, 5, 6],
                'PHE':[0, 1, 9, 10],
                'PRO':[0, 4, 5, 6],
                'SER':[0, 1, 4, 5],
                'THR':[0, 1, 5, 6],
                'TRP':[0, 1, 12, 13],
                'TYR':[0, 1, 10, 11],
                'VAL':[0, 1, 5, 6]}
    return backbone

"""
Initiators

"""
def initiator_info():
    initiators = {'APS':{'ATT':['O','O'],
                         'IND':[ 4 , 5 ]},
                  # 'ACPA': {'ATT': ['O', 'O'], #AZO
                  #         'IND': [4, 5]},
                  }
    return initiators


"""
Energy Paramters

atom_2_numbers -- number imbedding for atom types
atom_masses -- mass of each atom
atom_eng    -- nonbonded energies for probing

"""

def atom2PT():
    # Acpype doesn't recognize these, so just change them over to generic
    # This should also probably be done algorithmically, but low priority
    # 2025/08/29 - re.sub used to remove any numbers from HETATMS
    # https://www.random.org/strings/?num=100&len=3&digits=on&upperalpha=on&unique=on&format=html&rnd=new
    atom2PT = {
        "CA": "C",
        "CB": "C",
        "CD": "C",
        "CD1": "C",
        "CD2": "C",
        "CE": "C",
        "CE1": "C",
        "CE2": "C",
        "CE3": "C",
        "CG": "C",
        "CG1": "C",
        "CG2": "C",
        "CH2": "C",
        "CZ": "C",
        "CZ2": "C",
        "CZ3": "C",
        "H1": "H",
        "H2": "H",
        "H3": "H",
        "HA": "H",
        "HA1": "H",
        "HA2": "H",
        "HB": "H",
        "HB1": "H",
        "HB2": "H",
        "HB3": "H",
        "HD": "H",
        "HD1": "H",
        "HD11":"H",
        "HD12":"H",
        "HD13":"H",
        "HD2": "H",
        "HD21":"H",
        "HD22":"H",
        "HD23":"H",
        "HD3": "H",
        "HD31": "H",
        "HD32": "H",
        "HD33":"H",
        "HE": "H",
        "HE1": "H",
        "HE2": "H",
        "HE21": "H",
        "HE22": "H",
        "HG": "H",
        "HG1": "H",
        "HG11": "H",
        "HG12": "H",
        "HG13": "H",
        "HG2": "H",
        "HG21": "H",
        "HG22": "H",
        "HG23": "H",
        "HH": "H",
        "HH11": "H",
        "HH12": "H",
        "HH22": "H",
        "HH21": "H",
        "HN": "H",
        "HZ": "H",
        "1HZ": "H",
        "2HZ": "H",
        "3HZ": "H",
        "HZ1": "H",
        "HZ2": "H",
        "HZ3": "H",
        "H1": "H",
        "H2": "H",
        "H3": "H",
        "H4": "H",
        "ND1": "N",
        "ND2": "N",
        "NE": "N",
        "NE1": "N",
        "NE2": "N",
        "NH1": "N",
        "NH2": "N",
        "NZ": "N",
        "N1": "N",
        "N2": "N",
        "N3": "N",
        "N4": "N",
        "OC1": "O",
        "OC2": "O",
        "OD1": "O",
        "OD2": "O",
        "OE1": "O",
        "OE2": "O",
        "OG": "O",
        "OG1": "O",
        "OH": "O",
        "OT1": "O",
        "OT2": "O",
        "OXT": "O",
        "O1": "O",
        "O2": "O",
        "SD": "S",
        "SG": "S",
        "S1": "S",
        "S2": "S",
    }
    return atom2PT


def compressibility():
    # isothermal compressibility, bar^(-1), 1 atm and 300K
    compress = {'spc':'4.5E-5',
                'dmso':"",
                'dmf':"",
                'mecn':"1.15E-4"} #target density ~.777, https://cdnsciencepub.com/doi/pdf/10.1139/v76-398
    return compress

def initiator_radius():
    # initiator radius to be done by QCM at some point
    initrad = {"APS":5,
               "AZO":5}
    return initrad