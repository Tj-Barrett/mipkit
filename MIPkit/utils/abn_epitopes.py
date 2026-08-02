from MIPkit.constants.constants import res2oneshot
from MIPkit.utils.read_pdb import read_pdb
from MIPkit.utils.parser import sanitize_inputs

def res2ABN():
    res2ABN = {
        "ALA": "ALIPHATIC",
        "ARG": "BASIC",
        "ASN": "POLAR",
        "ASP": "ACIDIC",
        "CYS": "POLAR",
        "GLN": "POLAR",
        "GLU": "ACIDIC",
        "GLY": "NEUTRAL",
        "HIS": "BASIC",
        "ILE": "ALIPHATIC",
        "LEU": "ALIPHATIC",
        "LYS": "BASIC",
        "MET": "ALIPHATIC",
        "PHE": "AROMATIC",
        "PRO": "NEUTRAL",
        "SER": "POLAR",
        "THR": "POLAR",
        "TRP": "AROMATIC",
        "TYR": "AROMATIC",
        "VAL": "ALIPHATIC",
    }
    return res2ABN

def oneshot2ABN():
    oneshot2ABN = {
        "A": "ALIPHATIC",
        "R": "BASIC",
        "N": "POLAR",
        "D": "ACIDIC",
        "C": "POLAR",
        "Q": "POLAR",
        "E": "ACIDIC",
        "G": "NEUTRAL",
        "H": "BASIC",
        "I": "ALIPHATIC",
        "L": "ALIPHATIC",
        "K": "BASIC",
        "M": "ALIPHATIC",
        "F": "AROMATIC",
        "P": "NEUTRAL",
        "S": "POLAR",
        "T": "POLAR",
        "W": "AROMATIC",
        "Y": "AROMATIC",
        "V": "ALIPHATIC",
    }
    return oneshot2ABN


# https://www.sigmaaldrich.com/DE/en/technical-documents/technical-article/protein-biology/protein-structural-analysis/amino-acid-reference-chart
def determine_abn(args):

    pargs = sanitize_inputs(args)

    if not pargs.protein:
        print('Acidic/Basic/Neutral breakdown needs a protein structure.')
        exit()
    else:
        template = read_pdb(pargs.protein)
        
        res2one = res2oneshot()
        o2abn = oneshot2ABN()

        oneshotstring = ""
        for aa in template:
            if aa['name'] not in res2one:
                print(f"- Warning: unrecognized residue \"{aa['name']}\" skipped in ABN breakdown.")
                continue
            oneshotstring+= res2one[aa['name']]

        acidic = 0
        basic  = 0
        neutral= 0
        aliphatic = 0
        aromatic = 0
        polar = 0
        for R in oneshotstring:
            if o2abn[R] == "ACIDIC":
                acidic+= 1
            elif o2abn[R] == "BASIC":
                basic += 1
            elif o2abn[R] == "AROMATIC":
                aromatic += 1
            elif o2abn[R] == "ALIPHATIC":
                aliphatic += 1
            elif o2abn[R] == "POLAR":
                polar += 1
            else:
                neutral += 1

        print(f"Template : {oneshotstring}\n \
        \rACIDIC    : {acidic}\n \
        \rBASIC     : {basic} \n \
        \rALIPHATIC : {aliphatic} \n \
        \rAROMATIC  : {aromatic} \n \
        \rPOLAR     : {polar} \n \
        \rNEUTRAL   : {neutral}\n\n\
        \rTOTAL     : {len(oneshotstring)}")

def print_abn():
    res2one = res2oneshot()
    res2abn = res2ABN()


    _strout = 'Amino Acid'.ljust(7)+'L'.center(7)+'Type'.ljust(13)
    print(_strout)
    print('---------------------------')
    for R in res2one:
        _strout = R.ljust(10)+res2one[R].center(7)+res2abn[R].ljust(13)
        print(_strout)
    print()