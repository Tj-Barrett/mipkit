import yaml
import os
from MIPkit.utils.utils import isBool

def exceptions():
    exs = [ "ALA",
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
            # Template
            "UNK",
            "UNL",
            # Non standard amino acids
            "AIB",
            # Solvents
            "DMF",
            "DMS",
            "CHL",
            "SPC",
            ]
    return exs

def alphas():
    a = '0123456789ABCDEFGHIJKLMNOPQRSTUVWYZ'
    exs = exceptions()
    al = []
    for i in a:
        for j in a:
            if 'F'+i+j not in al and 'A'+i+j not in exs:
                al.append('F'+i+j)
    return al

def import_fm_lib(pargs=None, fmlib=None):
    if fmlib:
        fmlibfile = fmlib
        constantsdir = ""
    else:
        constantsdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../constants/"))
        fmlibfile = f'{constantsdir}/fm-list.yaml'

    if os.path.isfile(fmlibfile):
        with open(fmlibfile) as stream:
            yamlin = yaml.safe_load(stream)

            if yamlin["fm-list"]:
                fm_smiles = {}
                fm_fullname = {}
                fm_mw = {}
                fm_cost = {}
                fm_available = {}
                fm_type = {}
                fm_reaction = {}
                fmlist = yamlin["fm-list"]
                for fm in fmlist:
                    fm_smiles[fm] = fmlist[fm]['smiles']
                    fm_fullname[fm] = fmlist[fm]['fullname']

                    if fmlist[fm]['mw'] and fmlist[fm]['cost'] :
                        fm_mw[fm]   = float(fmlist[fm]['mw'])
                        fm_eur_g = float(fmlist[fm]['cost'])
                        fm_cost[fm] = fm_mw[fm] * fm_eur_g
                    else:
                        fm_cost[fm] = None

                    if fmlist[fm]['type']:
                        fm_type[fm] = fmlist[fm]['type']
                    else:
                        fm_type[fm] = None
                    
                    if isBool(fmlist[fm]['available']) :
                        fm_available[fm] = fmlist[fm]['available']
                    else:
                        fm_available[fm] = None

                    if fmlist[fm]['reaction']:
                        if fmlist[fm]['reaction'] == 'vinyl':
                            fm_reaction[fm] = fmlist[fm]['reaction']
                        else:
                            fm_reaction[fm] = None
                    else:
                        fm_reaction[fm] = None

            if yamlin["fm-equivalency"]:
                fm_equal = yamlin["fm-equivalency"]


        with open(f'{constantsdir}/fm2smiles.yaml', 'w') as outfile:
            yaml.dump(fm_smiles, outfile)
        
        with open(f'{constantsdir}/fm2names.yaml', 'w') as outfile:
            yaml.dump(fm_fullname, outfile)

        with open(f'{constantsdir}/fm2equal.yaml', 'w') as outfile:
            yaml.dump(fm_equal, outfile)

        with open(f'{constantsdir}/fm2costs.yaml', 'w') as outfile:
            yaml.dump(fm_cost, outfile)
        
        with open(f'{constantsdir}/fm2mws.yaml', 'w') as outfile:
            yaml.dump(fm_mw, outfile)

        with open(f'{constantsdir}/fm2available.yaml', 'w') as outfile:
            yaml.dump(fm_available, outfile)

        with open(f'{constantsdir}/fm2type.yaml', 'w') as outfile:
            yaml.dump(fm_type, outfile)
        
        with open(f'{constantsdir}/fm2reaction.yaml', 'w') as outfile:
            yaml.dump(fm_reaction, outfile)

        encode = {}
        decode = {}
        a = alphas()
        for i, fm in enumerate(fmlist):
            encode[fm] = a[i]
            decode[a[i]] = fm

        with open(f'{constantsdir}/encodefms.yaml', 'w') as outfile:
            yaml.dump(encode, outfile)
        
        with open(f'{constantsdir}/decodefms.yaml', 'w') as outfile:
            yaml.dump(decode, outfile)

    else:
        print(f"No FM lib file {fmlibfile} available. Exiting... ")
        quit()