import yaml
import os

def import_amino_lib(pargs=None, aminolib=None):
    if aminolib:
        aminolibfile = aminolib
        constantsdir = ""
    else:
        constantsdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../aminos/"))
        aminolibfile = f'{constantsdir}/aminos-list.yaml'

    if os.path.isfile(aminolibfile):
        with open(aminolibfile) as stream:
            yamlin = yaml.safe_load(stream)

            if yamlin["amino-list"]:
                aminolist = yamlin["amino-list"]
        
        return aminolist

    else:
        print(f"No Amino lib file {aminolibfile} available. Exiting... ")
        quit()