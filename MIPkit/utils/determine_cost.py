import yaml, os, sys
import numpy as np
from MIPkit.constants.constants import decode_fms, encode_fms
from MIPkit.utils.parser import sanitize_inputs
from MIPkit.utils.utils import  isInt
from MIPkit.utils.recommendation import suggest_new_fms


def determine_recipe_cost(pargs, suggest=False):
    constantsdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../constants/"))
    fmcostfile = f'{constantsdir}/fm2costs.yaml'
    if os.path.isfile(fmcostfile):
        with open(fmcostfile) as stream:
            fmcostdic = yaml.safe_load(stream)
    else:
        print(f'No cost file {fmcostfile} available. Exiting...')
        exit()

    fmmwfile = f'{constantsdir}/fm2mws.yaml'
    if os.path.isfile(fmmwfile):
        with open(fmmwfile) as stream:
            fmmwdic = yaml.safe_load(stream)
    else:
        print(f'No molecular weight file {fmmwfile} available. Exiting...')
        exit()

    defm = decode_fms()
    enfm = encode_fms()

    recipe = []
    fms = {}
    for i, f in enumerate(pargs.fms):
        if f in defm:
            recipe.append(defm[f])
            if defm[f] in fms:
                fms[defm[f]]+=1
            else:
                fms[defm[f]]=1
        else:
            # if the number is even
            if i % 2 == 0 and f in enfm:
                try:
                    fms[f] = int(pargs.fms[i + 1])
                    recipe.append(f)
                except (IndexError, ValueError, TypeError) as e:
                    print(f"Error in FMs. Are you trying to determine cost without a specific quantity?")
                    sys.exit()
            elif isInt(f):
                continue
    
    totalcost = 0.0
    totalmw = 0.0
    totalcostmw = 0.0
    omitted = []
    n = 0
    for fm in fms:
        cost = fmcostdic.get(fm)
        mw = fmmwdic.get(fm)
        if cost is not None and mw is not None:
            cost = round(fms[fm]*cost/1000, 2)
            totalmw += mw
            totalcost+=cost
            n+=1
        else:
            omitted.append(fm)

    if n == 0:
        print('No functional monomers detected. Check complex name or -fms entries.')
    else:
        totalcost = totalcost/n
        totalcostmw = totalcost/totalmw*1000

    if suggest:
        suggest_new_fms(fms, fmcostdic)
    
    print(f'€ Estimated Recipe Cost (eur/mMol) : {round(totalcost,2)}')
    print(f'€                          (eur/g) : {round(totalcostmw,2)}')
    if len(omitted) > 0 :
        print(f'  FMs {omitted} omitted from cost analysis. These can be updated by modifying "MIPkit/constants/fm-list.yaml".')
    print(f'\n')


def recipe_cost(args, suggest=False):

    pargs = sanitize_inputs(args)

    determine_recipe_cost(pargs, suggest)