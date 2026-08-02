import yaml, os
from MIPkit.constants.constants import decode_fms, encode_fms
from MIPkit.utils.utils import  isInt

def check_availability(pargs):
    constantsdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../constants/"))
    fmavailablefile = f'{constantsdir}/fm2available.yaml'
    if os.path.isfile(fmavailablefile):
        with open(fmavailablefile) as stream:
            fmavailabledic = yaml.safe_load(stream)
    else:
        print(f'No availability file {fmavailablefile} available. Exiting...')
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
                    print(f'- Warning : Could not parse a quantity for FM entry "{f}" in -fms. Skipping it.\n')
                    continue
            elif isInt(f):
                continue
    
    for fm in fms:
        if fm not in fmavailabledic or fmavailabledic[fm] is None:
            print(f'- Warning : According to the fm list, {fm} is not available commercially. If you can synthesis it in-house, you can ignore this message.\n')