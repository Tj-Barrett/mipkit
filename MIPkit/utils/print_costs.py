import yaml, os

def print_costs():
    constantsdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../constants/"))
    fmcostfile = f'{constantsdir}/fm2costs.yaml'
    if os.path.isfile(fmcostfile):
        with open(fmcostfile) as stream:
            fmcostdic = yaml.safe_load(stream)
    else:
        print(f'No cost file {fmcostfile} available. Exiting...')
        exit()

    print('Functional Monomer | Cost ')
    print('                   | (eur/mMol)      (eur/Mol)')
    print('----------------------------------------------')
    for fm in fmcostdic:
        cost = fmcostdic[fm]
        if cost is not None:
            fcost = round(cost, 5)
            mcost = round(cost/1000, 5)
            fcost = "{:.2f}".format(fcost)
            mcost = "{:.2f}".format(mcost)
        else:
            fcost = str(cost)
            mcost = str(cost)
        print(f'{fm}'.ljust(19)+mcost.rjust(11)+fcost.rjust(15))