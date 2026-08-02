# our work
from MIPkit.utils.parser import sanitize_inputs
from MIPkit.utils.read_pdb import (
    read_pdb,
)
from MIPkit.run.run_acpype import acpype_do
from MIPkit.run.run_gmx import gmx_do
from MIPkit.utils.utils import generic_top, generate_recipe
from MIPkit.utils.read_yaml import read_yaml

from MIPkit.run.run_obabel import obabel_do, obabel_combine
from MIPkit.utils.identify import identify
################################################################
#
# Main
#
################################################################
def interact(args):
    if args.id:
        args = read_yaml(args,recipe_only=True)
    pargs = sanitize_inputs(args)

    if args.id:
        pargs.recipe = generate_recipe(pargs.fms)

    if pargs.template_name and pargs.cplx:
        precomplex = read_pdb(pargs.cplx, hydrogens=True)
    else:
        print("Interaction needs both a template and complex.")
        exit()

    fms = []
    int_complex = []
    if pargs.wash:
        for mol in precomplex:
            if mol['name'][0] == "X":
                int_complex.append(mol)
                fms.append(mol['name'])
    else:
        for mol in precomplex:
            int_complex.append(mol)
            fms.append(mol['name'])

    pargs.fms = fms

    if len(fms) < 1:
        print("There are no functional monomers in the system. There is nothing to do.")
        exit()

    # Obabel
    logname = f"obabel-interaction.log"
    out_complex = obabel_do("interact", -1, int_complex, pargs, logname)

    # Acpype
    logname = f"acpype-interaction.log"
    acpype_do("interact", -1, int_complex, out_complex, pargs, logname)

    # Combined
    logname = f"combination-interaction.log"
    obabel_combine("interact", -1, out_complex, pargs, logname)

    # GMX
    logname = f"gmx-interaction.log"
    gmx_run = gmx_do()
    gmx_run.params("interact", -1, -1, pargs, logname)
    gmx_run.run()

def regen(args):
    pargs = sanitize_inputs(args)

    if pargs.template_name:
        precomplex = read_pdb(pargs.cplx, hydrogens=False)
    else:
        print("Regeneration needs either a template and complex.")
        exit()

    # Obabel
    logname = f"obabel-interaction.log"
    out_complex = obabel_do("interact", -1, precomplex, pargs, logname)

    # Acpype
    logname = f"acpype-interaction.log"
    acpype_do("interact", -1, precomplex, out_complex, pargs, logname)
