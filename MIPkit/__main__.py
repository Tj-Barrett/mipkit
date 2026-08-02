#!/usr/bin/env python3

from MIPkit import __version__
from MIPkit.constants.handle_fms import import_fm_lib
from MIPkit.constants.print_fms import print_fms, print_smiles
from MIPkit.run.complex import complex
from MIPkit.run.precomplex import precomplex
from MIPkit.run.interaction import interact, regen
from MIPkit.run.screen_template import screen_template
from MIPkit.utils.abn_epitopes import determine_abn, print_abn
from MIPkit.utils.determine_cost import recipe_cost
from MIPkit.utils.print_costs import print_costs
from MIPkit.utils.print_feats import print_feats
from MIPkit.utils.parser import parser
from MIPkit.utils.parser import read_config
from MIPkit.utils.rdkit_utils import generate_sdfs, generate_amino_sdfs, write_fm_pdb
from MIPkit.utils.recommendation import initialized_feature_matrix
from MIPkit.visualize.polyanimation import polyanimation


import time

def splash():
    return "\n=====================================================================\n\tBMBT MIPkit v2026.7\n\tDOI:\n\tBarrett et al.\n\n\tIssues? : tbar@tf.uni-kiel.de\n=====================================================================\n"


def main():
    args = parser()

    if args.timer:
        start = time.time()

    ################################################################
    #
    # init
    #
    ################################################################
    print(splash())

    if args.init:
        print("Initializing MIPkit and generating necessary files.\n")

        print("Reading global config file...")
        read_config(args)

        print("Generating FM Library...")
        import_fm_lib()

        print("Generating Functional Monomers...")
        generate_sdfs()

        print("Generating Amino Acids...")
        generate_amino_sdfs()
        
        print("Generating Functional Monomer Features...")
        print_feats(init=True)

        print("Generating Functional Monomer Similarity Matrix...")
        initialized_feature_matrix()
    
    if args.print_costs:
        print_costs()
    
    if args.print_fms:
        print_fms()

    if args.print_feats:
        print_feats()

    if args.print_smiles:
        print_smiles()

    if args.ABN:
        determine_abn(args)
    
    if args.print_abn:
        print_abn()

    if args.cost:
        recipe_cost(args, suggest=True)


    ################################################################
    #
    # gen_precomplex
    #
    ################################################################

    if args.dock or args.screen:
        precomplex(args)

    if args.screen_fm:
        screen_template(args)

    ################################################################
    #
    # gen_complex
    #
    ################################################################

    if args.react:
        complex(args)

    ################################################################
    #
    # visualize
    #
    ################################################################

    if args.interact:
        interact(args)

    ################################################################
    #
    # visualize
    #
    ################################################################

    if args.polygif:
        polyanimation(args)

    ################################################################
    #
    # Other
    #
    ################################################################

    if args.write_pdb:
        write_fm_pdb(args)

    if args.regen:
        regen(args)

    if args.timer:
        end = time.time()
        total = end-start
        print(
            "%s %s %s"
            % ( "Total time in MIPkit : ",
                str("%3.1f" % total).ljust(12),
                " seconds"
            )
        )

if __name__ == "__main__":
    main()
