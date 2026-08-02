import argparse, os, glob, multiprocessing, sys, shlex, subprocess
import numpy as np
from MIPkit.utils.utils import (
    isBool,
    isFloat,
    isInt,
    isString,
)
from MIPkit.utils.read_config import read_config
from MIPkit.utils.read_yaml import read_yaml

def parser():
    parser = argparse.ArgumentParser(
        prog='MIPkit', usage='%(prog)s [-screen/-dock/-react/-interact/...] [options]',
        description="MIPkit : Automated screening, docking, precomplexation, and complexation."
    )

    #===========================================================
    #
    # Main Args
    #
    #===========================================================
    group_required = parser.add_argument_group('Required Commands', 'These are the basic MIPkit functions.')
    required = group_required.add_mutually_exclusive_group()

    required.add_argument(
        "-init", 
        help="Initialize and generate pdbs.", 
        action="store_true"
    )
    required.add_argument(
        "-dock", 
        help="Dock recipe using vina or gnina.", 
        action="store_true"
    )
    required.add_argument(
        "-react", 
        help="Complexation reaction.", 
        action="store_true"
    )
    required.add_argument(
        "-screen", 
        help="Screen for functional monomers via docking.", 
        action="store_true"
    )
    required.add_argument(
        "-screen_fm", 
        "-gmxscreen",
        help="Screen for functional monomers via docking and molecular dynamics.", 
        action="store_true"
    )
    required.add_argument(
        "-interact",
        help="Run MIP/NIP interaction.",
        action="store_true"
    )


    #===========================================================
    #
    # Templates and Complex
    #
    #===========================================================
    group_features = parser.add_argument_group('Template and Complex', 'Provide a template molecule for docking, screening, or complexation. Use -protein for proteins in general PDB format (uses pdb2gmx), -gmxprotein for proteins in GROMACS PDB format, and -template for non-proteins. Use -template/-protein and -complex for MIP production, and just -complex for NIP production.')
    template = group_features.add_mutually_exclusive_group()

    template.add_argument(
        "-protein",
        help="Template protein in PDB format.",
        nargs=1,
        default=False,
        type=str
    )
    template.add_argument(
        "-gmxprotein",
        help="Template protein in PDB format that has already been processed by GROMACS.",
        nargs=1,
        default=False,
        type=str
    )
    template.add_argument(
        "-template",
        help="Template molecule.",
        nargs=1, 
        default=False,
        type=str
    )

    group_features.add_argument(
        "-complex",
        "-cplx",
        help="Precomplexation pdb for complexation.",
        nargs=1,
        default=False,
        type=str
    )
    #===========================================================
    #
    # Docking and Screening
    #
    #===========================================================
    group_dock_screen = parser.add_argument_group('Docking and Screening', 'These options inform GNINA/Vina on how to operate.')

    group_dock_screen.add_argument(
        "-dockmethod", 
        help="Dock method using vina or gnina.", 
        default='gnina',
        choices = ['gnina', 'vina'],
        type=str,
    )

    group_dock_screen.add_argument(
        "-energy_range",
        help="Energy range for Vina.",
        nargs=1,
        default=[4],
        type=int
    )
    group_dock_screen.add_argument(
        "-exhaustiveness",
        help="Exhaustiveness for Vina and GNINA.",
        nargs=1,
        default=[10],
        type=int
    )
    group_dock_screen.add_argument(
        "-fms", 
        help="FM and the number of them.", 
        nargs="+",
        type=str,
        default = None
    )
    group_dock_screen.add_argument(
        "-scale", 
        help="Scale the recipe.", 
        nargs=1, 
        default=[1],
        type=int
    )
    group_dock_screen.add_argument(
        "-config", 
        help="Use the config file to read.", 
        nargs=1,
        default=None,
        type=str
    )
    group_dock_screen.add_argument("-shuffle_seed",
        help="Shuffle FM seed number.",
        nargs=1,
        default="xF00D",
        type=str
        )
    group_dock_screen.add_argument(
        "-shuffle",
        help="Shuffle FM order.",
        action="store_true"
        )

    #===========================================================
    #
    # Reaction
    #
    #===========================================================
    group_react = parser.add_argument_group('Reaction Options', 'These options are determine -react behaviors.')

    group_react.add_argument(
        "-box",
        help="Box dimensions for GROMACS. Default is 10.0 10.0 10.0",
        nargs=3,
        default=[10.0, 10.0, 10.0],
        type=float
    )
    group_react.add_argument(
        "-cutoff", 
        help="Cutoff distance for bonding. First is used for C-C distances, second is used for C-APS distances.",
        nargs='+',
        default=[3.3, 3.3],
        type=float
    )
    group_react.add_argument(
        "-cycles", 
        help="Number of compute cycles to run.", 
        nargs=1, 
        default=[10],
        type=int
    )
    group_react.add_argument(
        "-min", 
        help="Minimize input pdb before complexation.", 
        action="store_true"
    )
    group_react.add_argument(
        "-reactiontype",
        help="Polymerization type. Currently just vinylization is supported.",
        nargs=1,
        default="vinyl",
        type=str
    )
    group_react.add_argument(
        "-shell",
        help="Use a shell script (.sh) instead of default GROMACS steps.",
        nargs=1,
        default=False,
        type=str
    )
    group_react.add_argument(
        "-solvent",
        help="Solvent type. Defaults to spc216",
        nargs=1,
        default=None,
        choices = ['dmf', 'dms', 'chl', 'mcn', 'spc'],
        type=str
    )
    group_react.add_argument(
        "-posre",
        help="GMX position restraint in [kcal/mol].",
        nargs=1,
        default=[1000],
        type=int
    )
    group_react.add_argument("-explicit",
        help="Explicit polymerization flag. Takes solvent acronym and mol ratio. Options: APS, TEMED, ACPA, AIBN ",
        nargs="+",
        default=False,
        type=str
    )
    group_react.add_argument("-implicit",
        help="Implicit polymerization flag. Takes a reaction probability number between 0 and 1.",
        nargs=1,
        default=False,
        type=float
    )
    group_react.add_argument("-implicit_seed",
        help="Implicit polymerization seed number.",
        nargs=1,
        default="@TIME",
        type=str
    )
    group_react.add_argument(
        "-restart",
        help="Restart from a previous step. This must be a valid step.",
        action="store_true"
        )
    group_react.add_argument(
        "-accharge",
        help="Charge derivation method. 'bcc' is Bond Charge Correction and 'gas' is Gasteiger",
        nargs=1,
        default='gas',
        choices = ['gas', 'bcc'],
        type=str
    )
    group_react.add_argument(
        "-acff",
        help="Forcefield to use for with ACPYPE.",
        nargs=1,
        default='gaff2',
        choices = ['gaff', 'gaff2'], #amber and amber2 are also options but they dont work
        type=str
    )
    group_react.add_argument(
        "-dt",
        help="Timestep in femtoseconds for minimization, reactions, and interactions. Takes up to two arguments, either (rxn & min dt), or (min dt, rxn dt)", 
        nargs='+',
        default=['2'],
        type=str
    )
    group_react.add_argument(
        "-relax",
        help="Relax the protein structure if running into ACPYPE issues.",
        action="store_true"
    )
    group_react.add_argument(
        "-ncpu", 
        help="Number of cpu cores to use.", 
        nargs=1, 
        default=False,
        type=int
    )
    group_react.add_argument(
        "-gpu_id",
        help="Specify GPU id number",
        nargs=1,
        default=[0],
        type=int
    )
    group_react.add_argument(
        "-parallel",
        help="Modify the temp directories for parallel runs. Use A, B, C, etc. with -gpu_id 1, 2, 3, etc.",
        nargs=1,
        default="",
        type=str
    )
    group_react.add_argument(
        "-force",
        help="Force ACPYPE to progress if atoms are too close.",
        action="store_true"
    )
    group_react.add_argument(
        "-nocap",
        help="Don't cap with hydrogens during implicit reactions, allowing for faster polymerization.",
        action="store_true"
    )
    group_react.add_argument(
        "-temp",
        help="Temperature in K.",
        nargs=1,
        default=[300.0],
        type=float
    )

    #===========================================================
    #
    # Interaction Args
    #
    #===========================================================
    group_interact = parser.add_argument_group('Interact Options', 'These options are determine -interact behaviors.')

    group_interact.add_argument(
        "-offset",
        help="Epitope offset for interaction with the NIP.",
        nargs = 3,
        default = False,
        type=str
    )
    group_interact.add_argument(
        "-wash", 
        help="Wash complexes to get MIP and NIPs.", 
        action="store_true"
    )
    group_interact.add_argument(
        "-regen",
        help="Regenerate forcefields.", 
        action="store_true"
    )
    group_interact.add_argument(
        "-id",
        "-identify",
        help="Decompose Polymers into FMs",
        action="store_true"
    )
    group_interact.add_argument(
        "-time", 
        help="Interaction time length in ns. Takes up to two arguments for minimization and rxn or interaction", 
        nargs='+',
        default=['100'],
        type=str
    )


    #===========================================================
    #
    # General Utils
    #
    #===========================================================
    group_utils = parser.add_argument_group('Utils', 'These options are general utilities and dev tools.')

    group_utils.add_argument(
        "-basedir", 
        help="Base working directory.", 
        nargs=1, 
        default=os.getcwd(),
        type=str
    )

    group_utils.add_argument(
        "-ABN",
        help="Determine acidic/basic/neutral breakdown of protein molecule.",
        action="store_true"
        )
    
    group_utils.add_argument(
        "-cost",
        help="Determine cost of a recipe in eur/mMol using either a config or -fms.",
        action="store_true"
        )
    
    group_utils.add_argument(
        "-print_abn",
        help="Print Acidic/Basic/Neutral breakdown of Amino Acids.",
        action="store_true"
        )
    
    group_utils.add_argument(
        "-print_costs",
        help="Print estimated costs of all functional monomers in eur/mol.",
        action="store_true"
        )

    group_utils.add_argument(
        "-print_feats",
        help="Print all FM features",
        action="store_true"
        )

    group_utils.add_argument(
        "-print_fms",
        help="Print all fms",
        action="store_true"
        )

    group_utils.add_argument(
        "-print_smiles",
        help="Print all SMILES",
        action="store_true"
        )

    group_utils.add_argument("-noclean",
        help="Don't clean up. This is only suitable for finding errors in the first few cycles if you are modifying the package.",
        action="store_true"
    )
    
    group_utils.add_argument(
        "-timer",
        help="Time event.",
        action="store_true"
    )

    group_utils.add_argument(
        "-verbose", 
        help="Turn on verbose mode.", 
        action="store_true"
    )
    
    group_utils.add_argument(
        "-write_pdb",
        help="Write a pdb of the functional monomer.",
        nargs=1,
        default=False,
        type=str
    )
    
    #===========================================================
    #
    # Visualization Args
    #
    #===========================================================
    group_visualize = parser.add_argument_group('Visualization Options', 'These options are for the basic visualization tools included.')
    group_visualize.add_argument(
        "-polygif", 
        help="Polymerization gif.", 
        action="store_true"
    )
    group_visualize.add_argument(
        "-xstats",
        help="Polymerization statistics file. Defaults to work directory. ",
        nargs=1,
        default=False,
        type=str
    )
    group_visualize.add_argument(
        "-mwstats",
        help="Molecular Weight file. Defaults to work directory. ",
        nargs=1,
        default=False,
        type=str
    )
    group_visualize.add_argument(
        "-gifname", 
        help="Gif file to print.", 
        nargs=1, 
        default=False,
        type=str
    )
    group_visualize.add_argument(
        "-giftime", 
        help="Gif time in seconds.", 
        nargs=1, 
        default=[5],
        type=int
    )
    group_visualize.add_argument(
        "-plotstep",
        help="Plot step time (x axis step) in nanoseconds.",
        nargs=1,
        default=[1],
        type=int
    )

    args = parser.parse_args()

    # Add obabel and vina locations
    args = read_config(args)

    return args

def sanitize_inputs(args):
    # ===========================================================
    #
    # We need to do os.path.exists(...)
    #
    # ===========================================================

    pargs = parsed_args()

    #===========================================================
    #
    # Applications
    #
    #===========================================================
    if isString(args.gmx):
        pargs.gmx=args.gmx
    else:
        print("- GMX alias is not specified in the config file, defaulting to gmx.")
        pargs.gmx="gmx"

    if isString(args.acpype):
        pargs.acpype=args.acpype
    else:
        print("- Acpype alias is not specified in the config file, defaulting to acpype.")
        pargs.acpype="acpype"

    if isString(args.obabel):
        pargs.obabel=args.obabel
    else:
        print("- Obabel alias is not specified in the config file, defaulting to obabel.")
        pargs.obabel="obabel"

    if isString(args.vina):
        pargs.vina=args.vina
    else:
        print("- Vina alias is not specified in the config file, defaulting to vina.")
        pargs.vina="vina"

    if isString(args.gnina):
        pargs.gnina=args.gnina
    else:
        print("- gnina alias is not specified in the config file, defaulting to gnina.")
        pargs.gnina="gnina"

    # ===========================================================
    #
    # Main Run Options
    #
    # ===========================================================

    if args.dockmethod:
        pargs.dockmethod = args.dockmethod
    else:
        print('No docking method specified. This should not happen.')
        exit()

    if args.dock:
        pargs.dock = True
        print(f"+ Docking with {pargs.dockmethod}")
    else:
        pargs.dock = False

    if args.min:
        print(f"+ Precomplex minimization is selected.")
        pargs.min=True
    else:
        pargs.min=False

    if args.interact:
        pargs.interact = True
        print("+ Interaction run. If you want experimental mip/nip structures, use -wash.")
    else:
        pargs.interact = False

    if args.react:
        pargs.react = True
        print("+ Reaction run.")
    else:
        pargs.react = False

    if args.regen:
        pargs.regen=True
        print("+ Regenerate forcefields.")
    else:
        pargs.regen = False

    if args.screen:
        pargs.screen = True
        print(f"+ Screening with {pargs.dockmethod}")
    else:
        pargs.screen = False

    if args.wash and args.interact:
        pargs.wash=True
        print("+ System will be rinsed of any unreacted functional monomers.")
    else:
        pargs.wash = False

    if args.restart:
        pargs.restart = True
    else:
        pargs.restart = False

    print('\n\nSelected Options:\n')
    # ===========================================================
    #
    # Directories
    #
    # ===========================================================

    if args.ABN:
        print("+ Acidic/Basic/Neutral breakdown of template protein.")

    if args.basedir and isinstance(args.basedir, list) and isString(args.basedir[0]):
        pargs.basedir = args.basedir[0]
        basedir = args.basedir[0]
    elif args.basedir and isString(args.basedir):
        pargs.basedir = args.basedir
        basedir = args.basedir
    else:
        print('We need a base directory. This should not have failed. Maybe a string issue?')
        exit()

    if not basedir.endswith(os.sep):
        basedir += os.sep
    pargs.basedir = basedir
    basedir_q = shlex.quote(basedir)

    if isBool(args.noclean):
        if args.noclean:
            pargs.clean = False
        else:
            pargs.clean = True
    else:
        print("\nError : Something is wrong; -noclean is a store true condition and did not return a bool")
        exit()

    if pargs.clean and not pargs.interact and not pargs.restart:
        if glob.glob(f"{basedir}\#*"):
            subprocess.run(f"rm {basedir_q}\\#*", shell=True)
        if glob.glob(f"{basedir}step**"):
            subprocess.run(f"rm {basedir_q}step*", shell=True)
        if glob.glob(f"{basedir}posre*"):
            subprocess.run(f"rm {basedir_q}posre*", shell=True)

    if not os.path.isdir(f"{basedir}work"):
        os.mkdir(f"{basedir}work")
    else:
        if glob.glob(f"{basedir}work/*") and args.react and not pargs.restart:
            subprocess.run(f"rm {basedir_q}work/*", shell=True)
    pargs.workdir = f"{basedir}work/"

    if not os.path.isdir(f"{basedir}comp"):
        os.mkdir(f"{basedir}comp")
    else:
        if glob.glob(f"{basedir}comp/*") and args.react and not pargs.restart:
            subprocess.run(f"rm {basedir_q}comp/*", shell=True)
    pargs.compdir = f"{basedir}comp/"

    if not os.path.isdir(f"{basedir}logs"):
        os.mkdir(f"{basedir}logs")
    else:
        if glob.glob(f"{basedir}logs/*") and args.react and not pargs.restart:
            subprocess.run(f"rm {basedir_q}logs/*", shell=True)
    pargs.logdir = f"{basedir}logs/"

    if not os.path.isdir(f"{basedir}RXN.acpype"):
        os.mkdir(f"{basedir}RXN.acpype")
    else:
        if glob.glob(f"{basedir}RXN.acpype/*") and args.react and not pargs.restart:
            subprocess.run(f"rm -r {basedir_q}RXN.acpype", shell=True)

    if not os.path.isdir(f"{basedir}RES.acpype"):
        os.mkdir(f"{basedir}RES.acpype")
    else:
        if glob.glob(f"{basedir}RES.acpype/*") and args.react and not pargs.restart:
            subprocess.run(f"rm -r {basedir_q}RES.acpype", shell=True)

    # gromacs scripts
    pargs.gmxdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../gmx/"))

    # amber inputs
    pargs.amberdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../amber/"))

    # sulfur inputs
    pargs.sulfurdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../sulfurs/"))

    # initiator inputs
    pargs.initiatordir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../initiators/"))

    # solvent inputs
    pargs.solventdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../solvents/"))

    #===========================================================
    #
    # Targets
    #
    #===========================================================

    # Precomplex to be polymerized
    if isString(args.complex):
        if os.path.isfile(args.complex[0]):
            pargs.cplx = args.complex[0]
            pargs.incomplex = args.complex[0]
        else:
            print(
                f"No complex file available. Exiting... "
            )
            quit()
    elif pargs.react or pargs.interact:
        print(
            "\nError : This function requires a precomplex to polymerize.\nYou might have a typo, or use -dock flag to create one.\n"
        )
        quit()

    # Protein template, uses pdb2gmx
    if isString(args.protein):
        if os.path.isfile(args.protein[0]):
            pargs.protein = args.protein[0]
            pargs.gmxprotein = False
            pargs.template_name = args.protein[0]
            pargs.intemplate = args.protein[0]
        else:
            print(
                f"No template file available. Exiting... "
            )
            quit()
    elif isString(args.gmxprotein):
        if os.path.isfile(args.gmxprotein[0]):
            pargs.protein = args.gmxprotein[0]
            pargs.gmxprotein = True
            pargs.template_name = args.gmxprotein[0]
            pargs.intemplate = args.gmxprotein[0]
        else:
            print(
                f"No template file available. Exiting... "
            )
            quit()
    else:
        pargs.protein = False
        pargs.gmxprotein = False

    # -template
    # Non-protein template
    if isString(args.template):
        if os.path.isfile(args.template[0]):
            pargs.template = args.template[0]
            pargs.template_name = args.template[0]
            pargs.intemplate = args.template[0]
        else:
            print(
                f"No template file available. Exiting... "
            )
            quit()
    else:
        pargs.template = False

    if args.id:
        pargs.id = True
        print("+ Identify Monomers in Polymers")
    else:
        pargs.id = False
    #===========================================================
    #
    # Reaction Arguments
    #
    #===========================================================

    timesteps = {'interact':0.001, # 1 fs
                 'min':None,
                 'min_run':0.001,
                 'npt':0.001,
                 'nvt':0.0005,
                 'nvt_vacuum':0.0001,
                 'nvt_stage0':0.0001,
                 'nvt_stage1':0.0005,
                 'nvt_stage2':0.001,
                 'nvt_stage3':0.0015,
                 'rxn':0.001}
    
    nsteps =    {'interact': 50000000, # 50 ns
                 'min': 200000,
                 'min_run': 10000000,
                 'npt': 200000, # 0.2 ns
                 'nvt': 500000, # 0.2 ns
                 'nvt_vacuum':200000,
                 'nvt_stage0':500000,
                 'nvt_stage1':200000,
                 'nvt_stage2':200000,
                 'nvt_stage3':200000,
                 'rxn': 500000}

    # rxn : [rxn]
    # min, rxn : [min, rxn]
    # interact : [interact]
    # min, interact : [min, interact]

    if isinstance(args.dt, list):
        if args.react and not args.min:
            if len(args.dt) == 1:
                timesteps['rxn'] = float(args.dt[0])/1000
            else:
                print(f'Too many timesteps {args.dt}')
                timesteps['rxn'] = float(args.dt[0])/1000

        elif args.react and args.min:
            if len(args.dt) == 1:
                timesteps['rxn'] = float(args.dt[0])/1000
                timesteps['min_run'] = float(args.dt[0])/1000
            elif len(args.dt) == 2:
                timesteps['rxn'] = float(args.dt[1])/1000
                timesteps['min_run'] = float(args.dt[0])/1000
            else:
                print(f'Too many timesteps {args.dt}')
                timesteps['rxn'] = float(args.dt[1])/1000
                timesteps['min_run'] = float(args.dt[0])/1000

        elif (args.interact and not args.min) or (args.screen and not args.min) or (args.screen_fm and not args.min):
            if len(args.dt) == 1:
                timesteps['interact'] = float(args.dt[0])/1000
            else:
                print(f'Too many timesteps {args.dt}')
                timesteps['interact'] = float(args.dt[0])/1000
            
        elif (args.interact and args.min) or (args.screen and args.min) or (args.screen_fm and args.min):
            if len(args.dt) == 1:
                timesteps['interact'] = float(args.dt[0])/1000
                timesteps['min_run'] = float(args.dt[0])/1000
            elif len(args.dt) == 2:
                timesteps['interact'] = float(args.dt[1])/1000
                timesteps['min_run'] = float(args.dt[0])/1000
            else:
                print(f'Too many timesteps {args.dt}')
                timesteps['interact'] = float(args.dt[1])/1000
                timesteps['min_run'] = float(args.dt[0])/1000

        for _dtkey, _dtval in timesteps.items():
            if _dtval is not None and _dtval <= 0:
                print(f'\nError with -dt : {args.dt}. Timestep must be a positive number. Exiting ...')
                exit()

        pargs.dt = timesteps

    else:
        print(f'Timesteps loaded as a float {args.dt}. Using defaults.')
        pargs.dt = timesteps

    # *1000 for conversion to go from fs of dt to ns of time
    if isinstance(args.time, list):
        if args.react and not args.min:
            if len(args.time) == 1:
                nsteps['rxn'] = int(float(args.time[0])/timesteps['rxn']*1000)
            else:
                nsteps['rxn'] = int(float(args.time[0])/timesteps['rxn']*1000)
                print(f'Too many times {args.time[0]}')

        elif args.react and args.min:
            if len(args.time) == 1:
                nsteps['rxn'] = int(float(args.time[0])/timesteps['rxn']*1000)
                nsteps['min_run'] = int(float(args.time[0])/timesteps['min_run']*1000)
            elif len(args.time) == 2:
                nsteps['rxn'] = int(float(args.time[1])/timesteps['rxn']*1000)
                nsteps['min_run'] = int(float(args.time[0])/timesteps['min_run']*1000)
            else:
                nsteps['rxn'] = int(float(args.time[1])/timesteps['rxn']*1000)
                nsteps['min_run'] = int(float(args.time[0])/timesteps['min_run']*1000)
                print(f'Too many times {args.time}')

        elif (args.interact and not args.min) or (args.screen and not args.min) or (args.screen_fm and not args.min):
            if len(args.time) == 1:
                nsteps['interact'] = int(float(args.time[0])/timesteps['interact']*1000)
            else:
                nsteps['interact'] = int(float(args.time[0])/timesteps['interact']*1000)
                print(f'Too many times {args.time}')
            
        elif (args.interact and args.min)  or (args.screen and args.min) or (args.screen_fm and args.min):
            if len(args.time) == 1:
                nsteps['interact'] = int(float(args.time[0])/timesteps['interact']*1000)
                nsteps['min_run'] = int(float(args.time[0])/timesteps['min_run']*1000)
            elif len(args.time) == 2:
                nsteps['interact'] = int(float(args.time[1])/timesteps['interact']*1000)
                nsteps['min_run'] = int(float(args.time[0])/timesteps['min_run']*1000)
            else:
                nsteps['interact'] = int(float(args.time[1])/timesteps['interact']*1000)
                nsteps['min_run'] = int(float(args.time[0])/timesteps['min_run']*1000)
                print(f'Too many times {args.time}')

        for _tskey, _tsval in nsteps.items():
            if _tsval is not None and _tsval <= 0:
                print(f'\nError with -time : {args.time}. Time must be a positive number. Exiting ...')
                exit()

        pargs.nsteps = nsteps

    else:
        pargs.nsteps = {'interact': 50000000, # 50 ns
                        'min': 100000,
                        'min_run': 10000000,
                        'npt': 500000, # 0.5 ns
                        'npt_vacuum':100000,
                        'nvt': 500000, # 0.5 ns
                        'nvt_vacuum':100000,
                        'rxn': 500000}

    # number of cycles to run
    if isInt(args.cycles[0]) and pargs.react:
        if not args.restart:
            try:
                # in nanoseconds?
                if int(args.cycles[0]) <= 0:
                    print("- Cycles must be a positive integer. Defaulting to 10 cycles.")
                    pargs.ncycles = 10
                else:
                    pargs.ncycles = int(args.cycles[0])
            except (ValueError, IndexError, TypeError) as e:
                print("- Could not convert cycles to int. Defaulting to 10 cycles.")
                pargs.ncycles = 10
        else:
            # get number from the complex, e.g. "step23-fullcomplex-prod.pdb" -> 23
            basefilename = os.path.basename(pargs.cplx)
            renum = ""
            endrenum = False
            for i in basefilename:
                if i.isnumeric() and not endrenum:
                    renum += i
                elif i=='-'and not endrenum:
                    endrenum = True
                elif endrenum:
                    break

            try:
                renum = int(renum)
            except ValueError:
                print(
                    f"\nError : Could not determine the step number to restart from in "
                    f"'{basefilename}'.\n-restart expects -cplx to be a previous step's "
                    f"output, e.g. 'stepN-fullcomplex-prod.pdb'. Exiting ..."
                )
                exit()

            if int(args.cycles[0]) <= 0:
                print("- Cycles must be a positive integer. Defaulting to 10 cycles.")
                pargs.ncycles = 10
            else:
                pargs.ncycles = int(args.cycles[0])

            if pargs.ncycles <= renum + 1:
                print(
                    f"\nError : -restart resumes after step {renum}, but -cycles "
                    f"{pargs.ncycles} does not leave any further steps to run. "
                    f"Pass a larger -cycles value. Exiting ..."
                )
                exit()

            pargs.restart_number = renum
    
    #  polymerization types
    pargs.initiator = []

    if not args.implicit and not args.explicit and pargs.react:
        # Default, not flagging either
        print("+ Both Explicit and Implicit flags are false. Assuming you want implicit.")
        pargs.explicit = False
        pargs.implicit = True
        pargs.rxn_prob = 1
    elif args.implicit and not args.explicit and pargs.react:
        # Implicit, tell it, it's not explicit
        pargs.explicit = False
        pargs.implicit = True
    elif args.explicit and not args.implicit and pargs.react:
        # Explicit, reiterate its not implicit
        pargs.explicit = True
        pargs.implicit = False
        ###
        if len(args.explicit) % 2 == 1:
            print('- Initiator count is not even. Both an initiator and the count are needed.')
            exit()

        inits = {}
        aae = args.explicit
        for i, f in enumerate(aae):
            # if the number is even
            if i % 2 == 0:
                try:
                    count = int(aae[i + 1])
                except (ValueError, IndexError):
                    print(f'- Initiator count for {f} is not a valid integer. Exiting...')
                    exit()
                if count <= 0:
                    print(f'- Initiator count for {f} must be a positive integer. Exiting...')
                    exit()

                if f == "APS" or f == "aps":
                    inits["aps"] = 2*count
                    inits["nh4"] = 2*count
                elif f == "TMD" or f == "tmd" or f == "TEMED" or f == "temed" or f == "TMEDA" or f == "tmeda":
                    inits["tmd"] = count
                elif f == "AZO" or f == "azo" or f == "ACPA" or f == "acpa":
                    inits["azo"] = 2 * count
                    inits["nnn"] = count
                elif f == "AIBN" or f == "aibn" or f == 'ABN' or f == "abn":
                    inits["abn"] = 2 * count # AIB is a non-proteinogenic amino acid
                    inits["nnn"] = count
                elif f == "NN" or f == "nn":
                    print('- NN is computed automatically.')
                else:
                    print(f'- Warning: unrecognized initiator "{f}" in -explicit. It will be ignored.')

        pargs.initiator = inits
        pargs.ininitiator = []

    elif args.explicit and args.implicit and args.react:
        print("- Both Explicit and Implicit flags are found. Please select one.")
        exit()
    else:
        pargs.explicit = False
        pargs.implicit = False


    if args.cutoff and pargs.implicit and pargs.react:
        if isFloat(args.cutoff[0]) and float(args.cutoff[0]) > 0:
            pargs.cutoff_CC = float(args.cutoff[0])
            pargs.cutoff_CO = 3.3
        else:
            print("- Could not convert cutoff to a positive float. Defaulting to 3.3 A.")
            pargs.cutoff_CC  = 3.3
            pargs.cutoff_CO = 3.3

        if args.implicit and isFloat(args.implicit[0]):
            if 0 <= float(args.implicit[0]) <= 1:
                pargs.rxn_prob  = float(args.implicit[0])
                print(f"+ Implicit reaction with probability {pargs.rxn_prob} and cutoff {pargs.cutoff_CC}.")
            else:
                print("- Implicit reaction probability is not between 0 and 1. Defaulting to 1.")
                pargs.rxn_prob = 1
        else:
            print("- Implicit reaction will be set to default of reaction probability 1.")
            pargs.rxn_prob = 1
    elif args.cutoff and pargs.explicit:
        if isinstance(args.cutoff, list):
            if len(args.cutoff) == 2 and isFloat(args.cutoff[0]) and isFloat(args.cutoff[1]) and float(args.cutoff[0]) > 0 and float(args.cutoff[1]) > 0 and pargs.react:
                pargs.cutoff_CC = float(args.cutoff[0])
                pargs.cutoff_CO = float(args.cutoff[1])
                print(f"+ Explicit reaction with cutoffs {pargs.cutoff_CC} A and {pargs.cutoff_CO} A.")
            elif len(args.cutoff) == 1 and isFloat(args.cutoff[0]) and float(args.cutoff[0]) > 0 and pargs.react:
                pargs.cutoff_CC = float(args.cutoff[0])
                pargs.cutoff_CO = float(args.cutoff[0])
                print(f"+ Explicit reaction using cutoff {pargs.cutoff_CC} A.")
            elif pargs.react:
                print("- Could not convert both cutoffs to a positive float. Defaulting to 3.3 A and 3.3 A.")
                pargs.cutoff_CC = 3.3
                pargs.cutoff_CO = 3.3
        else:
            print("- Could not convert both cutoffs to a positive float. Defaulting to 3.3 A and 3.3 A.")
            pargs.cutoff_CC = 3.3
            pargs.cutoff_CO = 3.3
    else:
        pass

    # Temperature
    if isinstance(args.temp, list):
        pargs.temp = args.temp[0]
    elif isinstance(args.temp, float):
        pargs.temp = args.temp
    else:
        print(f"+ Defaulting to 300 K")
        pargs.temp = 300

    if not isFloat(pargs.temp) or float(pargs.temp) <= 0:
        print(f"- Temperature {pargs.temp} K is not a positive number. Defaulting to 300 K")
        pargs.temp = 300

    # ===========================================================
    #
    # MD Variables - GROMACS
    #
    # ===========================================================

    # acpype charge derivation
    if isinstance(args.accharge, list):
        pargs.charge = args.accharge[0]
    elif isinstance(args.accharge, str):
        pargs.charge = args.accharge
    else:
        print(f"\nNo charge method specified for ACPYPE. This should not have happened. Exiting ...")
        exit()
    
    # acpype forcefield
    if isinstance(args.acff, list):
        pargs.acff = args.acff[0]
    elif isinstance(args.acff, str):
        pargs.acff = args.acff
    else:
        print(f"\nNo forcefield specified for ACPYPE. This should not have happened. Exiting ...")
        exit()

    # -box
    if args.box:
        if (
            isFloat(args.box[0]) and isFloat(args.box[1]) and isFloat(args.box[2])
            and float(args.box[0]) > 0 and float(args.box[1]) > 0 and float(args.box[2]) > 0
        ):
            box_x = float(args.box[0])
            box_y = float(args.box[1])
            box_z = float(args.box[2])
            box_str = f"-box {box_x} {box_y} {box_z}"
            pargs.box = box_str
        else:
            print(f"Box for input {args.box} failed (must be three positive numbers). Falling back to -box 10 10 10.")
            box_str = f"-box 10 10 10"
            pargs.box = box_str

    if args.posre:
        if isInt(args.posre[0]) and int(args.posre[0]) > 0:
            pargs.posre = int(args.posre[0])
        else:
            pargs.posre = 1000
            print(f"Position restraint is not a positive integer. Defaulting to 1000.")

    if isString(args.solvent):
        pargs.solvent = args.solvent[0]
        ensol = {'dmf':"Dimethylformamide (DMF)",
                 'chl':"Chloroform (CHL)",
                 'dms':"Dimethylsulfoxide (DMSO)",
                 'mcn':"Acetonitrile (MeCN)",
                 'spc':"SPC Water"}
        print(f"+ System is using { ensol[pargs.solvent] } as the solvent.")
    else:
        pargs.solvent = False

    if args.offset:
        if isFloat(args.offset[0]) and isFloat(args.offset[1]) and isFloat(args.offset[2]):
            ox = float(args.offset[0])
            oy = float(args.offset[1])
            oz = float(args.offset[2])
            pargs.offset = [ox, oy, oz]
        else:
            print(f"\nError with offset : {args.offset}. Exiting ...")
            exit()
    else:
        pargs.offset = False

    # ===========================================================
    #
    # Runtime Arguments - GROMACS
    #
    # ===========================================================
    if args.gpu_id:
        if isInt(args.gpu_id[0]) and int(args.gpu_id[0]) >= 0:
            pargs.gpu_id = f"-gpu_id {int(args.gpu_id[0])}"
        else:
            print("GPU ID failing for %s. Falling back to -gpu_id 0" % args.gpu_id)
            pargs.gpu_id = f"-gpu_id 0"
    else:
        pargs.gpu_id = ""

    if args.ncpu:
        if isInt(args.ncpu[0]) and int(args.ncpu[0]) > multiprocessing.cpu_count():
            print(f"\nError with ncpu : {args.ncpu}. Too many of CPU cores are requested. CPU only reports {multiprocessing.cpu_count()}. Exiting ...")
            exit()
        elif isInt(args.ncpu[0]) and int(args.ncpu[0]) > 0 and int(args.ncpu[0]) <= multiprocessing.cpu_count():
            pargs.ncpu = int(args.ncpu[0])
        else:
            print(f"\nError with ncpu : {args.ncpu}. Number of CPU cores needs to be a positive integer. Exiting ...")
            exit()
    else:
        # use a quarter of the cpu cores
        count = multiprocessing.cpu_count() / 4
        corelist = np.array([1, 2, 4, 8, 16, 32, 64, 128])
        n = np.abs(corelist - count).argmin()
        pargs.ncpu = corelist.flat[n]

    if len(args.parallel) > 0 and isString(args.parallel[0]):
        pargs.parallel = "-" + args.parallel[0]
    else:
        pargs.parallel = ""

    # ===========================================================
    #
    # Do Later
    #
    # ===========================================================
    if args.config:
        if args.config and os.path.isfile(args.config[0]):
            args = read_yaml(args)
            pargs.config = args.config[0]
        else:
            print(f'Recipe config {args.config[0]} does not exist.')
            sys.exit()
    
    if args.nocap:
        pargs.nocap = True
    else:
        pargs.nocap = False

    if args.force:
        pargs.force = args.force

    if args.fms:
        pargs.fms = args.fms

    if args.screen_fm:
        pargs.screen_fms = True

    if isinstance(args.scale, list):
        pargs.scale = args.scale[0]
    else:
        pargs.scale = args.scale

    if not isInt(pargs.scale) or int(pargs.scale) <= 0:
        print(f"- Scale {pargs.scale} is not a positive integer. Defaulting to 1.")
        pargs.scale = 1

    if args.shuffle:
        pargs.shuffle = True
    else:
        pargs.shuffle = False
        

    if args.shuffle_seed:
        if isinstance(args.shuffle_seed, list):
            pargs.shuffle_seed = args.shuffle_seed[0]
        else:
            pargs.shuffle_seed = args.shuffle_seed

    # -shell
    # For not using built in gmx method
    if args.shell:
        pargs.shell = args.shell[0]
    else:
        pargs.shell = False

    if isinstance(args.energy_range, list):
        pargs.energy_range = args.energy_range[0]
    else:
        pargs.energy_range = args.energy_range
    if not isInt(pargs.energy_range) or int(pargs.energy_range) <= 0:
        print(f"- Energy range {pargs.energy_range} is not a positive integer. Defaulting to 4.")
        pargs.energy_range = 4

    if isinstance(args.exhaustiveness, list):
        pargs.exhaustiveness = args.exhaustiveness[0]
    else:
        pargs.exhaustiveness = args.exhaustiveness
    if not isInt(pargs.exhaustiveness) or int(pargs.exhaustiveness) <= 0:
        print(f"- Exhaustiveness {pargs.exhaustiveness} is not a positive integer. Defaulting to 10.")
        pargs.exhaustiveness = 10

    # Just a spacer
    print("\n")
    return pargs

class parsed_args:
    def __init__(self):
        # Applications
        self.acpype = None
        self.gmx = None
        self.obabel = None
        self.vina = None
        self.gnina = None
        self.dockmethod = None

        # Run Options
        self.dock = None
        self.interact = None
        self.react = None
        self.regen = None
        self.screen = None
        self.restart = None
        self.restart_number = None
        self.force = None

        self.nocap = None

        # Screening FMS
        self.screen_fms = None
        self.ns = None

        # Directories
        self.amberdir = None
        self.compdir = None
        self.gmxdir = None
        self.initiatordir = None
        self.logdir = None
        self.solventdir = None
        self.sulfurdir = None
        self.workdir = None

        # Objects
        self.fms = None
        self.cplx = None
        self.protein = None
        self.gmxprotein = None
        self.template = None
        self.solvent = None
        self.initiator = None

        self.template_name = None

        self.config = None

        self.id = None
        self.recipe = None

        # Gromacs MDP 
        self.dt = None
        self.nsteps = None
        self.temp = None

        # Reaction Elements
        self.explicit = None
        self.implicit = None

        self.cutoff_CC = None
        self.cutoff_CO = None
        self.rxn_prob = None

        self.ncycles = None
        self.gmxtype = None
        self.mintype = None

        self.box = None
        self.offset = None
        self.posre = None

        self.incomplex = None
        self.intemplate = None
        self.ininitiator = None

        self.charge = None
        self.acff = None

        # Interact Elements
        self.wash = None

        # Runtime
        self.ncpu = None
        self.gpu_id = None
        self.parallel = None
        self.retry = True

        self.screen_list = None

def sanitize_visual(args):
    # ===========================================================
    #
    # We need to do os.path.exists(...)
    #
    # ===========================================================

    pargs = parsed_viz()

    if args.basedir and isString(args.basedir[0]):
        pargs.basedir = args.basedir[0]
        basedir = pargs.basedir
    else:
        print('We need a base directory. This should not have failed. Maybe a string issue?')
        exit()

    if not basedir.endswith(os.sep):
        basedir += os.sep
    pargs.basedir = basedir

    if not os.path.isdir(f"{basedir}viz"):
        os.mkdir(f"{basedir}viz")
    pargs.workdir = f"{basedir}viz/"

    return pargs

class parsed_viz:
    def __init__(self):
        self.basedir = None
        self.workdir = None

        self.wildcard = None

        self.filename = None
