
import random
import shlex
import subprocess
import numpy as np
import os
from rdkit.Chem import GetPeriodicTable
from rich.progress import Progress, TimeElapsedColumn, BarColumn, TextColumn
# our work
from MIPkit.constants.constants import fm2smiles
from MIPkit.utils.check_availability import check_availability
from MIPkit.utils.check_reaction import check_reaction
from MIPkit.utils.determine_cost import determine_recipe_cost
from MIPkit.utils.utils import generate_recipe, generate_config, isnumber
from MIPkit.utils.read_yaml import read_yaml
from MIPkit.utils.read_config import read_config
from MIPkit.utils.parser import sanitize_inputs
from MIPkit.utils.print_molecule import print_molecule, final_print, final_fm_print
from MIPkit.vina.vina import autodock_vina, preprocessing, update_vina
from MIPkit.vina.gnina import gnina, update_gnina


def precomplex(args):
    columns = [
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TimeElapsedColumn(),
    ]

    sdfdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../sdfs/"))

    with Progress(*columns) as pb:
        ################################################################
        #
        # Docking
        #
        ################################################################
        if args.dock:

            if args.config:
                args = read_yaml(args)
            pargs = sanitize_inputs(args)

            if pargs.fms:

                # read recipe from input
                recipe = generate_recipe(pargs.fms)

                check_availability(pargs)
                check_reaction(pargs)
                determine_recipe_cost(pargs)

                serialrecipe = []
                for fm in recipe:
                    for i in range(recipe[fm] * pargs.scale):
                        serialrecipe.append(fm)

                if pargs.shuffle:
                    random.seed(pargs.shuffle_seed)
                    random.shuffle(serialrecipe)

                PT = GetPeriodicTable()

                success = 0

                # generate the config
                vina_args = generate_config(pargs.template_name, args)
            
                if pargs.dockmethod == 'vina':

                    if not os.path.isdir("tmpvina"):
                        os.mkdir("tmpvina")
                    else:
                        subprocess.run(
                            f"rm tmpvina/*",
                            shell=True,
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.STDOUT,
                        )
                    
                    if pargs.template_name:
                        templatepdb = pargs.template_name
                        templatepdbqt = preprocessing(args, receptor=templatepdb)
                    else:
                        print("No template specified. Must have one for docking.")
                        exit()
                        
                    print(
                    "%s %s %s %s %s"
                        % (
                            "Interval".center(9),
                            "FM".center(7),
                            "Mean RMSD".center(12),
                            "Affinity".center(9),
                            "Score".center(9),
                        )
                    )

                    pb1 = pb.add_task(
                        "[steel_blue1]Precomplexation", total=len(serialrecipe)
                    )
                    
                    for i, fm in enumerate(serialrecipe):
                        logdir = "tmpvina"
                        logname = "vina.log"

                        ligand =f'tmpvina/{fm}.pdbqt'

                        subprocess.run(
                            f"{shlex.quote(pargs.obabel)} -isdf {shlex.quote(sdfdir)}/{shlex.quote(fm)}.sdf -opdbqt -O {shlex.quote(ligand)}",
                            shell=True,
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.STDOUT,
                        )

                        out, stats = autodock_vina(log=logname, 
                                    step = i,
                                    pargs=pargs, 
                                    template=templatepdbqt, 
                                    fm=ligand, 
                                    vina_args=vina_args)

                        # take the lowest score, most negative (highest absolute value) affitinity and lowest mean rmsd
                        affinity = []
                        mean_rmsd = []
                        weighted_rmsd = []
                        sanda_score = []
                        for j, stat in enumerate(stats):
                            if j == 0:
                                # excluded (self-referential rmsd), but keep the
                                # index aligned with `out`/`stats` via a NaN
                                # placeholder rather than shifting later entries
                                sanda_score.append(np.nan)
                                mean_rmsd.append(np.nan)
                                weighted_rmsd.append(np.nan)
                                affinity.append(np.nan)
                            else:
                                lb_rmsd = float(stat["vina_rmsd_lb"])
                                ub_rmsd = float(stat["vina_rmsd_ub"])
                                aff = float(stat["vina_affinity"])
                                w_rmsd = 0.7 * lb_rmsd + 0.3 * ub_rmsd

                                sanda_score.append(-np.abs(aff) / w_rmsd)

                                mean_rmsd.append((0.5 * lb_rmsd + 0.5 * ub_rmsd))
                                weighted_rmsd.append(w_rmsd)
                                affinity.append(aff)

                        # avg for 50/50, weights 70/30 lower/upper

                        if not np.isnan(sanda_score).all():
                            index = np.nanargmin(sanda_score)

                            confrmsd = mean_rmsd[index]
                            confscore = sanda_score[index]
                            confaff = affinity[index]

                            ligandmol = out[index]

                            templatepdb, templatepdbqt = update_vina(
                                pargs, success, fm, ligandmol, pargs.template_name
                            )
                            success += 1

                        else:
                            confrmsd = np.nan
                            confaff = np.nan
                            confscore = np.nan

                        print(
                            "%s %s %s %s %s"
                            % (
                                str(i).center(9),
                                fm.ljust(7),
                                str("%3.3f" % confrmsd).center(9),
                                str("%3.3f" % confaff).center(12),
                                str("%3.3f" % confscore).center(9),
                            )
                        )
                        pb.update(pb1, advance=1)

                    # export final output
                    final_print(templatepdb, dockmethod=pargs.dockmethod)

                    return pargs, templatepdb
                
                else:
                    if pargs.template_name:
                        templatepdb = pargs.template_name
                    else:
                        print("No template specified. Must have one for docking.")
                        exit()
                        
                    print(
                    "%s %s %s %s %s %s"
                        % (
                            "Interval".center(9),
                            "FM".center(7),
                            "gnina".center(12),
                            "gnina".center(12),
                            "Affinity".center(12),
                            "Intramol".center(12),
                        )
                    )
                    print(
                    "%s %s %s %s %s %s"
                        % (
                            " ".center(9),
                            " ".center(7),
                            "pose score".center(12),
                            "affinity".center(12),
                            "(kcal/mol)".center(12),
                            "(kcal/mol)".center(12),
                        )
                    )

                    pb1 = pb.add_task(
                        "[steel_blue1]Precomplexation", total=len(serialrecipe)
                    )

                    if not os.path.isdir("tmpgnina"):
                        os.mkdir("tmpgnina")
                    else:
                        subprocess.run(
                            f"rm tmpgnina/*",
                            shell=True,
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.STDOUT,
                        )

                    logdir = "tmpgnina"
                    logname = "gnina.log"
                    for i, fm in enumerate(serialrecipe):

                        ligand =f'tmpgnina/{fm}.pdb'

                        subprocess.run(
                            f"{shlex.quote(pargs.obabel)} -isdf {shlex.quote(sdfdir)}/{shlex.quote(fm)}.sdf -opdb -O {shlex.quote(ligand)}",
                            shell=True,
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.STDOUT,
                        )

                        ligandmol, stats = gnina(log=logname, 
                                    step = i,
                                    pargs=pargs, 
                                    template=templatepdb, 
                                    fm=ligand, 
                                    vina_args=vina_args)

                        templatepdb = update_gnina(
                            pargs, success, fm, ligandmol, basename=pargs.template_name
                        )

                        success += 1

                        print(
                            "%s %s %s %s %s %s"
                            % (
                                str(i).center(9),
                                fm.ljust(7),
                                str("%3.3f" % stats['posescore']).center(12),
                                str("%3.3f" % stats['cnnaffinity']).center(12),
                                str("%3.3f" % stats['affinity']).center(12),
                                str("%3.3f" % stats['intramol']).rjust(9),
                            )
                        )
                        pb.update(pb1, advance=1)

                    # export final output

                    final_print(templatepdb, dockmethod=pargs.dockmethod)

                    return pargs, templatepdb

            else:
                print("Docking requires FMs to be specified")

                generate_config(pargs.template_name, pargs)

                exit()

        ################################################################
        #
        # Main
        #
        ################################################################
        if args.screen:
            pargs = sanitize_inputs(args)
            
            if pargs.template_name:
                templatepdb = pargs.template_name
                templatepdbqt = preprocessing(args, receptor=templatepdb)
            else:
                print("No template specified. Must have one for screening.")
                exit()

            vina_args = generate_config(pargs.template_name, pargs)

            if pargs.dockmethod == 'vina':
                print(
                    "%s %s %s %s"
                    % (
                        "FM".ljust(12),
                        "Mean RMSD".center(12),
                        "Affinity".center(9),
                        "Score".center(9),
                    )
                )
                print(
                    "%s %s %s %s"
                    % (
                        "".ljust(12),
                        "(A)".center(12),
                        "(kcal/mol)".center(9),
                        "".center(9),
                    )
                )
            else:
                print(
                "%s %s %s %s %s"
                    % (
                        "FM".ljust(12),
                        "gnina".center(12),
                        "gnina".center(12),
                        "Affinity".center(12),
                        "Intramol".center(12),
                    )
                )
                print(
                "%s %s %s %s %s"
                    % (
                        " ".ljust(12),
                        "pose score".center(12),
                        "affinity".center(12),
                        "(kcal/mol)".center(12),
                        "(kcal/mol)".center(12),
                    )
                )

            fm2sm = fm2smiles()

            if pargs.fms:
                fmlist = {}
                for fm in pargs.fms:
                    fmlist[fm] = fm2sm[fm]
            else:
                fmlist = fm2sm

            fm_scores = []
            fm_rmsd = []
            fm_aff = []

            docked_mols = []

            for i, fm in enumerate(fmlist):
                
                if pargs.dockmethod == 'vina':
                    
                    logdir = "tmpvina"
                    logname = "vina.log"

                    ligand =f'tmpvina/{fm}.pdbqt'

                    subprocess.run(
                        f"{shlex.quote(pargs.obabel)} -isdf {shlex.quote(sdfdir)}/{shlex.quote(fm)}.sdf -opdbqt -O {shlex.quote(ligand)}",
                        shell=True,
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.STDOUT,
                    )

                    out, stats = autodock_vina(log=logname, 
                                    step = i,
                                    pargs=pargs, 
                                    template=templatepdbqt, 
                                    fm=ligand, 
                                    vina_args=vina_args)
                    
                    # take the lowest score, most negative (highest absolute value) affitinity and lowest mean rmsd
                    affinity = []
                    mean_rmsd = []
                    weighted_rmsd = []
                    sanda_score = []
                    for j, stat in enumerate(stats):
                        if j == 0:
                            # excluded (self-referential rmsd), but keep the
                            # index aligned with `out`/`stats` via a NaN
                            # placeholder rather than shifting later entries
                            sanda_score.append(np.nan)
                            mean_rmsd.append(np.nan)
                            weighted_rmsd.append(np.nan)
                            affinity.append(np.nan)
                        else:
                            lb_rmsd = float(stat["vina_rmsd_lb"])
                            ub_rmsd = float(stat["vina_rmsd_ub"])
                            aff = float(stat["vina_affinity"])
                            w_rmsd = 0.7 * lb_rmsd + 0.3 * ub_rmsd

                            sanda_score.append(-np.abs(aff) / w_rmsd)

                            mean_rmsd.append((0.5 * lb_rmsd + 0.5 * ub_rmsd))
                            weighted_rmsd.append(w_rmsd)
                            affinity.append(aff)
                    
                    # avg for 50/50, weights 70/30 lower/upper

                    if not np.isnan(sanda_score).all():
                        index = np.nanargmin(sanda_score)

                        confrmsd = mean_rmsd[index]
                        confaff = affinity[index]
                        confscore = sanda_score[index]

                        ligandmol = out[index]
                        ligandmol['name'] = fm
                        docked_mols.append(ligandmol)

                        fm_scores.append(confscore)
                        fm_aff.append(confaff)
                        fm_rmsd.append(confrmsd)
                    else:
                        confrmsd = np.nan
                        confaff = np.nan
                        confscore = np.nan

                        fm_scores.append(confscore)
                        fm_aff.append(confaff)
                        fm_rmsd.append(confrmsd)

                    print(
                        "%s %s %s %s"
                        % (
                            fm.ljust(12),
                            str("%3.3f" % confrmsd).center(12),
                            str("%3.3f" % confaff).center(9),
                            str("%3.3f" % confscore).center(9),
                        )
                    )
                else:

                    if not os.path.isdir("tmpgnina"):
                        os.mkdir("tmpgnina")
                    else:
                        subprocess.run(
                            f"rm tmpgnina/*",
                            shell=True,
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.STDOUT,
                        )
                    
                    logdir = "tmpgnina"
                    logname = "gnina.log"

                    ligand =f'tmpgnina/{fm}.pdb'

                    subprocess.run(
                        f"{shlex.quote(pargs.obabel)} -isdf {shlex.quote(sdfdir)}/{shlex.quote(fm)}.sdf -opdb -O {shlex.quote(ligand)}",
                        shell=True,
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.STDOUT,
                    )

                    ligandmol, stats = gnina(log=logname, 
                                step = i,
                                pargs=pargs, 
                                template=templatepdbqt, 
                                fm=ligand, 
                                vina_args=vina_args)
                    
                    ligandmol['name'] = fm
                    docked_mols.append(ligandmol)

                    print(
                        "%s %s %s %s %s"
                        % (
                            fm.ljust(12),
                            str("%3.3f" % stats['posescore']).center(12),
                            str("%3.3f" % stats['cnnaffinity']).center(12),
                            str("%3.3f" % stats['affinity']).center(12),
                            str("%3.3f" % stats['intramol']).center(12),
                        )
                    )

                    fm_scores.append(stats['posescore'])
                    fm_aff.append(stats['affinity'])
            
            if pargs.dockmethod == 'vina':
                zipped = zip(fmlist, fm_scores, fm_rmsd, fm_aff)
                nansrt = sorted(zipped, key=lambda t: t[1], reverse=False)
                srt = []
                for nn in nansrt:
                    if str(nn[1]) == 'nan':
                        continue
                    else:
                        srt.append(nn)

                print("\n----------------------------------\n")
                print("Top scoring functional monomers\n")
                print("----------------------------------\n")
                print(
                    "%s %s %s %s"
                    % (
                        "FM".ljust(12),
                        "Mean RMSD".center(12),
                        "Affinity".center(9),
                        "Score".center(9),
                    )
                )
                print(
                    "%s %s %s %s"
                    % (
                        "".ljust(12),
                        "(A)".center(12),
                        "(kcal/mol)".center(9),
                        "".center(9),
                    )
                )

                if len(srt) > 10:
                    srtlen = 10
                else:
                    srtlen = len(srt)

                for i in range(0, srtlen):
                    print(
                        "%s %s %s %s"
                        % (
                            srt[i][0].ljust(12),
                            str("%3.3f" % srt[i][2]).center(12),
                            str("%3.3f" % srt[i][3]).center(9),
                            str("%3.3f" % srt[i][1]).center(9),
                        )
                    )
            else:
                zipped = zip(fmlist, fm_scores, fm_aff)
                nansrt = sorted(zipped, key=lambda t: t[2], reverse=False)
                srt = []
                for nn in nansrt:
                    if str(nn[1]) == 'nan':
                        continue
                    else:
                        srt.append(nn)

                print("\n----------------------------------\n")
                print("Top scoring functional monomers\n")
                print("----------------------------------\n")
                print(
                    "%s %s %s"
                    % (
                        "FM".ljust(12),
                        "Score".center(9),
                        "Affinity".center(9),
                    )
                )

                if len(srt) > 10:
                    srtlen = 10
                else:
                    srtlen = len(srt)

                for i in range(0, srtlen):
                    print(
                        "%s %s %s"
                        % (
                            srt[i][0].ljust(12),
                            str("%3.3f" % srt[i][1]).center(9),
                            str("%3.3f" % srt[i][2]).center(9),
                        )
                    )

            print_molecule(docked_mols, "screened_fm_locations.pdb")