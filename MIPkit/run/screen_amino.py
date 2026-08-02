import os, subprocess, sys
import shlex
import numpy as np
from rich.progress import Progress, TimeElapsedColumn, BarColumn, TextColumn

# our work
from MIPkit.constants.constants import (
    fm2smiles,
    decode_fms,
    encode_fms,
    sulfur_exception,
)
from MIPkit.utils.parser import sanitize_inputs
from MIPkit.utils.read_pdb import (
    read_pdb,
)
from MIPkit.run.precomplex import precomplex
from MIPkit.run.run_acpype import acpype_do
from MIPkit.run.run_gmx import gmx_do
from MIPkit.utils.print_molecule import print_molecule
from MIPkit.utils.utils import generic_top, generate_recipe, generate_config
from MIPkit.utils.read_yaml import read_yaml

from MIPkit.run.run_obabel import obabel_do, obabel_combine
from MIPkit.utils.identify import identify

from MIPkit.vina.vina import autodock_vina, preprocessing, update_vina
from MIPkit.vina.gnina import gnina, update_gnina

################################################################
#
# Screen fms against individual using a list of FMS
#
# This is not accessable nor completed
#
################################################################
def local_dock(pargs, args, fm):

    sdfdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../sdfs/"))
    sdfdir_q = shlex.quote(sdfdir)
    fm_q = shlex.quote(fm)
    obabel_q = shlex.quote(pargs.obabel)

    vina_args = generate_config(pargs.template_name, args)

    if pargs.dockmethod == 'vina':
        
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
        
        logdir = "tmpvina"
        logname = "vina.log"

        ligand =f'tmpvina/{fm}.pdbqt'

        subprocess.run(
            f"{obabel_q} -isdf {sdfdir_q}/{fm_q}.sdf -opdbqt -O {shlex.quote(ligand)}",
            shell=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT,
        )

        out, stats = autodock_vina(log=logname, 
                    step = 0,
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
                # excluded (self-referential rmsd), but keep the index
                # aligned with `out`/`stats` via a NaN placeholder rather
                # than shifting later entries
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

            print(
                "%s %s %s %s %s"
                % (
                    str(0).center(9),
                    fm.ljust(7),
                    str("%3.3f" % confrmsd).center(9),
                    str("%3.3f" % confaff).center(12),
                    str("%3.3f" % confscore).center(9),
                )
            )

            ligandmol['name'] = fm
            
            return ligandmol

        else:
            print("Docking failed. Exiting ...")
            exit()
        
    
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

        logdir = "tmpgnina"
        logname = "gnina.log"

        ligand =f'tmpgnina/{fm}.pdb'

        subprocess.run(
            f"{obabel_q} -isdf {sdfdir_q}/{fm_q}.sdf -opdb -O {shlex.quote(ligand)}",
            shell=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT,
        )

        ligandmol, stats = gnina(log=logname, 
                    step = 0,
                    pargs=pargs, 
                    template=templatepdb, 
                    fm=ligand, 
                    vina_args=vina_args)

        print(
            "%s %s %s %s %s %s"
            % (
                str(0).center(9),
                fm.ljust(7),
                str("%3.3f" % stats['posescore']).center(12),
                str("%3.3f" % stats['cnnaffinity']).center(12),
                str("%3.3f" % stats['affinity']).center(12),
                str("%3.3f" % stats['intramol']).rjust(9),
            )
        )

        # export final output

        ligandmol['name'] = fm
        
        return ligandmol

def amino_screen(args):
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
        # Docking individual amino acids
        #
        ################################################################

        if not os.path.isdir("tmpfms"):
            os.mkdir("tmpfms")
        else:
            subprocess.run(
                f"rm tmpfms/*",
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )

        # --------------------------------------------------------------
        defm = decode_fms()
        enfm = encode_fms()
        fmsmiles = fm2smiles()
        sulfurs = sulfur_exception()

        pargs = sanitize_inputs(args)

        for fm in pargs.fms:
            if fm not in fmsmiles:
                print(f"Functional Monomer : {fm} is not in the library.\nCheck for misspellings, or add it to MIPkit/constants/fm-list.yaml.")
                sys.exit()

        workdir = pargs.workdir
        compdir = pargs.compdir
        workdir_q = shlex.quote(workdir)
        compdir_q = shlex.quote(compdir)

        if pargs.dockmethod=='vina':
            if not os.path.isdir("tmpvina"):
                    os.mkdir("tmpvina")
            else:
                subprocess.run(
                    f"rm tmpvina/*",
                    shell=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.STDOUT,
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

        screening_fms = pargs.fms

        screening_fms_docked = []
        
        for fm in screening_fms:
            fm_q = shlex.quote(fm)
            
            subprocess.run(f"rm {compdir_q}/*",
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )
            subprocess.run(f"rm {workdir_q}/*",
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )

            pargs.fms = [fm]

            docked_fm = local_dock(pargs, args, fm)
            screening_fms_docked.append(docked_fm)

            print_molecule([docked_fm], f"tmpfms/{fm}.pdb")
            
            # For FM identification and system setup
            pargs.config = "generated_config.yaml"
            pargs.complex = f"tmpfms/{fm}.pdb"

            # Obabel
            logname = f"obabel-{fm}.log"
            out_fm = obabel_do("screen", 0, [docked_fm], pargs, logname)

            # Acpype
            logname = f"acpype-{fm}.log"
            acpype_do("screen", 0, [docked_fm], out_fm, pargs, logname)

            # Combined
            logname = f"combination-{fm}.log"
            obabel_combine("interact", 0, out_fm, pargs, logname)

            # GMX
            logname = f"gmx-{fm}.log"
            gmx_run = gmx_do()
            gmx_run.params("interact", 0, 0, pargs, logname)
            gmx_run.run()

            # ENERGY
            subprocess.run(f"cp {workdir_q}template-{fm_q}_lj.xvg tmpfms/template-{fm_q}_lj.xvg", shell=True)
            subprocess.run(f"cp {workdir_q}template-{fm_q}_coul.xvg tmpfms/template-{fm_q}_coul.xvg", shell=True)

            # H-BONDS
            subprocess.run(f"cp {workdir_q}template-{fm_q}_hbond_num.xvg tmpfms/template-{fm_q}_hbond_num.xvg", shell=True)
            subprocess.run(f"cp {workdir_q}template-{fm_q}_hbond_dist.xvg tmpfms/template-{fm_q}_hbond_dist.xvg", shell=True)

            # RDF
            subprocess.run(f"cp {workdir_q}template-{fm_q}_rdf.xvg tmpfms/template-{fm_q}_rdf.xvg", shell=True)
        
        print_molecule(screening_fms_docked, f"screen_fms_docked_locations.pdb")


