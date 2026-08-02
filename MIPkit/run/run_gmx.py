import subprocess, os, sys, time, copy
import shlex
from MIPkit.utils.utils import generic_top
from MIPkit.utils.read_pdb import read_pdb
from MIPkit.utils.print_molecule import print_molecule
from MIPkit.utils.read_ndx import clean_ndx, merge_ndx_groups
from MIPkit.constants.constants import encode_fms
from MIPkit.run.run_obabel import obabel_solvent


class gmx_do():
	def __init__(self):
		self.run_type = None
		self.i = None
		self.maxi = None
		self.pargs = None
		self.logname = None
		self.rolled = False

	def params(self, run_type, i, maxi, pargs, logname):
		self.run_type = run_type
		self.i = i
		self.maxi = maxi
		self.pargs = copy.deepcopy(pargs)
		self.logname = logname
		self.init_dt = copy.deepcopy(self.pargs.dt['rxn'])
	
	def run(self):

		run_type = self.run_type
		i = self.i
		maxi = self.maxi
		pargs = self.pargs
		logname = self.logname

		if run_type == "min":
			prefix = "min"
		elif run_type == "rxn":
			prefix = f"step{i}"
		elif run_type == "interact":
			prefix = "interact"
		else:
			print("Incorrect run type, this should not happen.")
			exit()

		prefix_q = shlex.quote(prefix)

		enfm = encode_fms()

		if pargs.explicit:
			include_initiators = True
		else:
			include_initiators = False

		# directories
		amberdir = pargs.amberdir
		compdir = pargs.compdir
		gmxdir = pargs.gmxdir
		logdir = pargs.logdir
		solventdir = pargs.solventdir
		workdir = pargs.workdir

		# gmx
		gmx = pargs.gmx
		ncpu = pargs.ncpu
		gpuid = pargs.gpu_id
		bonded_flag = "-bonded auto"
		# 1 thread-MPI rank per GPU is recommended for single-GPU runs;
		# without a GPU, let -nt alone decide the thread-MPI/OpenMP split.
		thread_flags = f"-ntmpi 1 -ntomp {ncpu}" if gpuid else f"-nt {ncpu}"

		box_str = pargs.box
		gmxtype = pargs.gmxtype
		solvent = pargs.solvent

		# quoted variants for safe interpolation into shell=True commands
		compdir_q = shlex.quote(compdir)
		gmxdir_q = shlex.quote(gmxdir)
		solventdir_q = shlex.quote(solventdir)
		workdir_q = shlex.quote(workdir)
		gmx_q = shlex.quote(gmx)

		# ===========================================================
		#
		# Main Do
		#
		# ===========================================================

		with open(os.path.join(logdir, f"{logname}"), "w+") as f:

			if pargs.template_name:
				gtop = generic_top(
					pargs,
					'#include "template-clean.itp"',
					"Template".ljust(16) + "1",
					'#include "RXN_merged.itp"',
					"CLX".ljust(16) + "1",
					"template_posres.itp",
					initiator=include_initiators,
					solvents=pargs.solvent
				)
			else:
				gtop = generic_top(
					pargs,
					"",
					"",
					'#include "RXN_merged.itp"',
					"CLX".ljust(16) + "1",
					initiator=include_initiators,
					solvents=pargs.solvent
				)

			with open(f"{workdir}rxn_topol.top", "w+") as _top:
				_top.write(gtop)

			subprocess.run(
				f'cp {workdir_q}rxn_topol.top {compdir_q}topol.top',
				shell=True,
				stdout=f,
				stderr=f,
			)

			if pargs.shell:
				shell = pargs.shell
				shell_q = shlex.quote(shell)
				subprocess.run(
					f"bash {shell_q} {workdir_q}{prefix_q}-fullcomplex.pdb {i} ",
					shell=True,
				)  # , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

			else:
				maxwarn = "-maxwarn 1"
				if not pargs.template_name:
					# turn on max warning to ignore posres
					maxwarn = "-maxwarn 2"
				# preparation
				subprocess.run(
					f'{gmx_q} editconf -f {workdir_q}{prefix_q}-fullcomplex.pdb  -o {workdir_q}rxn_box.gro {box_str} -c',
					shell=True,
					stdout=f,
					stderr=f,
				)

				if (pargs.explicit and i == 0 and not pargs.min) or (pargs.explicit and pargs.min and run_type == "min"):
					# grab initiators and count to add by insert molecules
					initiatordir = pargs.initiatordir
					initiatordir_q = shlex.quote(initiatordir)
					initiators = pargs.initiator
					for initiator in initiators:
						init_count = initiators[initiator]

						subprocess.run(
							f'{gmx_q} insert-molecules -f {workdir_q}rxn_box.gro -nmol {init_count} -ci {initiatordir_q}/{shlex.quote(initiator.lower())}.pdb -o {workdir_q}rxn_box.gro ',
							shell=True,
							stdout=f,
							stderr=f,
						)

					# grab initiators and count to add to topol file
					tmp_top = ""
					initiators = pargs.initiator
					for __i, initiator in enumerate(initiators):
						nl = ""
						if __i < len(initiators)-1:
							nl = "\n"
						tmp_top += initiator.upper().ljust(16) + str(initiators[initiator]) + nl

					subprocess.run(
						f'echo {shlex.quote(tmp_top)}  | tee -a {workdir_q}rxn_topol.top',
						shell=True,
						stdout=f,
						stderr=f,
					)

				elif pargs.explicit and run_type == "rxn" and (i > 0 or pargs.min):
					# grab initiators and count to add to topol file (updated counts after reactions)
					tmp_top = ""
					initiators = pargs.initiator
					for __i, initiator in enumerate(initiators):
						nl = ""
						if __i < len(initiators)-1:
							nl = "\n"
						tmp_top += initiator.upper().ljust(16) + str(initiators[initiator]) + nl

					subprocess.run(
						f'echo {shlex.quote(tmp_top)}  | tee -a {workdir_q}rxn_topol.top',
						shell=True,
						stdout=f,
						stderr=f,
					)

				'''
				Finalize System
				- 1 - Minimize to get tpr
				- 2 - Solvent to get solvent locations for ions
				- 3 - Replace ions
				- 4 - Reminimize
				'''

				# Vacuum minimization
				gmx_make_mdp(gmxdir,workdir,'min', pargs)

				subprocess.run(
					f'{gmx_q} grompp -f {workdir_q}/fm_min.mdp -c {workdir_q}rxn_box.gro -r {workdir_q}rxn_box.gro -p {workdir_q}rxn_topol.top -o {workdir_q}rxn_min.tpr -maxwarn 2',
					shell=True,
					stdout=f,
					stderr=f,
				)
				output = subprocess.run(
					f'{gmx_q} mdrun {thread_flags} {gpuid} {bonded_flag} -v -deffnm {workdir_q}rxn_min',
					shell=True,
					stdout=f,
					stderr=subprocess.PIPE,
					text=True
				)
				if self.gmx_error_check(output,"rxn_min",177):
					return

				'''
				NVT Run
				- rxn, min
				- gpu, cpu
				'''
				gmx_make_mdp(gmxdir,workdir,'nvt_vacuum', pargs)

				subprocess.run(
						f'{gmx_q} grompp -f {workdir_q}/fm_nvt_vacuum.mdp -c {workdir_q}rxn_min.gro -r {workdir_q}rxn_min.gro -p {workdir_q}rxn_topol.top -o {workdir_q}rxn_min.tpr  {maxwarn}',
						shell=True,
						stdout=f,
						stderr=f,
					)
				output = subprocess.run(
					f'{gmx_q} mdrun {thread_flags} {gpuid} {bonded_flag} -v -deffnm {workdir_q}rxn_min ',
					shell=True,
					stdout=subprocess.PIPE,
					stderr=subprocess.PIPE,
					text = True
				)
				if self.gmx_error_check(output,"rxn_vacuum_nvt",213):
					return

				# convert to 

				# if run_type == 'min' or run_type == 'interact' or (run_type == 'rxn' and i==0 and not pargs.min):
				if pargs.solvent and pargs.solvent != "spc":
					subprocess.run(
						f'{gmx_q} solvate -cp {workdir_q}rxn_box.gro -cs {solventdir_q}/{shlex.quote(solvent)}.gro -p {workdir_q}rxn_topol.top -o {workdir_q}rxn_solv_box.gro -radius 1.5',
						shell=True,
						stdout=f,
						stderr=f,
					)
				else:
					subprocess.run(
						f'{gmx_q} solvate -cp {workdir_q}rxn_box.gro -cs "spc216.gro" -p {workdir_q}rxn_topol.top -o {workdir_q}rxn_solv_box.gro -radius 1.5',
						shell=True,
						stdout=f,
						stderr=f,
					)
				# else:
				# 	subprocess.run(
				# 		f'{gmx} editconf -f "{workdir}rxn_min.gro" -o "{workdir}{prefix}-fullcomplex.pdb"',
				# 		shell=True,
				# 		stdout=f,
				# 		stderr=f,
				# 	)

				# 	tmp_top = ""
				# 	if pargs.solvent:
				# 		tmp_top += pargs.solvent.upper().ljust(16) + str(solvent_count)
				# 	else:
				# 		tmp_top += 'SOL'.ljust(16) + str(solvent_count)

				# 	subprocess.run(
				# 		f'echo "{tmp_top}"  | tee -a "{workdir}rxn_topol.top"',
				# 		shell=True,
				# 		stdout=f,
				# 		stderr=f,
				# 	)

				# 	obabel_solvent("rxn",i, pargs, 'rxn_solv_box', 'pdb', 'solvent_add.log')

				# 	subprocess.run(
				# 		f'{gmx} editconf -f "{workdir}rxn_solv_box.pdb"  -o "{workdir}rxn_solv_box.gro"  {box_str} -c',
				# 		shell=True,
				# 		stdout=f,
				# 		stderr=f,
				# 	)

				# generate tpr
				subprocess.run(
					f'{gmx_q} grompp -f {gmxdir_q}/fm_min.mdp -c {workdir_q}rxn_solv_box.gro -r {workdir_q}rxn_solv_box.gro -p {workdir_q}rxn_topol.top -o {workdir_q}rxn_min.tpr -maxwarn 2',
					shell=True,
					stdout=f,
					stderr=f,
				)

				if pargs.solvent == "spc" or not pargs.solvent:
					subprocess.run(
						f'echo "SOL" | {gmx_q} genion -s {workdir_q}rxn_min.tpr -p {workdir_q}rxn_topol.top -o {workdir_q}rxn_solv_box.gro -pname NA -nname CL -neutral',
						shell=True,
						stdout=f,
						stderr=f,
					)
				else:
					uSOL = pargs.solvent.upper()
					subprocess.run(
						f'echo {shlex.quote(uSOL)} | {gmx_q} genion -s {workdir_q}rxn_min.tpr -p {workdir_q}rxn_topol.top -o {workdir_q}rxn_solv_box.gro -pname NA -nname CL -neutral',
						shell=True,
						stdout=f,
						stderr=f,
					)
				'''
				Minimize Run
				- rxn, min
				- gpu, cpu
				
				- report this out to log; poorly equilibrated solvents can throw an exception error
				- also, during dev, if centos goes to sleep cuda might not restart on wake
				'''

				gmx_make_mdp(gmxdir,workdir,'min', pargs)

				subprocess.run(
					f'{gmx_q} grompp -f {workdir_q}/fm_min.mdp -c {workdir_q}rxn_solv_box.gro -r {workdir_q}rxn_solv_box.gro -p {workdir_q}rxn_topol.top -o {workdir_q}rxn_min.tpr  {maxwarn}',
					shell=True,
					stdout=f,
					stderr=f,
				)
				output = subprocess.run(
					f'{gmx_q} mdrun {thread_flags} {gpuid} {bonded_flag} -v -deffnm {workdir_q}rxn_min',
					shell=True,
					stdout=f,
					stderr=subprocess.PIPE,
					text=True
				)
				if self.gmx_error_check(output,"rxn_min",310):
					return

				# for s in [0,1,2]:
				# 	gmx_make_mdp(gmxdir,workdir,f'nvt_stage{s}', pargs)

				# 	if s == 0:
				# 		subprocess.run(
				# 				f'{gmx} grompp -f "{workdir}/fm_nvt_stage{s}.mdp" -c "{workdir}rxn_min.tpr" -r "{workdir}rxn_min.tpr" -p "{workdir}rxn_topol.top" -o "{workdir}rxn_nvt.tpr"  {maxwarn}',
				# 				shell=True,
				# 				stdout=f,
				# 				stderr=f,
				# 			)
				# 	else:
				# 		subprocess.run(
				# 				f'{gmx} grompp -f "{workdir}/fm_nvt_stage{s}.mdp" -c "{workdir}rxn_nvt.tpr" -r "{workdir}rxn_nvt.tpr" -p "{workdir}rxn_topol.top" -o "{workdir}rxn_nvt.tpr"  {maxwarn}',
				# 				shell=True,
				# 				stdout=f,
				# 				stderr=f,
				# 			)
				# 	output = subprocess.run(
				# 		f'{gmx} mdrun -nt {ncpu}  {gpuid} -pin auto -v -deffnm "{workdir}rxn_nvt" ',
				# 		shell=True,
				# 		stdout=subprocess.PIPE,
				# 		stderr=subprocess.PIPE,
				# 		text = True
				# 	)
				# 	self.gmx_error_check(output,f"rxn_nvt_stage{s}",213)


				'''
				NVT Run
				- rxn, min
				- gpu, cpu
				'''
				if run_type == 'min':
					gmx_make_mdp(gmxdir,workdir,'nvt', pargs)
					
					subprocess.run(
						f'{gmx_q} grompp -f {workdir_q}/fm_nvt.mdp -c {workdir_q}rxn_min.gro -r {workdir_q}rxn_min.gro -p {workdir_q}rxn_topol.top -o {workdir_q}rxn_nvt.tpr  {maxwarn}',
						shell=True,
						stdout=f,
						stderr=f,
					)
				else:
					gmx_make_mdp(gmxdir,workdir,'nvt', pargs)

					subprocess.run(
						f'{gmx_q} grompp -f {workdir_q}/fm_nvt.mdp -c {workdir_q}rxn_min.gro -r {workdir_q}rxn_min.gro -p {workdir_q}rxn_topol.top -o {workdir_q}rxn_nvt.tpr  {maxwarn}',
						shell=True,
						stdout=f,
						stderr=f,
					)
				output = subprocess.run(
					f'{gmx_q} mdrun {thread_flags} {gpuid} {bonded_flag} -v -deffnm {workdir_q}rxn_nvt ',
					shell=True,
					stdout=subprocess.PIPE,
					stderr=subprocess.PIPE,
					text = True
				)
				if self.gmx_error_check(output,"rxn_nvt",213):
					return

				'''
				NPT Run
				- rxn, min
				- gpu, cpu
				'''
				if run_type == 'min':
					gmx_make_mdp(gmxdir,workdir,'npt', pargs)

					subprocess.run(
						f'{gmx_q} grompp -f {workdir_q}/fm_npt.mdp -c {workdir_q}rxn_nvt.gro -r {workdir_q}rxn_nvt.gro -t {workdir_q}rxn_nvt.cpt -p {workdir_q}rxn_topol.top -o {workdir_q}rxn_npt.tpr  {maxwarn}',
						shell=True,
						stdout=f,
						stderr=f,
					)
				else:
					gmx_make_mdp(gmxdir,workdir,'npt', pargs)

					subprocess.run(
						f'{gmx_q} grompp -f {workdir_q}/fm_npt.mdp -c {workdir_q}rxn_nvt.gro -r {workdir_q}rxn_nvt.gro -t {workdir_q}rxn_nvt.cpt -p {workdir_q}rxn_topol.top -o {workdir_q}rxn_npt.tpr  {maxwarn}',
						shell=True,
						stdout=f,
						stderr=f,
					)
				output = subprocess.run(
					f'{gmx_q} mdrun {thread_flags} {gpuid} {bonded_flag} -v -deffnm {workdir_q}rxn_npt ',
					shell=True,
					stdout=subprocess.PIPE,
					stderr=subprocess.PIPE,
					text=True,
				)
				if self.gmx_error_check(output, "rxn_npt", 226):
					return

				'''
				Production Run
				- rxn, min
				- gpu, cpu
				'''
				if run_type == "rxn":
					gmx_make_mdp(gmxdir,workdir,'rxn', pargs)

					subprocess.run(
						f'{gmx_q} grompp -f {workdir_q}/fm_rxn.mdp -c {workdir_q}rxn_npt.gro -r {workdir_q}rxn_npt.gro -t {workdir_q}rxn_npt.cpt -p {workdir_q}rxn_topol.top -o {workdir_q}rxn_prod.tpr -maxwarn 1',
						shell=True,
						stdout=f,
						stderr=f,
					)
				elif run_type == "min":
					gmx_make_mdp(gmxdir,workdir,'min_run', pargs)

					subprocess.run(
						f'{gmx_q} grompp -f {workdir_q}/fm_min_run.mdp -c {workdir_q}rxn_npt.gro -r {workdir_q}rxn_npt.gro -t {workdir_q}rxn_npt.cpt -p {workdir_q}rxn_topol.top -o {workdir_q}rxn_prod.tpr -maxwarn 2',
						shell=True,
						stdout=f,
						stderr=f,
					)

				if not run_type == "interact":
					output = subprocess.run(
						f'{gmx_q} mdrun {thread_flags} {gpuid} {bonded_flag} -v -deffnm {workdir_q}rxn_prod ',
						shell=True,
						stdout=subprocess.PIPE,
						stderr=subprocess.PIPE,
						text = True,
					)
					if self.gmx_error_check(output,"rxn_prod",260):
						return


				# Copy
				if run_type == "rxn" or run_type == "min":
					# _prod.gro is not being produced? just export npt instead
					subprocess.run(
						f'cp {workdir_q}rxn_min.gro {compdir_q}{prefix_q}_prod.gro ',
						shell=True,
						stdout=subprocess.DEVNULL,
						stderr=subprocess.DEVNULL,
					)
					subprocess.run(
						f'cp {workdir_q}rxn_prod.xtc {compdir_q}{prefix_q}_prod.xtc ',
						shell=True,
						stdout=subprocess.DEVNULL,
						stderr=subprocess.DEVNULL,
					)
					subprocess.run(
						f'cp {workdir_q}rxn_prod.tpr {compdir_q}{prefix_q}_prod.tpr ',
						shell=True,
						stdout=subprocess.DEVNULL,
						stderr=subprocess.DEVNULL,
					)
					subprocess.run(
						f'cp {workdir_q}rxn_prod.tpr {compdir_q}{prefix_q}_prod.edr ',
						shell=True,
						stdout=subprocess.DEVNULL,
						stderr=subprocess.DEVNULL,
					)
					subprocess.run(
						f'cp {workdir_q}rxn_topol.top {compdir_q}{prefix_q}_topol.top ',
						shell=True,
						stdout=subprocess.DEVNULL,
						stderr=subprocess.DEVNULL,
					)
				else:
					subprocess.run(
						f'cp {workdir_q}rxn_topol.top {compdir_q}{prefix_q}_topol.top ',
						shell=True,
						stdout=subprocess.DEVNULL,
						stderr=subprocess.DEVNULL,
					)
				subprocess.run(
						f'cp {workdir_q}rxn_topol.top {compdir_q}topol.top ',
						shell=True,
						stdout=subprocess.DEVNULL,
						stderr=subprocess.DEVNULL,
					)

			'''
			Build ndxs and export for next timestep
			- clean_ndx can take both inclusions and exclusions
			- merge_nd_groups to combine all functional monomer groups
			'''
			if not run_type == "interact":
				output = subprocess.run(
					f'(echo 0 ; echo "q") | {gmx_q} make_ndx -f {workdir_q}rxn_prod.tpr -o {workdir_q}index.ndx', # "System"
					shell=True,
					stdout=subprocess.PIPE,
					stderr=subprocess.PIPE,
					text = True,
				)
				if self.gmx_error_check(output,"rxn_prod_indexing",391):
					return
			else:
				output = subprocess.run(
					f'(echo 0 ; echo "q") | {gmx_q} make_ndx -f {workdir_q}rxn_npt.tpr -o {workdir_q}index.ndx', # "System"
					shell=True,
					stdout=subprocess.PIPE,
					stderr=subprocess.PIPE,
					text = True,
				)
				if self.gmx_error_check(output,"rxn_prod_indexing",391):
					return

			if os.path.isfile(f"{workdir}index.ndx"):
				clean_ndx(f"{workdir}index.ndx", f"{workdir}fms.ndx", include=pargs.fms, exclude=False)
				merge_ndx_groups(f"{workdir}fms.ndx",f"{workdir}fms-all.ndx")
			else:
				print('\nError: No index file created. Check the gmx step log but this is usually due to segfaults during minimization.\nLine 351 run_gmx')
				exit()

			if pargs.template or pargs.protein:

				_clean = read_pdb(f"{workdir}template.pdb", hydrogens=False, initiators=False)
				clean_res = []
				for mol in _clean:
					clean_res.append(mol['name'])

				clean_ndx(f"{workdir}index.ndx", f"{workdir}template.ndx", include=clean_res, exclude=False)
				template_groups = list(set(clean_res))

				for fm in pargs.fms:
					clean_res.append(fm)
				clean_ndx(f"{workdir}index.ndx", f"{workdir}template-fms.ndx", include=clean_res, exclude=False)
				eng_groups = list(set(clean_res))

			'''
			Parse information and separate into pdbs
			- protein and not solvent -- dont need to use indexes
			- template and not solvent -- need to use indexes for template and fms
			- protein and solvent -- dont need to use indexes
			- template and solvent -- need to use indexes for template and fms
			- solvent -- need to use indexes for fms
			- else -- ? the old method, have it print a warning and ask them to contact me
			'''
			if not run_type == "interact":

				if pargs.protein and not pargs.solvent:
					# center system
					subprocess.run(
						f'( echo "Protein" ; echo "Non-Water" ) | {gmx_q} trjconv -f {workdir_q}rxn_prod.xtc -s {workdir_q}rxn_prod.tpr -o {workdir_q}{prefix_q}-fullsystem.pdb -dump 1000000 -center -pbc mol -ur compact',
						shell=True,
						stdout=f,
						stderr=f,
					)
					if not os.path.exists(f"{workdir}{prefix}-fullsystem.pdb"):
						print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
						print(f'\nGMX Error : GROMACS failed to produce "{workdir}{prefix}-fullsystem.pdb".\nSee logs/gmx-step... for details.\n\n')
						print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
						exit()

					subprocess.run(
						f'{gmx_q} editconf -f {workdir_q}{prefix_q}-fullsystem.pdb -o {workdir_q}{prefix_q}-fullsystem.pdb ' , #-center 0 0 0',
						shell=True,
						stdout=f,
						stderr=f,
					)

					subprocess.run(
						f'( echo "Protein" ) | {gmx_q} editconf -f {workdir_q}{prefix_q}-fullsystem.pdb -o {workdir_q}{prefix_q}-protein-prod.pdb -ndef',
						shell=True,
						stdout=f,
						stderr=f,
					)
					subprocess.run(
						f'{gmx_q} editconf -f {workdir_q}{prefix_q}-fullsystem.pdb -n {workdir_q}fms-all.ndx -o {workdir_q}{prefix_q}-fullcomplex-prod.pdb -ndef',
						shell=True,
						stdout=f,
						stderr=f,
					)

					if not os.path.exists(f"{workdir}{prefix}-fullcomplex-prod.pdb"):
						print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
						print(f'\nGMX Error : GROMACS failed to produce "{workdir}{prefix}-fullcomplex-prod.pdb".\nSee logs/gmx-step... for details.\n\n')
						print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
						exit()

				elif pargs.template and not pargs.solvent:

					subprocess.run(
						f'( echo "Other" ; echo "Other" ) | {gmx_q} trjconv -f {workdir_q}rxn_prod.xtc -s {workdir_q}rxn_prod.tpr -o {workdir_q}{prefix_q}-fullsystem.pdb -dump 1000000 -center -pbc mol -ur compact',
						shell=True,
						stdout=f,
						stderr=f,
					)

					if not os.path.exists(f"{workdir}{prefix}-fullsystem.pdb"):
						print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
						print(f'\nGMX Error : GROMACS failed to produce "{workdir}{prefix}-fullsystem.pdb".\nSee logs/gmx-step... for details.\n\n')
						print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
						exit()

					subprocess.run(
						f'{gmx_q} editconf -f {workdir_q}{prefix_q}-fullsystem.pdb -o {workdir_q}{prefix_q}-fullsystem.pdb ' , #-center 0 0 0',
						shell=True,
						stdout=f,
						stderr=f,
					)

					output = subprocess.run(
						f'{gmx_q} editconf -f {workdir_q}{prefix_q}-fullsystem.pdb -n {workdir_q}template.ndx -o {workdir_q}{prefix_q}-template-prod.pdb -ndef',
						shell=True,
						stdout=f,
						stderr=subprocess.PIPE,
						text=True,
					)
					if self.gmx_error_check(output, "template_ndx", 491):
						return


					subprocess.run(
						f'{gmx_q} editconf -f {workdir_q}{prefix_q}-fullsystem.pdb -n {workdir_q}fms-all.ndx -o {workdir_q}{prefix_q}-fullcomplex-prod.pdb -ndef',
						shell=True,
						stdout=f,
						stderr=f,
					)

					if not os.path.exists(f"{workdir}{prefix}-fullcomplex-prod.pdb"):
						print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
						print(f'\nGMX Error : GROMACS failed to produce "{workdir}{prefix}-fullcomplex-prod.pdb".\nSee logs/gmx-step... for details.\n\n')
						print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
						exit()

				elif pargs.protein and pargs.solvent:

					subprocess.run(
						f'( echo "Protein" ; echo "Other" ) | {gmx_q} trjconv -f {workdir_q}rxn_prod.xtc -s {workdir_q}rxn_prod.tpr -o {workdir_q}{prefix_q}-fullsystem.pdb -dump 1000000 -center -pbc mol -ur compact',
						shell=True,
						stdout=f,
						stderr=f,
					)
					if not os.path.exists(f"{workdir}{prefix}-fullsystem.pdb"):
						print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
						print(f'\nGMX Error : GROMACS failed to produce "{workdir}{prefix}-fullsystem.pdb".\nSee logs/gmx-step... for details.\n\n')
						print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
						exit()

					subprocess.run(
						f'{gmx_q} editconf -f {workdir_q}{prefix_q}-fullsystem.pdb -o {workdir_q}{prefix_q}-fullsystem.pdb ' , #-center 0 0 0',
						shell=True,
						stdout=f,
						stderr=f,
					)

					subprocess.run(
						f'( echo "Protein" ) | {gmx_q} editconf -f {workdir_q}{prefix_q}-fullsystem.pdb -o {workdir_q}{prefix_q}-protein-prod.pdb -ndef',
						shell=True,
						stdout=f,
						stderr=f,
					)
					subprocess.run(
						f'{gmx_q} editconf -f {workdir_q}{prefix_q}-fullsystem.pdb -n {workdir_q}fms-all.ndx -o {workdir_q}{prefix_q}-fullcomplex-prod.pdb -ndef',
						shell=True,
						stdout=f,
						stderr=f,
					)

					if not os.path.exists(f"{workdir}{prefix}-fullcomplex-prod.pdb"):
						print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
						print(f'\nGMX Error : GROMACS failed to produce "{workdir}{prefix}-fullcomplex-prod.pdb".\nSee logs/gmx-step... for details.\n\n')
						print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
						exit()

				elif pargs.template and pargs.solvent:

					subprocess.run(
						f'( echo "Other" ; echo "Other" ) | {gmx_q} trjconv -f {workdir_q}rxn_prod.xtc -s {workdir_q}rxn_prod.tpr -o {workdir_q}{prefix_q}-fullsystem.pdb -dump 1000000 -center -pbc mol -ur compact',
						shell=True,
						stdout=f,
						stderr=f,
					)
					if not os.path.exists(f"{workdir}{prefix}-fullsystem.pdb"):
						print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
						print(f'\nGMX Error : GROMACS failed to produce "{workdir}{prefix}-fullsystem.pdb".\nSee logs/gmx-step... for details.\n\n')
						print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
						exit()

					subprocess.run(
						f'{gmx_q} editconf -f {workdir_q}{prefix_q}-fullsystem.pdb -o {workdir_q}{prefix_q}-fullsystem.pdb ' , #-center 0 0 0',
						shell=True,
						stdout=f,
						stderr=f,
					)

					output = subprocess.run(
						f'{gmx_q} editconf -f {workdir_q}{prefix_q}-fullsystem.pdb -n {workdir_q}template.ndx -o {workdir_q}{prefix_q}-template-prod.pdb -ndef',
						shell=True,
						stdout=f,
						stderr=subprocess.PIPE,
						text=True,
					)
					if self.gmx_error_check(output, "template_ndx", 501):
						return

					subprocess.run(
						f'{gmx_q} editconf -f {workdir_q}{prefix_q}-fullsystem.pdb -n {workdir_q}fms-all.ndx -o {workdir_q}{prefix_q}-fullcomplex-prod.pdb -ndef',
						shell=True,
						stdout=f,
						stderr=f,
					)

					if not os.path.exists(f"{workdir}{prefix}-fullcomplex-prod.pdb"):
						print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
						print(f'\nGMX Error : GROMACS failed to produce "{workdir}{prefix}-fullcomplex-prod.pdb".\nSee logs/gmx-step... for details.\n\n')
						print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
						exit()

				# Not this one though, this is for nips
				elif pargs.cplx and not pargs.template_name and not pargs.solvent:
					# NIP Condition
					subprocess.run(
						f'( echo "Other" ; echo "FMs" ) | {gmx_q} trjconv -f {workdir_q}rxn_prod.xtc -s {workdir_q}rxn_prod.tpr -n {workdir_q}fms-all.ndx -o {workdir_q}{prefix_q}-fullcomplex-prod.pdb -dump 1000000 -center -pbc mol -ur compact',
						shell=True,
						stdout=f,
						stderr=f,
					)
					
					subprocess.run(
						f'{gmx_q} editconf -f {workdir_q}{prefix_q}-fullcomplex-prod.pdb -o {workdir_q}{prefix_q}-fullcomplex-prod.pdb',
						shell=True,
						stdout=f,
						stderr=f,
					)

					if not os.path.exists(f"{workdir}{prefix}-fullcomplex-prod.pdb"):
						print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
						print(f'\nGMX Error : GROMACS failed to produce "{workdir}{prefix}-fullcomplex-prod.pdb".\nSee logs/gmx-step... for details.\n\n')
						print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
						exit()

				elif pargs.cplx and not pargs.template_name and pargs.solvent:

					subprocess.run(
						f'( echo "Other" ; echo "Other" ) | {gmx_q} trjconv -f {workdir_q}rxn_prod.xtc -s {workdir_q}rxn_prod.tpr -o {workdir_q}{prefix_q}-fullsystem.pdb -dump 1000000 -center -pbc mol -ur compact',
						shell=True,
						stdout=f,
						stderr=f,
					)
					if not os.path.exists(f"{workdir}{prefix}-fullsystem.pdb"):
						print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
						print(f'\nGMX Error : GROMACS failed to produce "{workdir}{prefix}-fullsystem.pdb".\nSee logs/gmx-step... for details.\n\n')
						print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
						exit()

					subprocess.run(
						f'{gmx_q} editconf -f {workdir_q}{prefix_q}-fullsystem.pdb -o {workdir_q}{prefix_q}-fullsystem.pdb ' , #-center 0 0 0',
						shell=True,
						stdout=f,
						stderr=f,
					)

					subprocess.run(
						f'{gmx_q} editconf -f {workdir_q}{prefix_q}-fullsystem.pdb -n {workdir_q}fms-all.ndx -o {workdir_q}{prefix_q}-fullcomplex-prod.pdb -ndef',
						shell=True,
						stdout=f,
						stderr=f,
					)

					if not os.path.exists(f"{workdir}{prefix}-fullcomplex-prod.pdb"):
						print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
						print(f'\nGMX Error : GROMACS failed to produce "{workdir}{prefix}-fullcomplex-prod.pdb".\nSee logs/gmx-step... for details.\n\n')
						print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
						exit()

				else:
					print("Failed solvent/fm/template separation? Report any errors from this on GitHub. If it works just ignore this.")
					subprocess.run(
						f'( echo "Other" ; echo "FMs" ) | {gmx_q} trjconv -f {workdir_q}rxn_prod.xtc -s {workdir_q}rxn_prod.tpr -n {workdir_q}fms-all.ndx -o {workdir_q}{prefix_q}-fullcomplex-prod.pdb -dump 1000000 -center -pbc mol -ur compact',
						shell=True,
						stdout=f,
						stderr=f,
					)

					subprocess.run(
						f'{gmx_q} editconf -f {workdir_q}{prefix_q}-fullcomplex-prod.pdb -o {workdir_q}{prefix_q}-fullcomplex-prod.pdb',
						shell=True,
						stdout=f,
						stderr=f,
					)


				if not pargs.solvent or pargs.solvent == 'spc':

					if pargs.protein:
						# center system
						subprocess.run(
							f'( echo "Protein" ; echo "Water" ) | {gmx_q} trjconv -f {workdir_q}rxn_prod.xtc -s {workdir_q}rxn_prod.tpr -o {workdir_q}spc.pdb -dump 1000000 -center -pbc mol -ur compact',
							shell=True,
							stdout=f,
							stderr=f,
						)
					elif pargs.template or pargs.cplx:
						# center system
						subprocess.run(
							f'( echo "Other" ; echo "Water" ) | {gmx_q} trjconv -f {workdir_q}rxn_prod.xtc -s {workdir_q}rxn_prod.tpr -o {workdir_q}spc.pdb -dump 1000000 -center -pbc mol -ur compact',
							shell=True,
							stdout=f,
							stderr=f,
						)

					if not os.path.exists(f"{workdir}spc.pdb"):
						print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
						print(f'\nGMX Error : GROMACS failed to produce "{workdir}spc.pdb".\nSee logs/gmx-step... for details.\n\n')
						print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
						exit()

			'''
			Maintaining initiator locations
			- separate out APS/TEMED/AZO etc. 
			- might need to do the same with solvents and ions at some point if it is 
				energetically important. Only issue is dealing with new bonds, which 
				in edge cases cause issues		
			'''
			if pargs.explicit:
				for initiator in pargs.initiator:
					subprocess.run(
						f'( echo "Other" ; echo {shlex.quote(initiator.upper())} ) | {gmx_q} trjconv -f {workdir_q}rxn_prod.xtc -s {workdir_q}rxn_prod.tpr -o {workdir_q}{prefix_q}-{shlex.quote(initiator)}-prod.pdb -dump 1000000 -center -pbc mol -ur compact',
						shell=True,
						stdout=f,
						stderr=f,
					)

			if run_type == "interact":

				gmx_make_mdp(gmxdir,workdir,'interact', pargs)

				#   Run for Energy Groups
				subprocess.run(
					f'{gmx_q} grompp -f {workdir_q}/fm_interact.mdp -c {workdir_q}rxn_npt.gro -r {workdir_q}rxn_npt.gro -t {workdir_q}rxn_npt.cpt -p {workdir_q}rxn_topol.top -o {workdir_q}interact_md.tpr -maxwarn 1',
					shell=True,
					stdout=f,
					stderr=f,
				)
				output = subprocess.run(
					f"{gmx_q} mdrun {thread_flags} {gpuid} {bonded_flag} -v -s {workdir_q}interact_md.tpr -deffnm {workdir_q}interact_md",
					shell=True,
					stdout=subprocess.PIPE,
					stderr=subprocess.PIPE,
					text=True,
				)
				if self.gmx_error_check(output, "interact_md", 569):
					return

				subprocess.run(
					f'cp {workdir_q}interact_md.gro {compdir_q}{prefix_q}_prod.gro ',
					shell=True,
					stdout=subprocess.DEVNULL,
					stderr=subprocess.DEVNULL,
				)
				subprocess.run(
					f'cp {workdir_q}interact_md.xtc {compdir_q}{prefix_q}_prod.xtc ',
					shell=True,
					stdout=subprocess.DEVNULL,
					stderr=subprocess.DEVNULL,
				)
				subprocess.run(
					f'cp {workdir_q}interact_md.tpr {compdir_q}{prefix_q}_prod.tpr ',
					shell=True,
					stdout=subprocess.DEVNULL,
					stderr=subprocess.DEVNULL,
				)
				subprocess.run(
					f'cp {workdir_q}interact_md.edr {compdir_q}{prefix_q}_prod.edr ',
					shell=True,
					stdout=subprocess.DEVNULL,
					stderr=subprocess.DEVNULL,
				)

				# Full Interaction Run
				#
				# This is for FMs
				#
				eng_groups_str = "energygrps = "
				if pargs.id:
					eng_groups = []
					for mol in pargs.recipe:
						if mol != "UNK":
							eng_groups_str += f"{enfm[mol]} "
							eng_groups.append(enfm[mol])
				else:
					for mol in eng_groups:
						if mol != "UNK":
							eng_groups_str += f"{mol} "

				if pargs.protein:
					eng_groups_str += f"Protein "
				elif pargs.template:
					eng_groups_str += f"UNK "

				# subprocess.run(
				# 	f"cp {gmxdir}/fm_eng.mdp {workdir}/fm_eng.mdp",
				# 	shell=True,
				# 	stdout=f,
				# 	stderr=f,
				# )

				gmx_make_mdp(gmxdir,workdir,'eng', pargs)

				subprocess.run(
					f"echo {shlex.quote(eng_groups_str)} >> {workdir_q}/fm_eng.mdp",
					shell=True,
					stdout=f,
					stderr=f,
				)

				subprocess.run(
					f"{gmx_q} grompp -f {workdir_q}fm_eng.mdp -c {workdir_q}interact_md.tpr -r {workdir_q}interact_md.tpr -p {workdir_q}rxn_topol.top -o {workdir_q}interact_eng_groups.tpr -maxwarn 1",
					shell=True,
					stdout=f,
					stderr=f,
				)

				output = subprocess.run(
					f"{gmx_q} mdrun {thread_flags} {gpuid} {bonded_flag} -v -s {workdir_q}interact_eng_groups.tpr -rerun {workdir_q}interact_md.xtc -e {workdir_q}interact_eng_groups.edr",
					shell=True,
					stdout=f,
					stderr=subprocess.PIPE,
					text=True,
				)
				if self.gmx_error_check(output, "interact_eng", 760):
					return

				if pargs.protein:
					
					# if hbond fails for protein, no reason to screen all
					hbond_protein_failure = True
					
					for fm in eng_groups:
						if fm != "UNK":
							if hbond_protein_failure:
								output = subprocess.run(
									f'( echo "Protein" ; echo {shlex.quote(fm.upper())} ) | {gmx_q} hbond-legacy -f {workdir_q}interact_md.xtc -s {workdir_q}interact_md.tpr -num {workdir_q}template-{shlex.quote(fm)}_hbond_num.xvg -dist {workdir_q}template-{shlex.quote(fm)}_hbond_dist.xvg',
									shell=True,
									stdout=f,
									stderr=subprocess.PIPE,
									text=True,
								)
								hbond_protein_failure = gmx_hbond_check(output, f"{fm}")

							output = subprocess.run(
								f'( echo "Protein" ; echo {shlex.quote(fm.upper())} ) | {gmx_q} rdf -f {workdir_q}interact_md.xtc -s {workdir_q}interact_md.gro -o {workdir_q}template-{shlex.quote(fm)}_rdf.xvg',
								shell=True,
								stdout=f,
								stderr=subprocess.PIPE,
								text=True,
							)


					for fm in eng_groups:
						if fm != "UNK":
							subprocess.run(
								f"echo LJ-SR:Protein-{shlex.quote(fm)} | {gmx_q} energy -f {workdir_q}interact_eng_groups.edr -s {workdir_q}interact_eng_groups.tpr -o {workdir_q}template-{shlex.quote(fm)}_lj.xvg",
								shell=True,
								stdout=f,
								stderr=f,
								)
							subprocess.run(
								f"echo Coul-SR:Protein-{shlex.quote(fm)} | {gmx_q} energy -f {workdir_q}interact_eng_groups.edr -s {workdir_q}interact_eng_groups.tpr -o {workdir_q}template-{shlex.quote(fm)}_coul.xvg",
								shell=True,
								stdout=f,
								stderr=f,
								)
							subprocess.run(
								f"echo LJ-SR:{shlex.quote(fm)}-Protein | {gmx_q} energy -f {workdir_q}interact_eng_groups.edr -s {workdir_q}interact_eng_groups.tpr -o {workdir_q}template-{shlex.quote(fm)}_lj.xvg",
								shell=True,
								stdout=f,
								stderr=f,
							)
							subprocess.run(
								f"echo Coul-SR:{shlex.quote(fm)}-Protein | {gmx_q} energy -f {workdir_q}interact_eng_groups.edr -s {workdir_q}interact_eng_groups.tpr -o {workdir_q}template-{shlex.quote(fm)}_coul.xvg",
								shell=True,
								stdout=f,
								stderr=f,
							)
				elif pargs.template:
					for fm in eng_groups:
						for mol in template_groups:
							if mol != fm:
								subprocess.run(
									f"echo LJ-SR:{shlex.quote(mol)}-{shlex.quote(fm)} | {gmx_q} energy -f {workdir_q}interact_eng_groups.edr -s {workdir_q}interact_eng_groups.tpr -o {workdir_q}{shlex.quote(mol)}-{shlex.quote(fm)}_lj.xvg",
									shell=True,
									stdout=f,
									stderr=f,
								)
								subprocess.run(
									f"echo Coul-SR:{shlex.quote(mol)}-{shlex.quote(fm)} | {gmx_q} energy -f {workdir_q}interact_eng_groups.edr -s {workdir_q}interact_eng_groups.tpr -o {workdir_q}{shlex.quote(mol)}-{shlex.quote(fm)}_coul.xvg",
									shell=True,
									stdout=f,
									stderr=f,
								)
								#
								# This is dumb but i havent figured out how to maintain an order with gromacs yet
								# might have something to do with order of energy groups from set()
								#
								subprocess.run(
									f"echo LJ-SR:{shlex.quote(fm)}-{shlex.quote(mol)} | {gmx_q} energy -f {workdir_q}interact_eng_groups.edr -s {workdir_q}interact_eng_groups.tpr -o {workdir_q}{shlex.quote(fm)}-{shlex.quote(mol)}_lj.xvg",
									shell=True,
									stdout=f,
									stderr=f,
								)
								subprocess.run(
									f"echo Coul-SR:{shlex.quote(fm)}-{shlex.quote(mol)} | {gmx_q} energy -f {workdir_q}interact_eng_groups.edr -s {workdir_q}interact_eng_groups.tpr -o {workdir_q}{shlex.quote(fm)}-{shlex.quote(mol)}_coul.xvg",
									shell=True,
									stdout=f,
									stderr=f,
								)

		incomplex = f"{workdir}{prefix}-fullcomplex-prod.pdb"
		if pargs.protein:
			intemplate = f"{workdir}{prefix}-protein-prod.pdb"
		elif pargs.template:
			intemplate = f"{workdir}{prefix}-template-prod.pdb"
		else:
			intemplate = None

		if pargs.explicit:
			ininitiator = f"{workdir}{prefix}-initiator-prod.pdb"

			initiators = []
			for initiator in pargs.initiator:
				initiators.extend(read_pdb(f"{workdir}{prefix}-{initiator}-prod.pdb", hydrogens=True, initiators=True))
			
			print_molecule(initiators, ininitiator)

		else:
			ininitiator = None

		if pargs.clean and run_type=="rxn":
			# if it is the last step make a copy before cleaning
			if (i + 1) == maxi:
				# first template and complex
				subprocess.run(f"cp {shlex.quote(incomplex)} {workdir_q}final-complex-prod.pdb", shell=True)
				if pargs.protein:
					subprocess.run(f"cp {shlex.quote(intemplate)} {workdir_q}final-template-prod.pdb", shell=True)
					subprocess.run(f"cp {workdir_q}template-clean.itp {workdir_q}final-template.itp", shell=True)
					subprocess.run(f"cp {workdir_q}template_posres.itp {workdir_q}final-template-posres.itp", shell=True)

				# forcefield parameters
				# subprocess.run(f"cp {workdir}rxn_topol.top {workdir}final-topol.top", shell=True)
				subprocess.run(f"cp {workdir_q}RXN_merged.itp {workdir_q}final-complex.itp", shell=True)
				if pargs.protein:
					gtop = generic_top(
						pargs,
						'#include "final-template.itp"',
						"Template".rjust(7) + "1".rjust(18),
						'#include "final-complex.itp"',
						"CLX".rjust(7) + "1".rjust(18),
						"final-template-posres.itp",
					)
				else:
					gtop = generic_top(
						pargs,
						"",
						"",
						'#include "final-complex.itp"',
						"CLX".rjust(7) + "1".rjust(18),
					)

				_top = open(f"{workdir}final-topol.top", "w+")
				_top.write(gtop)
				_top.close()

			subprocess.run(f"rm {workdir_q}rxn*", shell=True)
			subprocess.run(f"rm {workdir_q}RXN*", shell=True)
			subprocess.run(f"rm {workdir_q}\\#*", shell=True)

		return incomplex, intemplate, ininitiator

	def gmx_error_check(self, output, loc, lineno):
		# this only typeerrors for the interaction one?
		# try:
		output = output.stderr.split('\n')
		full_text = '\n'.join(output)
		error = None
		error_code = ""
		if "Error in user input:" in full_text:
			error_code+="\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
			error_code+=f"\nError in user input at GMX : {loc}\n\nLine {lineno}\n\nStacktrace:\n"
			for line in output:
				error_code+=line+'\n'
			error_code+=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
			error = True		
		elif "Range checking error (possible bug):" in full_text:
			error_code+="\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
			error_code+=f"\nError in user input at GMX : {loc}\n\nThis is potentially an issue with your CUDA driver. Check `nvidia-smi`\n\nLine {lineno}\n\nStacktrace:\n"
			for line in output:
				error_code+=line+'\n'
			error_code+=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
			error = True
		elif "Segmentation fault (core dumped)" in full_text:
			error_code+="\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
			error_code+=f"\nError GMX (likely energy related): {loc}\n\nLine {lineno}\n\nStacktrace:\n"
			for line in output:
				error_code+=line+'\n'
			error_code+=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
			error = True
		elif "Aborted (core dumped)" in full_text:
			error_code+="\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
			error_code+=f"\nError GMX (likely energy related): {loc}\n\nLine {lineno}\n\nStacktrace:\n"
			for line in output:
				error_code+=line+'\n'
			error_code+=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
			error = True
		elif "Fatal error:" in full_text:
			error_code+="\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
			error_code+=f"\nError GMX indexing: {loc}\nDouble check if -template was used instead of -protein for a protein PDB.\n\nLine {lineno}\n\nStacktrace:\n"
			for line in output:
				error_code+=line+'\n'
			error_code+=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
			error = True
		elif "Warning: Only triclinic boxes with the first vector parallel to the x-axis and the second vector in the xy-plane are supported." in full_text:
			error_code+="\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
			error_code+=f"\nError GMX in pressure scaling: {loc}\n\nLine {lineno}\n\nStacktrace:\n"
			for line in output:
				error_code+=line+'\n'
			error_code+=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
			error = True
		elif "Assertion failed:" in full_text:
			error_code+="\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
			error_code+=f"\nError GMX (likely energy related): {loc}\n\nLine {lineno}\n\nStacktrace:\n"
			for line in output:
				error_code+=line+'\n'
			error_code+=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
			error = True

		if error:
			if self.pargs.retry and self.pargs.dt['rxn'] > 1E-4:
				self.pargs.dt['rxn'] = self.pargs.dt['rxn'] / 2
				print("   - Failed at %s Retrying with dt : %s" % (loc, self.pargs.dt['rxn']))
				self.run()
				# signal the caller (run()) to stop -- self.run() above already
				# re-executed the whole pipeline from scratch; letting the
				# caller keep going would run a second, overlapping pipeline
				# against the same output files
				return True
			# elif not self.rolled:
			# 	print("Failed with dt : %s" % (self.pargs.dt['rxn']))
			# 	print("Rolling back one step")
			# 	self.rolled = True
			# 	self.i += -1
			# 	self.pargs.dt['rxn'] = self.init_dt
			# 	self.run()
			print("\t Failed with dt : %s" % (self.pargs.dt['rxn']))
			print(error_code)
			sys.exit()
		return False

def gmx_hbond_check(output, fm):
	output = output.stderr.split('\n')
	if "Selection 'Protein' has no donors AND has no acceptors! Nothing to be done." in output:
		print(f"No hydrogen bonds found for {fm}.")
		print("Selection 'Protein' has no donors AND has no acceptors! Nothing to be done.")
		return False
	else:
		return True

def gmx_make_mdp(gmxdir,workdir,run_type,pargs):

	gmxdir_q = shlex.quote(gmxdir)
	workdir_q = shlex.quote(workdir)
	run_type_q = shlex.quote(run_type)

	# print(run_type)
	try:
		if run_type == 'eng':
			dt = pargs.dt['interact']
			nsteps = pargs.nsteps['interact']
		else:
			dt = pargs.dt[run_type]
			nsteps = pargs.nsteps[run_type]
	except KeyError as e:
		print(f'Make .mdp has an error: missing configuration for run_type "{run_type}". Exiting.')
		sys.exit(1)


	if run_type in ['nvt_vacuum', 'nvt_stage0', 'nvt_stage1', 'nvt_stage2', 'nvt_stage3']:
		subprocess.run(
			f'cp {gmxdir_q}/fm_nvt.mdp {workdir_q}/fm_{run_type_q}.mdp',
			shell=True,
		)
	else:
		subprocess.run(
			f'cp {gmxdir_q}/fm_{run_type_q}.mdp {workdir_q}/fm_{run_type_q}.mdp',
			shell=True,
		)
	subprocess.run(
		f"echo '\n;MIPkit Run Parameters' >> {workdir_q}/fm_{run_type_q}.mdp",
		shell=True,
	)
	subprocess.run(
		f"echo 'nsteps\t\t\t\t\t= {nsteps}' >> {workdir_q}/fm_{run_type_q}.mdp",
		shell=True,
	)
	if run_type != "min":

		subprocess.run(
			f"echo 'dt\t\t\t\t\t\t= {dt}' >> {workdir_q}/fm_{run_type_q}.mdp",
			shell=True,
		)
		subprocess.run(
			f"echo 'ref_t\t\t\t\t\t= {pargs.temp}' >> {workdir_q}/fm_{run_type_q}.mdp",
			shell=True,
		)