import os, subprocess
import shlex
from MIPkit.utils.print_molecule import print_molecule
from MIPkit.utils.extract_heavies import print_heavies, print_backbone
from MIPkit.utils.read_itp import (
	merge_itps,
	merge_atomtypes,
	replace_residues,
)
from MIPkit.utils.read_pdb import read_pdb

def acpype_do(run_type, i, step_complex, out_complex, pargs, logname):
	if run_type == "min":
		prefix = "min"
	elif run_type == "rxn":
		prefix = f"step{i}"
	elif run_type == "interact":
		prefix = "interact"
	elif run_type == "screen":
		prefix = "interact"
	else:
		print("Incorrect run type, this should not happen.")
		exit()

	if pargs.force:
		force = "-f"
	else:
		force = ""

	acpype = pargs.acpype
	gmx = pargs.gmx
	obabel = pargs.obabel

	acpype_q = shlex.quote(acpype)
	gmx_q = shlex.quote(gmx)
	obabel_q = shlex.quote(obabel)

	ll = pargs.parallel
	ll_q = shlex.quote(ll)

	charge_q = shlex.quote(pargs.charge)
	acff_q = shlex.quote(pargs.acff)

	if os.path.exists(f'RXN{ll}.acpype'):
		subprocess.run(
			f"rm RXN{ll_q}.acpype/*", 
			shell=True, 
			stdout = subprocess.DEVNULL, 
			stderr = subprocess.DEVNULL)
	
	if os.path.exists(f'RES{ll}.acpype'):
		subprocess.run(
			f"rm RES{ll_q}.acpype/*", 
			shell=True, 
			stdout = subprocess.DEVNULL, 
			stderr = subprocess.DEVNULL)

	compdir = pargs.compdir
	logdir = pargs.logdir
	initiatordir = pargs.initiatordir
	solventdir = pargs.solventdir
	workdir = pargs.workdir

	compdir_q = shlex.quote(compdir)
	workdir_q = shlex.quote(workdir)

	numparts = len(out_complex)

	# ===========================================================
	#
	# Main Do
	#
	# ===========================================================
	acpype_types = {}
	with open(os.path.join(logdir, logname), "w+") as f:
		for j in range(numparts):

			fm = step_complex[j]["name"]

			if fm in acpype_types:
				fmitp = acpype_types[fm] 
				subprocess.run(
					f"cp {workdir_q}{shlex.quote(fmitp)} {workdir_q}RXN_p{j}.itp", 
					shell=True, 
					stdout = subprocess.DEVNULL, 
					stderr = subprocess.DEVNULL)

			else:
				subprocess.run(
					f'{acpype_q} -i {compdir_q}{shlex.quote(prefix)}-part{j}.pdb -b RXN{ll_q} -q sqm -c {charge_q} -n 0 {force} -a {acff_q}',
					shell=True,
					stdout=f,
					stderr=f,
				)
				
				if not os.path.exists(f'RXN{ll}.acpype/RXN{ll}_GMX.itp'):
					print("----------------------------------------------------------------------")
					print(f'\nAcpype Error : Acpype failed on complex part {j}.\nSee logs/acpype-step... for details.\n\n')
					print("----------------------------------------------------------------------")
					exit()

				subprocess.run(
					f" cp RXN{ll_q}.acpype/RXN{ll_q}_GMX.itp {workdir_q}RXN_p{j}.itp  ",
					shell=True,
					stdout=f,
					stderr=f,
				)
				subprocess.run(
					f"rm RXN{ll_q}.acpype/RXN*", 
					shell=True, 
					stdout = subprocess.DEVNULL, 
					stderr = subprocess.DEVNULL)

				acpype_types[fm] = f"RXN_p{j}.itp"

		# Merge .itps
		itps = []
		fms = []
		for j in range(numparts):
			itps.append(f"{workdir}RXN_p{j}.itp")
			fms.append(step_complex[j]["name"])

		merge_itps(itps, f"{workdir}RXN_merged.itp", fms, molname="CLX")

		if pargs.id and pargs.interact:
			fms = []
			for j in range(numparts):
				part = read_pdb(f"{compdir}{prefix}-part{j}-id.pdb")
				for fm in part:
					for r in fm["atoms"]:
						fms.append(fm["name"])

			replace_residues(
				f"{workdir}RXN_merged.itp",
				f"{workdir}RXN_merged.itp",
				fms,
				molname="CLX",
			)


		if pargs.template_name:
			template_name = pargs.template_name
			posre = pargs.posre
			template_name_q = shlex.quote(template_name)

			# just do once
			if (i == 0 and not pargs.min) or run_type == "min" or run_type == "interact" or run_type == "screen":
				if pargs.protein and run_type != "interact":
					if not pargs.gmxprotein:
						output = subprocess.run(
							f"(echo 2 ; echo 7) | {gmx_q} pdb2gmx -f {template_name_q} -o {workdir_q}template.pdb -p {workdir_q}template.top -i {workdir_q}template_posres.itp -posrefc {posre}",
							shell=True,
							stdout=f,
							stderr=subprocess.PIPE,
							text=True
						)
						acpype_error_check(output,"pdb2gmx",112)

					else:
						template = read_pdb(f"{template_name}", hydrogens=True)

						subprocess.run(
							f"cp {template_name_q} {workdir_q}template.pdb",
							shell=True,
							stdout=f,
							stderr=f,
						)

						print_backbone(template,
							f"{workdir}backbone.pdb",
							conect=False,
							encode=False,
						)

						subprocess.run(
							f"(echo 0) | {gmx_q} genrestr -f {workdir_q}backbone.pdb -o {workdir_q}template_posres.itp -fc {posre}",
							shell=True,
							stdout=f,
							stderr=f,
						)
				
				elif pargs.protein and run_type == "interact":
					# generate constraints

					template = read_pdb(f"{template_name}", hydrogens=True)

					print_backbone(template,
						f"{workdir}backbone.pdb",
						conect=False,
						encode=False,
					)

					subprocess.run(
						f"(echo 0) | {gmx_q} genrestr -f {workdir_q}backbone.pdb -o {workdir_q}template_posres.itp -fc {posre}",
						shell=True,
						stdout=f,
						stderr=f,
					)

				else: # template only
					subprocess.run(
						f"cp {template_name_q} {workdir_q}template.pdb",
						shell=True,
						stdout=f,
						stderr=f,
					)
					
					template = read_pdb(f"{template_name}", hydrogens=True)

					print_heavies(template,
						f"{workdir}heavies.pdb",
						conect=False,
						encode=False,
					)

					subprocess.run(
						f"(echo 0) | {gmx_q} genrestr -f {workdir_q}heavies.pdb -o {workdir_q}template_posres.itp -fc {posre}",
						shell=True,
						stdout=f,
						stderr=f,
					)

				#
				# For Interactions with NIPs, offset
				#
				if run_type == "interact" and pargs.offset:
					ox,oy,oz = pargs.offset
					subprocess.run(
						f'{gmx_q} editconf -f {workdir_q}template.pdb -o {workdir_q}template.pdb -center {ox} {oy} {oz}',
						shell=True,
						stdout=f,
						stderr=f,
					)

				template = read_pdb(f"{workdir}template.pdb", hydrogens=True)

				single_res = []
				residues = []
				for res in template:
					for r in res["atoms"]:
						residues.append(res["name"])

					res["name"] = "UNK"
					single_res.append(res)

				print_molecule(
					single_res,
					f"{compdir}single_res_template.pdb",
					conect=False,
					encode=False,
				)
				subprocess.run(
					f'{acpype_q} -i {compdir_q}single_res_template.pdb -b RES{ll_q} -q sqm -c {charge_q} -n 0 {force} -a {acff_q}',
					shell=True,
					stdout=f,
					stderr=f,
				)
				subprocess.run(
					f" cp RES{ll_q}.acpype/RES{ll_q}_GMX.itp {workdir_q}single_res_template.itp  ",
					shell=True,
					stdout=f,
					stderr=f,
				)

				subprocess.run(
					f" cp {compdir_q}single_res_template.pdb {workdir_q}template.pdb ",
					shell=True,
					stdout=f,
					stderr=f,
				)

				replace_residues(
					f"{workdir}single_res_template.itp",
					f"{workdir}template.itp",
					residues,
					molname="Template",
				)

			else:
				if pargs.template:
					template = read_pdb(pargs.intemplate)
				elif pargs.protein:
					template = read_pdb(pargs.intemplate)

				intemplate = pargs.intemplate
				intemplate_q = shlex.quote(intemplate)

				residues = []
				for res in template:
					for r in res["atoms"]:
						residues.append(res["name"])
					res["name"] = "UNK"

				subprocess.run(
					f'{obabel_q} -ipdb {intemplate_q} -opdb -O {compdir_q}single_res_template.pdb ',
					shell=True,
					stdout=f,
					stderr=f,
				)

				subprocess.run(
					f" cp {compdir_q}single_res_template.pdb {workdir_q}template.pdb ",
					shell=True,
					stdout=f,
					stderr=f,
				)

				replace_residues(
					f"{workdir}single_res_template.itp",
					f"{workdir}template.itp",
					residues,
					molname="Template",
				)

			merge_atomtypes(
				itpfrom=f"{workdir}template.itp",
				itpinto=f"{workdir}RXN_merged.itp",
				cleaned=f"{workdir}template-clean.itp",
			)

		if pargs.explicit:
			merge_atomtypes(
				itpfrom=f"{initiatordir}/initiator_atoms.itp",
				itpinto=f"{workdir}RXN_merged.itp",
				cleaned=f"{workdir}RXN_tmp.itp",
			)

		if pargs.solvent and pargs.solvent != 'spc':
			merge_atomtypes(
				itpfrom=f"{solventdir}/solvent_atoms.itp",
				itpinto=f"{workdir}RXN_merged.itp",
				cleaned=f"{workdir}RXN_tmp.itp",
			)

def acpype_error_check(output, loc, line):
	# this only typeerrors for the interaction one?
	# try:
	output = output.stderr.split('\n')
	full_text = '\n'.join(output)
	if "Fatal error:" in full_text:
		print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
		print(f"\nError in user input at : {loc}\nAre you using a PDB that has already been processed by Gromacs? Gromacs formats exported PDBs differently than typical ones. Fix the file using -make_cartoon.\n\nLine {line}\n\nStacktrace:\n")
		for line in output:
			print(line)
		print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
		exit()