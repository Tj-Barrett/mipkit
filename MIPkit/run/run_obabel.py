import subprocess, os, shlex
from MIPkit.utils.print_molecule import print_molecule
from MIPkit.utils.read_pdb import (
	read_pdb,
	clean_pdb)

from MIPkit.utils.identify import identify


def obabel_initiators(i, aps_temed, pargs, logname):

	obabel = pargs.obabel
	obabel_q = shlex.quote(obabel)

	logdir = pargs.logdir
	workdir = pargs.workdir
	workdir_q = shlex.quote(workdir)

	if os.path.isfile(os.path.join(logdir, f"{logname}")):
		fileoperation = "a+"
	else:
		fileoperation = "w+"

	with open(os.path.join(logdir, f"{logname}"), fileoperation) as f:
		print_molecule(
			aps_temed,
			f"{workdir}step{i}-initiator-prod.pdb",
			conect=False,
			encode=False,
		)

		subprocess.run(
			f'{obabel_q} {workdir_q}step{i}-fullcomplex.pdb {workdir_q}step{i}-initiator-prod.pdb -ipdb -opdb -O {workdir_q}step{i}-fullcomplex.pdb',
			shell=True,
			stdout=f,
			stderr=f,
		)
		clean_pdb(f"{workdir}step{i}-fullcomplex.pdb")

def obabel_combine(run_type, i, out_complex, pargs, logname):

	obabel = pargs.obabel
	obabel_q = shlex.quote(obabel)

	logdir = pargs.logdir
	workdir = pargs.workdir
	workdir_q = shlex.quote(workdir)

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

	if pargs.template_name:
		proteinpdb = pargs.template_name
	else:
		proteinpdb = False

	with open(os.path.join(logdir, f"{logname}"), "w+") as f:
		print_molecule(
			out_complex,
			f"{workdir}{prefix}-fullcomplex.pdb",
			
			conect=False,
			encode=False,
		)

		if proteinpdb:
			# merge protein pdb with complex
			subprocess.run(
				f'{obabel_q} {workdir_q}template.pdb {workdir_q}{prefix_q}-fullcomplex.pdb -ipdb -opdb -O {workdir_q}{prefix_q}-fullcomplex.pdb',
				shell=True,
				stdout=f,
				stderr=f,
			)
			clean_pdb(f"{workdir}{prefix}-fullcomplex.pdb")

def obabel_solvent(run_type, i, pargs, outname, outtype, logname):

	obabel = pargs.obabel
	obabel_q = shlex.quote(obabel)

	logdir = pargs.logdir
	workdir = pargs.workdir
	workdir_q = shlex.quote(workdir)

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
	outname_q = shlex.quote(outname)
	outtype_q = shlex.quote(outtype)

	with open(os.path.join(logdir, f"{logname}"), "a+") as f:
		if pargs.solvent:
			solvent_q = shlex.quote(pargs.solvent)
			# merge protein pdb with complex
			subprocess.run(
				f'{obabel_q} {workdir_q}{prefix_q}-fullcomplex.pdb {workdir_q}{solvent_q}.pdb -ipdb -o{outtype_q} -O {workdir_q}{outname_q}.{outtype_q} ',
				shell=True,
				stdout=f,
				stderr=f,
			)
			if outtype == 'pdb':
				clean_pdb(f"{workdir}{outname}.{outtype}")
		else:
			# merge protein pdb with complex
			subprocess.run(
				f'{obabel_q} {workdir_q}{prefix_q}-fullcomplex.pdb {workdir_q}spc.pdb -ipdb -o{outtype_q} -O {workdir_q}{outname_q}.{outtype_q} ',
				shell=True,
				stdout=f,
				stderr=f,
			)
			if outtype == 'pdb':
				clean_pdb(f"{workdir}{outname}.{outtype}")

def obabel_do(run_type, i, min_complex, pargs, logname):

	obabel = pargs.obabel
	obabel_q = shlex.quote(obabel)

	compdir = pargs.compdir
	workdir = pargs.workdir
	workdir_q = shlex.quote(workdir)
	logdir = pargs.logdir

	if run_type == "min":
		prefix = "min"
	elif run_type == "interact":
		prefix = "interact"
	elif run_type == "screen":
		prefix = "interact"
	elif run_type == "rxn":
		prefix = f"step{i}"
	else:
		print("Incorrect run type, this should not happen.")
		exit()

	prefix_q = shlex.quote(prefix)

	with open(os.path.join(logdir, f"{logname}"), "w+") as f:
		out_complex = []
		for j, part in enumerate(min_complex):
			if part["name"][0:1] == "X":
				print_molecule(
					[part],
					f"{compdir}{prefix}-part{j}.pdb",
					
					conect=False,
					encode=False,
				)
				if pargs.id and pargs.interact:
					identify(pargs, pargs.recipe, f"{compdir}{prefix}-part{j}.pdb", None)
			else:
				print_molecule(
					[part],
					f"{compdir}{prefix}-part{j}.pdb",
					
					conect=False,
					encode=True,
				)
			tmp = read_pdb(f"{compdir}{prefix}-part{j}.pdb", hydrogens=True)[0]
			# read pdb doest read bonds, so get them from the step_complex
			tmp["bonds"] = part["bonds"]

			out_complex.append(tmp)

		print_molecule(
			out_complex,
			f"{workdir}{prefix}-fullcomplex.pdb",
			conect=False,
			encode=False,
		)
		subprocess.run(
			f'{obabel_q} -ipdb {workdir_q}{prefix_q}-fullcomplex.pdb -opdb -O {workdir_q}{prefix_q}-fullcomplex.pdb --minimize',
			shell=True,
			stdout=f,
			stderr=f,
		)
		return out_complex