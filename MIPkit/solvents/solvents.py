
import shlex
import subprocess

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, GetPeriodicTable, Descriptors, Draw
from rdkit.Geometry import Point3D
from rdkit.Chem.rdchem import BondType

gmx_q = shlex.quote("gmx")

def run_step(cmd, label):
	result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
	if result.returncode != 0:
		print(f"solvents.py: step '{label}' failed (exit code {result.returncode}):\n{result.stderr}")
		raise RuntimeError(f"solvents.py step '{label}' failed")
	return result

solvents = {
			"DMF" :"CN(C)C=O",
			"DMS" :"CS(=O)C",
			# "MCN" :"CC#N",
			# "CHL":"C(Cl)(Cl)Cl"
			}

for sol in solvents:

	mol = Chem.MolFromSmiles(solvents[sol])

	params = AllChem.ETKDGv3()
	params.useRandomCoords = True
	# params.maxAttempts = 10000
	params.enforceChirality = False


	mol = Chem.rdmolops.AddHs(mol, addCoords=True)
	Chem.SanitizeMol(mol)
	AllChem.EmbedMolecule(mol, params)

	# https://sourceforge.net/p/rdkit/mailman/message/36404394/
	mi = Chem.AtomPDBResidueInfo()
	mi.SetResidueName(sol)

	[a.SetMonomerInfo(mi) for a in mol.GetAtoms()]

	w = Chem.PDBWriter(
		f"{sol}.pdb"
	)
	w.write(mol)
	w.close()

	sol_q = shlex.quote(sol)

	subprocess.run(
		f"rm #*",
		shell=True
	)

	if sol == "DMF":
		run_step(f'{gmx_q} editconf -f {sol_q}.pdb -o {sol_q}.gro -box 0.6 0.6 0.6', "editconf")
	elif sol == "CHL":
		run_step(f'{gmx_q} editconf -f {sol_q}.pdb -o {sol_q}.gro -box 0.7 0.7 0.7', "editconf")
	else:
		run_step(f'{gmx_q} editconf -f {sol_q}.pdb -o {sol_q}.gro -box 0.5 0.5 0.5', "editconf")

	run_step(f'acpype -i {sol_q}.pdb -b {sol_q} -q sqm -c bcc ', "acpype")

	subprocess.run(
		f"mkdir {sol_q}-gmx",
		shell=True
	)

	subprocess.run(
		f" cp {sol_q}.acpype/{sol_q}_GMX.itp {sol_q}.itp  ",
		shell=True,
	)

	subprocess.run(
		f" cp {sol_q}.acpype/{sol_q}_GMX.itp {sol_q}-gmx/{sol_q}_GMX.itp  ",
		shell=True,
	)

	subprocess.run(
		f" cp {sol_q}.acpype/{sol_q}_GMX.top {sol_q}-gmx/topol.top  ",
		shell=True,
	)

	subprocess.run(
		f" cp {sol_q}.gro {sol_q}-gmx/{sol_q}.gro  ",
		shell=True,
	)

	workdir=f"{sol}-gmx"
	workdir_q = shlex.quote(workdir)

	subprocess.run(
		f"rm {workdir_q}/#*",
		shell=True
	)

	run_step(f'{gmx_q} solvate -cp {sol_q}.gro -cs {sol_q}.gro -p {workdir_q}/topol.top -o {workdir_q}/{sol_q}.gro -box 7 7 7', "solvate")

	run_step(f'{gmx_q} grompp -f "gmx/fm_min.mdp" -c {workdir_q}/{sol_q}.gro -p {workdir_q}/topol.top -o {workdir_q}/sol_min.tpr -maxwarn 1 ', "grompp (min)")

	run_step(f'{gmx_q} mdrun -nt 8 -gpu_id 0 -v -deffnm {workdir_q}/sol_min', "mdrun (min)")

	run_step(f'{gmx_q} grompp -f "gmx/fm_nvt.mdp" -c {workdir_q}/{sol_q}.gro -p {workdir_q}/topol.top -o {workdir_q}/sol_nvt.tpr  -maxwarn 1', "grompp (nvt)")

	run_step(f'{gmx_q} mdrun -nt 8 -gpu_id 0 -v -deffnm {workdir_q}/sol_nvt', "mdrun (nvt)")

	run_step(f'{gmx_q} grompp -f "gmx/fm_npt.mdp" -c {workdir_q}/{sol_q}.gro -p {workdir_q}/topol.top -o {workdir_q}/sol_npt.tpr  -maxwarn 1', "grompp (npt)")

	run_step(f'{gmx_q} mdrun -nt 8 -gpu_id 0 -v -deffnm {workdir_q}/sol_npt', "mdrun (npt)")

	subprocess.run(
		f'cp "gmx/fm_md_long.mdp" {workdir_q}/fm_md_long.mdp',
		shell=True,
	)

	run_step(f'{gmx_q} grompp -f {workdir_q}/fm_md_long.mdp -c {workdir_q}/{sol_q}.gro -p {workdir_q}/topol.top -o {workdir_q}/sol_prod.tpr  -maxwarn 1', "grompp (prod)")

	run_step(f'{gmx_q} mdrun -nt 8 -gpu_id 0 -v -deffnm {workdir_q}/sol_prod', "mdrun (prod)")

	subprocess.run(
		f'cp {workdir_q}/sol_prod.gro {sol_q}.gro',
		shell=True,
	)

	subprocess.run(
		f'(echo "Density") |{gmx_q} energy -f {sol_q}-gmx/sol_prod.edr -o {sol_q}-density.xvg',
		shell=True,
	)
	subprocess.run(
		f'(echo "Total-Energy") | {gmx_q} energy -f {sol_q}-gmx/sol_prod.edr -o {sol_q}-eng.xvg',
		shell=True,
	)

