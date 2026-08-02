workdir="CD20-MD-3-MIP/"
scriptdir="/media/tj/Gromacs/Polymerization/250528-Rebinding/Scripts/"

rm ${workdir}\#*
# MIP
(echo "Other") | gmx trjconv -f ${workdir}epitope_md.gro -s ${workdir}epitope_md.tpr -o ${workdir}mip.pdb -conect yes #-center yes
(echo "Protein") | gmx trjconv -f ${workdir}epitope_md.gro -s ${workdir}epitope_md.tpr -o ${workdir}template.pdb -conect yes #-center yes

#cat ${workdir}mip.pdb ${workdir}template.pdb > ${workdir}mip-template.pdb

#gmx editconf -f ${workdir}mip-template.gro -center 0 0 0 -aligncenter 0 0 1 -o ${workdir}mip-template.gro

vmd -m ${workdir}template.pdb ${workdir}mip.pdb -e ${scriptdir}/VMD_MIP.tcl
