#!/bin/bash

input_prot_pdb=$1
input_lig_pdb=$2
output_in=$3
output_top=$4
output_crd=$5
output_pdb=$6
num_pos_ion=$7
num_neg_ion=$8

cat > ${output_in} << EOF
source leaprc.gaff
source leaprc.protein.ff14SB
source leaprc.water.tip3p
loadamberparams /mnt/Tsunami_HHD/wendy/mmgbsa/3SVG/ODR.frcmod
loadoff /mnt/Tsunami_HHD/wendy/mmgbsa/3SVG/ODR.lib
protein = loadpdb ${input_prot_pdb}
lig = loadpdb ${input_lig_pdb}
model = combine {protein lig}
solvatebox model TIP3PBOX 10 iso
addIons model Na+ 0
addIons model Na+ ${num_pos_ion}
addIons model Cl- ${num_neg_ion}
saveamberparm model ${output_top} ${output_crd}
savepdb model ${output_pdb}
quit
EOF

tleap -f ${output_in}
