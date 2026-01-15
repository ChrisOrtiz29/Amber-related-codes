#!/bin/bash

input_prot_pdb=$1
output_in=$2
output_top=$3
output_crd=$4
output_pdb=$5
num_pos_ion=$6
num_neg_ion=$7

cat > ${output_in} << EOF
set default PBRadii mbondi3
source leaprc.gaff2
source leaprc.protein.ff14SB
source leaprc.RNA.OL3
source leaprc.DNA.0L24
source leaprc.water.tip3p
model = loadpdb ${input_prot_pdb}
solvatebox model TIP3PBOX 10 iso
addIons model Na+ 0
addIons model Na+ ${num_pos_ion}
addIons model Cl- ${num_neg_ion}
saveamberparm model ${output_top} ${output_crd}
savepdb model ${output_pdb}
quit
EOF

tleap -f ${output_in}
