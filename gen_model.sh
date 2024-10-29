#!/bin/bash

input_prot_pdb=$1
output_in=$2
output_top=$3
output_crd=$4
output_pdb=$5

cat > ${output_in} << EOF
source leaprc.gaff
source leaprc.protein.ff14SB
source leaprc.water.tip3p
model = loadpdb ${input_prot_pdb}
solvatebox model TIP3PBOX 10 iso
addIons model Na+ 0
addIons model Cl- 0
saveamberparm model ${output_top} ${output_crd}
savepdb model ${output_pdb}
quit
EOF

tleap -f ${output_in}
