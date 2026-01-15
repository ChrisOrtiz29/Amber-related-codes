#!/bin/bash
#SBATCH --job-name=mmgbsa
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=user@gmail.com
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --gres=gpu:1g.5gb:1
#SBATCH --output=mmgbsa.log
#SBATCH --partition=COMPUTE1Q
#SBATCH --account=yanglab

module load openmpi/2.1.1

mpirun -np 10 singularity exec --nv --bind /raid:/raid /raid/images/amber20.sif MMPBSA.py.MPI \
-O -i mmgbsa.in -o result_mmgbsa.dat -sp solvated_topol.prmtop \
-cp dry_cmp.parm7 -rp dry_rcp.parm7 -lp dry_lgd.parm7 -y solvated_md_NPT.nc  -do mmgbsa.do -eo mmgbsa.eo -deo mmgbsa.deo
