#!/bin/bash
#SBATCH -J gamboot
#SBATCH --account somerville_lab
#SBATCH --mem 20G
#SBATCH -p ncf
#SBATCH --cpus-per-task 10
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 1-00:00
#SBATCH -o log/%x_%A_%a.out
#SBATCH --mail-user=john_flournoy@fas.harvard.edu
#SBATCH --mail-type=ALL

module load gcc/7.1.0-fasrc01
module load R/3.5.1-fasrc01

cp ~/.R/Makevars.gcc ~/.R/Makevars

export R_LIBS_USER=/ncf/mclaughlin/users/jflournoy/R_3.5.1_GCC:$R_LIBS_USER

runme=/net/holynfs01/srv/export/mclaughlin/share_root/users/jflournoy/code/hcpd_gam_bootstrap/boot_gam_deriv.R

srun -c 10 `which Rscript` "${runme}"
