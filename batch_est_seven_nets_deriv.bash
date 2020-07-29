#!/bin/bash
#SBATCH --job-name=gam_deriv
#SBATCH --output=%x_%A_%a.out
#SBATCH --account=somerville_lab
#SBATCH --time=1-00:00:00
#SBATCH --partition=ncf
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
module load gcc/7.1.0-fasrc01
module load R/3.5.1-fasrc01
cp ~/.R/Makevars.gcc ~/.R/Makevars
export R_LIBS_USER=/ncf/mclaughlin/users/jflournoy/R_3.5.1_GCC:$R_LIBS_USER

cores=$SLURM_CPUS_ON_NODE

srun -c "${cores}" Rscript --verbose --vanilla gam_deriv_7nets.R