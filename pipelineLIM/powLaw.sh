#!/bin/bash
#SBATCH --job-name=dof_grid_powLaw
##SBATCH --output="./slurmout/slurm-%j.out"
##SBATCH --partition=thinkstation
##SBATCH --nodelist=worker7
##SBATCH --nodes=1
##SBATCH --ntasks=1
##SBATCH --cpus-per-task=10
##SBATCH --mem=0

#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --nodelist=worker4

srun matlab -nosplash -nodesktop -nodisplay -r "\
cd('./cluster'); \
c=parcluster('local'); \
numWorkers=min(c.NumWorkers, str2double(getenv('SLURM_CPUS_PER_TASK'))); \
delete(gcp('nocreate')); \
parpool('local', numWorkers); \
gridsearch_dof_powlaw_save_cases(); \
exit;"