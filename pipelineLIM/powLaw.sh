#!/bin/bash
#SBATCH --job-name=dof_grid_powLaw
#SBATCH --output="./pipelinelim/out/slurm-%j.out"
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
##SBATCH --mem=0

srun matlab -nosplash -nodesktop -nodisplay -r "\
cd('./cluster'); \
c=parcluster('local'); \
numWorkers=min(c.NumWorkers, str2double(getenv('SLURM_CPUS_PER_TASK'))); \
delete(gcp('nocreate')); \
parpool('local', numWorkers); \
gridsearch_dof_powlaw_save_cases(); \
exit;"