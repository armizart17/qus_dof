#!/bin/bash
#SBATCH --job-name=dof_grid_gaussian
#SBATCH --output=/mnt/nfs2/emiranda/proj/dof26/qus_dof/slurmout/slurm-%j.out
#SBATCH --partition=thinkstation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --nodelist=worker9

echo "=============================="
echo "SLURM JOB START"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "User: $(whoami)"
echo "PWD: $(pwd)"
echo "HOME: $HOME"
echo "Date: $(date)"
echo "=============================="

# Fix MATLAB home issues
export HOME=/mnt/nfs2/emiranda
export MATLAB_PREFDIR=/mnt/nfs2/emiranda/.matlab

mkdir -p $HOME
mkdir -p $MATLAB_PREFDIR
mkdir -p /mnt/nfs2/emiranda/proj/dof26/qus_dof/slurmout

echo "After setting HOME:"
echo "HOME: $HOME"
echo "MATLAB_PREFDIR: $MATLAB_PREFDIR"

echo "=============================="
echo "Starting MATLAB..."
echo "=============================="

srun matlab -nosplash -nodesktop -nodisplay -r "\
fprintf('--- MATLAB START ---\n'); \
fprintf('PWD before cd: %s\n', pwd); \
cd('/mnt/nfs2/emiranda/proj/dof26/qus_dof/cluster'); \
fprintf('PWD after cd: %s\n', pwd); \
fprintf('Starting parpool...\n'); \
c=parcluster('local'); \
numWorkers=min(c.NumWorkers, str2double(getenv('SLURM_CPUS_PER_TASK'))); \
fprintf('Workers requested: %d\n', numWorkers); \
delete(gcp('nocreate')); \
parpool('local', numWorkers); \
fprintf('Running main function...\n'); \
gridsearch_dof_gauss_save_cases(); \
fprintf('--- MATLAB END ---\n'); \
exit;"