#!/bin/bash
#SBATCH --job-name=p0930_prof3
#SBATCH --partition=pi_qmqi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=logs/par_prof3_N8_%j.out
#SBATCH --error=logs/par_prof3_N8_%j.err

mkdir -p logs
module load matlab/matlab-2026a

export NUM_PROCS=8
export DO_PROFILE=1
export RUN_TAG=par_prof3_N8

echo "=================================================="
echo "Slurm Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURMD_NODENAME}  CPUs: ${SLURM_CPUS_PER_TASK}  Mem: 64G  Time cap: 14:00:00"
echo "Queue Timestamp (UTC): $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
echo "Folder: $(basename $(pwd))   Script: run_0930_prof3.m"
echo "RUN_TAG=par_prof3_N8  NUM_PROCS=8  DO_PROFILE=1"
echo "=================================================="

matlab -nodisplay -nosplash -nodesktop -r "try, run('run_0930_prof3.m'), catch ME, disp(getReport(ME,'extended')); exit(1); end; exit(0);"
