#!/bin/bash
#SBATCH --job-name=p0930_verify
#SBATCH --partition=pi_qmqi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --output=logs/verify_parfor_%j.out
#SBATCH --error=logs/verify_parfor_%j.err

mkdir -p logs
module load matlab/matlab-2026a

export NUM_PROCS=0
export DO_PROFILE=0
export RUN_TAG=verify_parfor

echo "=================================================="
echo "Slurm Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURMD_NODENAME}  CPUs: ${SLURM_CPUS_PER_TASK}  Mem: 16G  Time cap: 00:45:00"
echo "Queue Timestamp (UTC): $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
echo "Folder: $(basename $(pwd))   Script: verify_parfor_0930.m"
echo "RUN_TAG=verify_parfor  NUM_PROCS=0  DO_PROFILE=0"
echo "=================================================="

matlab -nodisplay -nosplash -nodesktop -r "try, run('verify_parfor_0930.m'), catch ME, disp(getReport(ME,'extended')); exit(1); end; exit(0);"
