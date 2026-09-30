#!/bin/bash
#SBATCH --job-name=p0930_N4
#SBATCH --partition=pi_qmqi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=10:00:00
#SBATCH --output=logs/par_N4_%j.out
#SBATCH --error=logs/par_N4_%j.err

mkdir -p logs
module load matlab/matlab-2026a

export NUM_PROCS=4
export DO_PROFILE=0
export RUN_TAG=par_N4

echo "=================================================="
echo "Slurm Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURMD_NODENAME}  CPUs: ${SLURM_CPUS_PER_TASK}  Mem: 48G  Time cap: 10:00:00"
echo "Queue Timestamp (UTC): $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
echo "Folder: $(basename $(pwd))   Script: run_0930_bench.m"
echo "RUN_TAG=par_N4  NUM_PROCS=4  DO_PROFILE=0"
echo "=================================================="

matlab -nodisplay -nosplash -nodesktop -r "try, run('run_0930_bench.m'), catch ME, disp(getReport(ME,'extended')); exit(1); end; exit(0);"
