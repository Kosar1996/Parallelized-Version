#!/bin/bash
#SBATCH --job-name=w150_N8
#SBATCH --partition=pi_qmqi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/run_parallel_150_wide6_N8_%j.out
#SBATCH --error=logs/run_parallel_150_wide6_N8_%j.err

mkdir -p logs

module load matlab/matlab-2026a

export NUM_PROCS=8

echo "=================================================="
echo "Slurm Job ID: ${SLURM_JOB_ID}"
echo "Queue Timestamp (UTC): $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
echo "Script: run_parallel_150_wide6.m, NUM_PROCS=8"
echo "=================================================="

matlab -nodisplay -nosplash -nodesktop -r "try, run('run_parallel_150_wide6.m'), catch ME, disp(getReport(ME,'extended')); exit(1); end; exit(0);"
