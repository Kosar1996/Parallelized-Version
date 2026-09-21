#!/bin/bash
#SBATCH --job-name=7_q_restart_parallel_heatmaps
#SBATCH --partition=pi_qmqi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=logs/plot_parallel_%j.out
#SBATCH --error=logs/plot_parallel_%j.err

mkdir -p logs

module load matlab/matlab-2026a

echo "=================================================="
echo "Slurm Job ID: ${SLURM_JOB_ID}"
echo "Queue Timestamp (UTC): $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
echo "Folder: Leukocyte_Main_Files-0918-7_q_restart"
echo "Script: make_videos_parallel.m"
echo "=================================================="

# Pass SLURM cpus-per-task directly into MATLAB parallel pool
matlab -nodisplay -nosplash -nodesktop -r "parpool('local', ${SLURM_CPUS_PER_TASK}); run('make_videos_parallel.m');"
