#!/bin/bash
#SBATCH --job-name=profile_5hr
#SBATCH --partition=pi_qmqi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=05:00:00
#SBATCH --output=logs/run_profile_1timestep_5hr_%j.out
#SBATCH --error=logs/run_profile_1timestep_5hr_%j.err

mkdir -p logs

module load matlab/matlab-2026a

echo "=================================================="
echo "Slurm Job ID: ${SLURM_JOB_ID}"
echo "Queue Timestamp (UTC): $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
echo "Script: run_profile_1timestep_5hr.m"
echo "=================================================="

matlab -nodisplay -nosplash -nodesktop -r "try, run('run_profile_1timestep_5hr.m'), catch ME, disp(getReport(ME,'extended')); exit(1); end; exit(0);"
