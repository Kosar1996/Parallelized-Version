#!/bin/bash
#SBATCH --job-name=par_N16
#SBATCH --partition=pi_qmqi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=logs/run_parallel_1timestep_N16_%j.out
#SBATCH --error=logs/run_parallel_1timestep_N16_%j.err

mkdir -p logs

module load matlab/matlab-2026a

export NUM_PROCS=16

echo "=================================================="
echo "Slurm Job ID: ${SLURM_JOB_ID}"
echo "Queue Timestamp (UTC): $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
echo "Script: run_parallel_1timestep.m, NUM_PROCS=16"
echo "=================================================="

matlab -nodisplay -nosplash -nodesktop -r "try, run('run_parallel_1timestep.m'), catch ME, disp(getReport(ME,'extended')); exit(1); end; exit(0);"
