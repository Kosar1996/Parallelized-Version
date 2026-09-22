#!/bin/bash
#SBATCH --job-name=7_s_parallel_heatmaps
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
module load ffmpeg 2>/dev/null || spack load ffmpeg 2>/dev/null || true

echo "=================================================="
echo "Slurm Job ID: ${SLURM_JOB_ID}"
echo "Queue Timestamp (UTC): $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
echo "Folder: Leukocyte_Main_Files-0918-7_s_long"
echo "Script: make_videos_parallel.m"
echo "=================================================="

# Pass SLURM cpus-per-task directly into MATLAB parallel pool
matlab -nodisplay -nosplash -nodesktop -r "parpool('local', ${SLURM_CPUS_PER_TASK}); run('make_videos_parallel.m');"

# Automatically convert AVI / MJ2 outputs to standard MP4 using native mpeg4 encoder
cd videos
for f in case_7_s_*.avi case_7_s_*.mj2; do
    if [ -f "$f" ]; then
        outName="${f%.*}.mp4"
        echo "Converting $f to MP4 via native mpeg4 codec..."
        
        ffmpeg -y -i "$f" \
            -c:v mpeg4 -q:v 2 \
            -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" \
            "${f%.*}_clean.mp4"
            
        if [ -f "${f%.*}_clean.mp4" ]; then
            mv "${f%.*}_clean.mp4" "$outName"
            rm -f "$f"
        fi
    fi
done
