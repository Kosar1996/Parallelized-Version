#!/bin/bash
#SBATCH --job-name=plot_0923
#SBATCH --partition=pi_qmqi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=01:00:00
#SBATCH --output=logs/plot_parallel_%j.out
#SBATCH --error=logs/plot_parallel_%j.err

mkdir -p logs

module load matlab/matlab-2026a
module load ffmpeg 2>/dev/null || spack load ffmpeg 2>/dev/null || true

echo "=================================================="
echo "Slurm Job ID: ${SLURM_JOB_ID}"
echo "Queue Timestamp (UTC): $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
echo "Folder: Leukocyte_Main_Files-0923"
echo "Script: make_videos_parallel.m"
echo "Cores: ${SLURM_CPUS_PER_TASK} | Memory: 128G"
echo "=================================================="

# 1. Execute MATLAB parallel pool rendering job
matlab -nodisplay -nosplash -nodesktop -r "parpool('local', ${SLURM_CPUS_PER_TASK}); run('make_videos_parallel.m');"

# 2. Bash Post-Processing: Convert AVI to MP4 using built-in native MPEG-4 encoder
VIDEO_DIR="videos"

if [ -d "$VIDEO_DIR" ]; then
    echo "=================================================="
    echo "Checking for AVI files to convert in ${VIDEO_DIR}..."
    echo "=================================================="
    
    for avi_file in "${VIDEO_DIR}"/*.avi; do
        if [ -f "$avi_file" ]; then
            mp4_file="${avi_file%.avi}.mp4"
            echo "Converting: ${avi_file} -> ${mp4_file}"
            
            # Using native -c:v mpeg4 and high target bitrate -b:v 8M (no libx264 dependency)
            ffmpeg -y -i "${avi_file}" \
                -c:v mpeg4 -b:v 8M \
                -pix_fmt yuv420p \
                -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" \
                "${mp4_file}"
            
            if [ -f "${mp4_file}" ]; then
                echo "Successfully created ${mp4_file}. Removing original ${avi_file}."
                rm -f "${avi_file}"
            else
                echo "[Warning] Conversion failed for ${avi_file}."
            fi
        fi
    done
fi

echo "=================================================="
echo "All post-processing tasks completed."
echo "=================================================="
