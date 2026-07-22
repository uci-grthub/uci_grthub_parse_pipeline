#!/bin/bash
#SBATCH --job-name=geo_ftp_upload
#SBATCH -A SBSANDME_LAB
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --error=logs/ftp/slurm-%x.%j.err
#SBATCH --output=logs/ftp/slurm-%x.%j.out

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

python3 -c "import yaml" 2>/dev/null || pip install --user pyyaml

python3 src/upload_geo_ftp.py "$@"
