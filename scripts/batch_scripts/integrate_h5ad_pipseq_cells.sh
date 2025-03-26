#!/bin/bash
#SBATCH --job-name=pip_cells_integrate_h5ad
#SBATCH --output=integrate_h5ad_%j.out
#SBATCH --error=integrate_h5ad_%j.err
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
#SBATCH --partition=short

# Print some job information
echo "Starting at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Job name: $SLURM_JOB_NAME"
echo "Node: $SLURMD_NODENAME"


# Set input and output paths (modify these as needed)
INPUT_DIR="/private/groups/russelllab/jodie/scRNAseq/pipseq/cell_pipseq/h5ad_output"
OUTPUT_FILE="/private/groups/russelllab/jodie/scRNAseq/pipseq/cell_pipseq/integrated_data.h5ad"
PYTHON_SCRIPTS="/private/groups/russelllab/jodie/scRNAseq/scripts"

# Run the integration script with BBKNN and Harmony
python $PYTHON_SCRIPTS/integrate_h5ad.py \
    --input_dir "$INPUT_DIR" \
    --output_file "$OUTPUT_FILE" \
    --min_cells 3 \
    --min_genes 200 \
    --n_pcs 50 \
    --n_neighbors 20 \
    --method "both"  \
    --calculate_titer

# Print completion message
echo "Job completed at: $(date)"