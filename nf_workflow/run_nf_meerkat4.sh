#!/bin/bash
#SBATCH --partition=prod_med
#SBATCH --mem=2G
#SBATCH --time=24:00:00

##### EDIT HERE: ######
EXPERIMENT_NAME="MY_EXPERIMENT"
SAMPLESHEET="$(dirname "$0")/example_samplesheet.csv"

###################

PUBLISHDIR="/group/kidn4/ACTIVE/COLOCALISATION_PIPELINE/PIPELINE_OUTPUTS/${EXPERIMENT_NAME}"
OUTPUTDIR="/ibm/hpcfs1/tmp/${USER}/COLOCALISATION_NF_WORK/${EXPERIMENT_NAME}_results"

mkdir -p $PUBLISHDIR
module load apptainer
module load miniconda3

source /group/kidn3/GRAFT_REPOSITORY/scripts/graft_annotation_pipeline/conda_setup.sh
conda activate /group/kidn4/ACTIVE/PODOCIN_SCREEN_ANALYSIS/venvs/miniconda_m4
cd /ibm/hpcfs1/tmp/${USER}/COLOCALISATION_NF_WORK


nextflow run "$(dirname "$0")/colocalisation.nf" \
-profile meerkatv4 -w /ibm/hpcfs1/tmp/${USER}/COLOCALISATION_NF_WORK \
-o $OUTPUTDIR \
--samplesheet $SAMPLESHEET \
--publishDir $PUBLISHDIR \
--singularity_image /group/kidn4/ACTIVE/COLOCALISATION_PIPELINE/containers/colocalisation.sif \
-resume
