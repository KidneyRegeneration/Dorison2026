#!/bin/bash
#SBATCH --partition=prod_med
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --output=slurm_logs/slurm-%j.out

##### EDIT HERE: ######
EXPERIMENT_NAME="MY_EXPERIMENT"
SAMPLESHEET="$(dirname "$0")/samplesheet.csv"
WORKDIR="/ibm/hpcfs1/tmp/${USER}/${EXPERIMENT_NAME}"

###################
PUBLISHDIR="$(dirname "$0")/PIPELINE_OUTPUTS/${EXPERIMENT_NAME}"
OUTPUTDIR="/ibm/hpcfs1/tmp/${USER}/${EXPERIMENT_NAME}_results"
echo $PUBLISHDIR

mkdir -p "${PUBLISHDIR}"
cp $SAMPLESHEET $PUBLISHDIR
mkdir -p "${WORKDIR}"
module load apptainer
module load miniconda3

source /group/kidn3/GRAFT_REPOSITORY/scripts/graft_annotation_pipeline/conda_setup.sh
conda activate /group/kidn4/ACTIVE/PODOCIN_SCREEN_ANALYSIS/venvs/miniconda_m4
cd $WORKDIR

current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
nextflow run "$(dirname "$0")/pipeline.nf" \
-profile meerkatv4 -w "${WORKDIR}" \
-o "${OUTPUTDIR}" \
--samplesheet "${SAMPLESHEET}" \
--publishDir "${PUBLISHDIR}" \
--singularity_image /group/kidn4/ACTIVE/COLOCALISATION_PIPELINE/containers/colocalisation.sif \
-with-report "$(dirname "$0")/logs/${EXPERIMENT_NAME}-${current_datetime}.html" \
-with-trace "$(dirname "$0")/logs/${EXPERIMENT_NAME}-${current_datetime}.trace.txt" \
-with-timeline "$(dirname "$0")/logs/${EXPERIMENT_NAME}-${current_datetime}.timeline.html" \
-with-dag "$(dirname "$0")/logs/${EXPERIMENT_NAME}-${current_datetime}.dag.html" \
-resume

chmod 777 -R ${WORKDIR}
