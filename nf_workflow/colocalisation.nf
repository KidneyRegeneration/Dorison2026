#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define parameters
params.samplesheet = null
params.outdir = "results"

// Help message
def helpMessage() {
    log.info"""
    Usage:
    nextflow run main.nf --samplesheet samplesheet.csv --outdir results
    
    Required arguments:
    --samplesheet    Path to samplesheet CSV file with columns: input_folder, conditions, stain_name
    --outdir         Output directory (default: results)
    """.stripIndent()
}

// Show help message if no samplesheet provided
if (params.samplesheet == null) {
    helpMessage()
    exit 1
}


// Process 1: Find and filter .ims images
process FIND_IMS_IMAGES {
    tag "${input_folder}"
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "*.log"
    
    input:
    tuple val(input_folder), val(conditions), val(stain_name)
    
    output:
    path "found_images.txt", emit: images
    path "find_images_${stain_name}.log", emit: log
    val(stain_name), emit: stain_name
    script:
    conditions_list = conditions.split(',').collect { it.trim() }
    
    """
    #!/bin/bash
    
    # Create log file
    LOG_FILE="find_images_${stain_name}.log"
    echo "Processing folder: ${input_folder}" > \$LOG_FILE
    echo "Looking for conditions: ${conditions}" >> \$LOG_FILE
    echo "Stain name: ${stain_name}" >> \$LOG_FILE
    echo "---" >> \$LOG_FILE
    
    # Create output file for found images
    > found_images.txt
    
    # Find .ims files that don't contain "_20x_" in their path
    find "${input_folder}" -maxdepth 1 -name "*.ims" -type f | grep -v "_20x_" > temp_files.txt
    
    # Process each file
    while IFS= read -r filepath; do
        if [[ -f "\$filepath" ]]; then
            filename=\$(basename "\$filepath")
            
            # Check which condition this file belongs to
            condition_found=""
            for condition in ${conditions_list.join(' ')}; do
                if [[ "\$filename" == *"_\${condition}_"* ]]; then
                    condition_found="\$condition"
                    break
                fi
            done
            
            if [[ -n "\$condition_found" ]]; then
                date_timestamp=\$(echo "\$filename" | grep -o '^[0-9]\\{4\\}-[0-9]\\{2\\}-[0-9]\\{2\\}_[0-9]\\{2\\}\\.[0-9]\\{2\\}\\.[0-9]\\{2\\}')
                experiment=\$(echo "\$filename" | grep -o '_AD[0-9]\\{4\\}_' | sed 's/_//g')
                replicate=\$(echo "\$filename" | grep -o '_Rep[0-9]\\+_' | sed 's/_Rep\\([0-9]\\+\\)_/\\1/')
                echo "Found: condition=\$condition_found, experiment=\$experiment, replicate=\$replicate, date_timestamp=\$date_timestamp, file=\$filename" >> \$LOG_FILE
                echo "\$condition_found,\$experiment,\$replicate,\$date_timestamp,\$filename,\$filepath" >> found_images.txt
            else
                echo "No matching condition for: \$filename" >> \$LOG_FILE
            fi
        fi
    done < temp_files.txt
    
    rm -f temp_files.txt
    """
}

process FIND_NUCLEI {
    tag "${metadata.condition}_${metadata.experiment}_Rep${metadata.replicate}"
    publishDir "${params.publishDir}/NUCLEI_MASK/", mode: 'copy', pattern: "*.tif"
    publishDir "${params.publishDir}/NUCLEI_QC/", mode: 'copy', pattern: "*.png"
    container params.singularity_image

    input:
    tuple val(metadata), val(image_info), val(unique_id)
    
    output:
    tuple val(metadata), val(image_info), val(unique_id), path("NUCLEI_MASK_${unique_id}.tif"), path("NUCLEI_QC_${unique_id}.png"), emit: nuclei_mask    
    
    script:
    def (filename, filepath) = image_info
    """
    python /group/kidn4/ACTIVE/COLOCALISATION_PIPELINE/general_coloc/segmentation/find_nuclei.py -i "${filepath}" \
    -o NUCLEI_MASK_${unique_id}.tif \
    -q NUCLEI_QC_${unique_id}.png

    """
}

process FIND_GLOMS {
    tag "${metadata.condition}_${metadata.experiment}_Rep${metadata.replicate}"
    publishDir "${params.publishDir}/GLOM_MASK/", mode: 'copy'
    container params.singularity_image

    input:
    tuple val(metadata), val(image_info), val(unique_id)
    
    output:
    tuple val(metadata), val(image_info), val(unique_id), path("GLOM_MASK_${unique_id}.tif"),  emit: glom_mask    
    
    script:
    def (filename, filepath) = image_info
    """
    python /group/kidn4/ACTIVE/COLOCALISATION_PIPELINE/general_coloc/segmentation/find_gloms.py -i "${filepath}" \
    -o GLOM_MASK_${unique_id}.tif \

    """
}

process QC_GLOMS {
    tag "${metadata.condition}_${metadata.experiment}_Rep${metadata.replicate}"
    publishDir "${params.publishDir}/GLOM_MASK_QC/", mode: 'copy'
    container params.singularity_image

    input:
    tuple val(metadata), val(image_info), val(unique_id), val(glom_mask)
    
    output:
    tuple val(metadata), val(image_info), val(unique_id), path("GLOM_QC_${unique_id}.gif"),  emit: glom_qc    
    
    script:
    def (filename, filepath) = image_info
    """
    python /group/kidn4/ACTIVE/COLOCALISATION_PIPELINE/general_coloc/qc/qc_glom_masks.py -i "${filepath}" \
    -m ${glom_mask} \
    -o GLOM_QC_${unique_id}.gif \

    """

}

process FIND_PODOCYTES {
    tag "${metadata.condition}_${metadata.experiment}_Rep${metadata.replicate}"
    publishDir "${params.publishDir}/POD_MASK/", mode: 'copy'
    container params.singularity_image

    input:
    tuple val(metadata), val(image_info), val(unique_id), val(glom_mask), val(nuclei_mask)

    output:
    tuple val(metadata), val(image_info), val(unique_id), val(nuclei_mask), path("PODOCYTE_MASK_${unique_id}.tif"), emit: pod_mask
    
    script:
    """
    python /group/kidn4/ACTIVE/COLOCALISATION_PIPELINE/general_coloc/segmentation/find_podocytes.py \
    -n ${glom_mask} \
    -d ${nuclei_mask} \
    -o PODOCYTE_MASK_${unique_id}.tif

    """
}

process QUANTIFY_PODOCYTE_STAINS {
    tag "${metadata.condition}_${metadata.experiment}_Rep${metadata.replicate}"
    publishDir "${params.publishDir}/${stain_name}/", mode: 'copy'
    container params.singularity_image

    input:
    tuple val(metadata), val(image_info), val(unique_id), val(nuclei_mask),  val(pod_mask), val(stain_name)

    output:
    tuple val(metadata), val(image_info), val(unique_id), path("${stain_name}_MASK_${unique_id}.tif"), path("${stain_name}_${unique_id}.csv"), path("${stain_name}_${unique_id}.h5"), emit: pod_quantify
    
    script:
    def (filename, filepath) = image_info

    """
    python /group/kidn4/ACTIVE/COLOCALISATION_PIPELINE/general_coloc/quantify/get_stain_in_pods.py \
    -i "${filepath}" \
    -p ${pod_mask} \
    -d ${nuclei_mask} \
    -o ${stain_name}_MASK_${unique_id}.tif \
    -c ${stain_name}_${unique_id}.csv \
    -f ${stain_name}_${unique_id}.h5

    """
}

process GET_MEAN_CSVS {
    tag "CSV_COLLECT"
    publishDir "${params.publishDir}/CSV/", mode: 'copy'
    container params.singularity_image

    input:
    val(csvs)

    output:
    path "mean.csv", emit: mean_csv
    
    script:

    """
    python /group/kidn4/ACTIVE/COLOCALISATION_PIPELINE/general_coloc/quantify/combine_results.py \
    -i ${csvs.join(' ')} \
    -o mean.csv
    
    """

}



workflow {
    // Read samplesheet
    samplesheet_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, quote: '"')
        .map { row -> 
            tuple(row.input_folder, row.conditions, row.stain_name)
        }
    
    // Process 1: Find IMS images
    FIND_IMS_IMAGES(samplesheet_ch)
    
       
   // Convert the output file to a channel of tuples
    image_ch = FIND_IMS_IMAGES.out.images
    .splitCsv(header: false)
    .map { row ->
        // metadata as dictionary
        def metadata = [
            condition: row[0],
            experiment: row[1],
            replicate: row[2],
            date_timestamp: row[3]
        ]
        // image info: (filename, filepath) - now at positions 4,5
        def image_info = tuple(row[4], file(row[5]))
        // unique_id: concatenation including date_timestamp
        def unique_id = "${row[0]}_${row[1]}_${row[2]}_${row[3]}"
        
        tuple(metadata, image_info, unique_id)
    }
    stain_name = FIND_IMS_IMAGES.out.stain_name.unique()
    FIND_NUCLEI(image_ch)
    FIND_GLOMS(image_ch)
    
    QC_GLOMS(FIND_GLOMS.out.glom_mask)

    combined_masks = FIND_GLOMS.out.glom_mask
        .map {metadata, image_info, unique_id, glom_mask -> 
            tuple(metadata, image_info, unique_id,  glom_mask)
        }
        .join(
            FIND_NUCLEI.out.nuclei_mask
                .map { metadata, image_info, unique_id, nuclei_mask, nuclei_qc -> 
                    tuple(metadata, image_info, unique_id, nuclei_mask)
                }
        )
        .map { metadata, image_info, unique_id, glom_mask, image_info2, unique_id2, nuclei_mask ->
            tuple(metadata, image_info, unique_id, glom_mask, nuclei_mask)
        }

    FIND_PODOCYTES(combined_masks)
    
    stain_value = stain_name
    .unique()
    .first()

    pod_ch = FIND_PODOCYTES.out.pod_mask
    .map { metadata, image_info, unique_id, nuclei_mask,  pod_mask ->
        tuple(metadata, image_info, unique_id, nuclei_mask, pod_mask, stain_value.val)
    }

    
    QUANTIFY_PODOCYTE_STAINS(pod_ch)

    QUANTIFY_PODOCYTE_STAINS.out.pod_quantify

    csv_files = QUANTIFY_PODOCYTE_STAINS.out.pod_quantify
    .map { metadata, image_info, unique_id, stain_mask, output_csv, h5 ->
        output_csv
    }
    .collect()  // Collect all CSV files into a single list

    GET_MEAN_CSVS(csv_files)
}

// Workflow completion message
workflow.onComplete {
    log.info """
    Pipeline completed at: ${workflow.complete}
    Execution status: ${workflow.success ? 'OK' : 'failed'}
    Execution duration: ${workflow.duration}
    """
}