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
    --samplesheet    Path to samplesheet CSV file with columns: input_folder,conditions,marker_name,marker_idx,glom_idx,nuclei_idx
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
    publishDir "${params.publishDir}/FOUND_IMAGES_LOG/", mode: 'copy', pattern: "*.log"

    input:
    tuple val(input_folder), val(conditions), val(stain_name), val(channel_idx), val(glom_idx), val(nuclei_idx)
    
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
    find "${input_folder}" -maxdepth 1 -name "2*.ims" -type f | grep -v "_20x" | grep -v "ARTIFACT" > temp_files.txt

    # Process each file
    while IFS= read -r filepath; do
        if [[ -f "\$filepath" ]]; then
            filename=\$(basename "\$filepath")
            
            # Check which condition this file belongs to
            condition_found=""
            for condition in ${conditions_list.join(' ')}; do
                if [[ "\$filename" == *"_\${condition}"* ]]; then
                    condition_found="\$condition"
                    break
                fi
            done
            
            if [[ -n "\$condition_found" ]]; then
                date_timestamp=\$(echo "\$filename" | grep -o '^[0-9]\\{4\\}-[0-9]\\{2\\}-[0-9]\\{2\\}_[0-9]\\{2\\}\\.[0-9]\\{2\\}\\.[0-9]\\{2\\}')
                experiment=\$(echo "\$filename" | grep -o '_AD[0-9]\\{4\\}_' | sed 's/_//g')
                if_exp=\$(echo "\$filename" | grep -oE '_AD(SY)?IF[0-9]{2,4}_' | sed 's/_//g')
                replicate=\$(echo "\$filename" | grep -o '_[rR]ep[0-9]\\+[_.]' | sed 's/_[rR]ep\\([0-9]\\+\\)[_.]/\\1/')
                channel_idx=${channel_idx}
                glom_idx=${glom_idx}
                nuclei_idx=${nuclei_idx}
                echo "Found: condition=\$condition_found, experiment=\$experiment, if_exp=\$if_exp, replicate=\$replicate, date_timestamp=\$date_timestamp, channel_idx=\$channel_idx, glom_idx=\$glom_idx,nuceli_idx=\$nuclei_idx,file=\$filename" >> \$LOG_FILE
                echo "\$condition_found,\$experiment,\$if_exp,\$replicate,\$date_timestamp,\$filename,\$filepath",\$channel_idx,\$glom_idx,\$nuclei_idx >> found_images.txt
            else
                echo "No matching condition for: \$filename" >> \$LOG_FILE
            fi
        fi
    done < temp_files.txt
    
    rm -f temp_files.txt
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
    --nphs1_index ${metadata.glom_idx} 
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
    -n ${metadata.glom_idx}
    """

}

process QUANTIFY_GLOM_STAINS {
    tag "${metadata.condition}_${metadata.experiment}_Rep${metadata.replicate}"
    publishDir "${params.publishDir}/${stain_name}/", mode: 'copy'
    container params.singularity_image

    input:
    tuple val(metadata), val(image_info), val(unique_id),  val(glom_mask), val(stain_name)

    output:
    tuple val(metadata), val(image_info), val(unique_id), val(stain_name), path("${stain_name}_MASK_${unique_id}.tif"), path("${stain_name}_${unique_id}.csv"), path("${stain_name}_${unique_id}.h5"), emit: glom_quantify
    
    script:
    def (filename, filepath) = image_info

    """
    python /group/kidn4/ACTIVE/COLOCALISATION_PIPELINE/general_coloc/quantify/quantify_marker_in_gloms.py \
    -i "${filepath}" \
    -g ${glom_mask} \
    -o ${stain_name}_MASK_${unique_id}.tif \
    -c ${stain_name}_${unique_id}.csv \
    -f ${stain_name}_${unique_id}.h5 \
    -n ${metadata.marker_idx} \
    --stain_names ${stain_name} \
    -x ${metadata.nuclei_idx}

    """
}



process GET_MEAN_CSVS {
    tag "CSV_COLLECT"
    publishDir "${params.publishDir}/CSV/", mode: 'copy'
    container params.singularity_image

    input:
        tuple path(csvs), path(meta_files)

    output:
        path "mean.csv", emit: mean_csv

    script:
    """
    python /group/kidn4/ACTIVE/COLOCALISATION_PIPELINE/general_coloc/quantify/combine_csvs.py \
        -i ${csvs.join(' ')} \
        -m ${meta_files.join(' ')} \
        -o mean.csv
    """
}

process PLOT_RESULTS{
    tag "PLOT_RESULTS"
    publishDir "${params.publishDir}/PLOTS/", mode: 'copy'
    container params.singularity_image

    input: 
    val(mean_csv)

    output:
    path("boxplots/*")
    path("ridgeplots/*")

    script:
    """
    python /group/kidn4/ACTIVE/COLOCALISATION_PIPELINE/general_coloc/plotting/plot_in_gloms.py \
    -i ${mean_csv} \
    -b boxplots \
    -r ridgeplots

    """
}

process WRITE_METADATA {
    tag "${unique_id}"
    container params.singularity_image

    input:
        tuple val(metadata), val(image_info), val(unique_id), val(stain_name), path(stain_mask), path(output_csv), path(h5)

    output:
        tuple path(output_csv), path("${unique_id}_metadata.json")

    script:
    def sample_name = image_info[0].replace('.ims', '').replace(' ', '_')
    """
    cat > ${unique_id}_metadata.json << EOF
    {
        "sample_name":  "${sample_name}",
        "stain":        "${stain_name}",
        "condition":    "${metadata.condition}",
        "experiment":   "${metadata.experiment}",
        "replicate":    "${metadata.replicate}",
        "if_exp":       "${metadata.if_exp}",
        "unique_id":    "${unique_id}"
    }
    EOF
    """
}


workflow {
    // Read samplesheet
    samplesheet_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, quote: '"')
        .map { row -> 
            tuple(row.input_folder, row.conditions, row.marker_name, row.marker_idx, row.glom_idx, row.nuclei_idx)
        }
    
    // Process 1: Find IMS images
    FIND_IMS_IMAGES(samplesheet_ch)
    
    image_ch = FIND_IMS_IMAGES.out.images
        .splitCsv(header: false)
        .map { row ->
            // metadata as dictionary
            def metadata = [
                condition: row[0],
                experiment: row[1],
                if_exp: row[2],
                replicate: row[3],
                date_timestamp: row[4],
                marker_idx: row[7],
                glom_idx: row[8],
                nuclei_idx: row[9]
            ]
            // image info: (filename, filepath) - now at positions 4,5
            def image_info = tuple(row[5], file(row[6]))
            // unique_id: concatenation including date_timestamp
            def unique_id = "${row[0].replace('_', '-')}_${row[1]}_${row[2]}${row[3]}_${row[4]}"
            tuple(metadata, image_info, unique_id)
        }


    stain_name = FIND_IMS_IMAGES.out.stain_name.unique()
    
    stain_value = stain_name
    .unique()
    .first()

    FIND_GLOMS(image_ch)
    
    QC_GLOMS(FIND_GLOMS.out.glom_mask)
    
    glom_ch = FIND_GLOMS.out.glom_mask
    .map { metadata, image_info, unique_id,  glom_mask ->
        tuple(metadata, image_info, unique_id, glom_mask, stain_value.val)
    }
    
    QUANTIFY_GLOM_STAINS(glom_ch)


    // csv_files = QUANTIFY_GLOM_STAINS.out.glom_quantify
    // .map { metadata, image_info, unique_id, stain_name, stain_mask, output_csv, h5 ->
    //     output_csv
    // }
    // .collect()  // Collect all CSV files into a single list

    csv_files = WRITE_METADATA(QUANTIFY_GLOM_STAINS.out.glom_quantify)
    .collect()
    .map { flat ->
        def pairs        = flat.collate(2)
        def csvs         = pairs.collect { it[0] }
        def meta_files   = pairs.collect { it[1] }
        tuple(csvs, meta_files)
    }



    GET_MEAN_CSVS(csv_files)

    PLOT_RESULTS(GET_MEAN_CSVS.out.mean_csv)

}

// Workflow completion message
workflow.onComplete {
    log.info """
    Pipeline completed at: ${workflow.complete}
    Execution status: ${workflow.success ? 'OK' : 'failed'}
    Execution duration: ${workflow.duration}
    """
}