#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define parameters
params.samplesheet = null
params.outdir = "results"
params.manual_glom_folder = null

// Help message
def helpMessage() {
    log.info"""
    Usage:
    nextflow run main.nf --samplesheet samplesheet.csv --outdir results [--manual_glom_folder /path/to/masks]
    
    Required arguments:
    --samplesheet         Path to samplesheet CSV file with columns: input_folder,conditions,markerA_name,markerA_idx,markerB_name,markerB_idx,glom_idx,nuclei_idx,markerA_threshold,markerA_smooth,markerB_threshold,markerB_smooth
    --outdir              Output directory (default: results)
    
    Optional arguments:
    --manual_glom_folder  Path to folder containing manually generated GLOM masks (GLOM_MASK_*.tif files).
                          If provided, masks will be loaded from this folder where available.
                          For samples without matching masks, automatic GLOM detection will be performed.
                          Masks should be named as: GLOM_MASK_<unique_id>.tif where unique_id matches sample identifiers.
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
    tuple val(input_folder), val(conditions), val(conditions2), val(stainA_name), val(stainA_index), val(stainB_name), val(stainB_index), val(glom_idx), val(nuclei_idx), val(thresholdA), val(smoothA), val(thresholdB), val(smoothB)
    
    output:
    path "found_images.txt", emit: images
    path "find_images_${stainA_name}_${stainB_name}.log", emit: log
    script:
    conditions_list = conditions.split(',').collect { it.trim() }
    conditions2_list = conditions2 ? conditions2.split(',').collect { it.trim() } : []
    
    """
    #!/bin/bash
    
    # Create log file
    LOG_FILE="find_images_${stainA_name}_${stainB_name}.log"
    echo "Processing folder: ${input_folder}" > \$LOG_FILE
    echo "Looking for conditions: ${conditions}" >> \$LOG_FILE
    ${conditions2_list ? "echo \"Looking for condition2: ${conditions2}\" >> \$LOG_FILE" : ""}
    echo "Markers: ${stainA_name}_${stainB_name}" >> \$LOG_FILE
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
                if [[ "\$filename" == *"_\${condition}_"* ]]; then
                    condition_found="\$condition"
                    break
                fi
            done
            
            # Check for condition2 if provided
            condition2_found=""
            ${conditions2_list ? """if [[ -n "\$condition_found" ]]; then
                for condition2 in ${conditions2_list.join(' ')}; do
                    if [[ "\$filename" == *"_\${condition2}_"* ]] || [[ "\$filename" == *"_\${condition2}." ]]; then
                        condition2_found="\$condition2"
                        break
                    fi
                done
            fi""" : ""}
            
            if [[ -n "\$condition_found" ]]; then
                date_timestamp=\$(echo "\$filename" | grep -o '^[0-9]\\{4\\}-[0-9]\\{2\\}-[0-9]\\{2\\}_[0-9]\\{2\\}\\.[0-9]\\{2\\}\\.[0-9]\\{2\\}')
                experiment=\$(echo "\$filename" | grep -o '_AD[0-9]\\{4\\}_' | head -n 1 | sed 's/_//g')
                if_exp=\$(echo "\$filename" | grep -oE 'IF[0-9]{2,4}_' | sed 's/_//g')
                replicate=\$(echo "\$filename" | grep -o '_[rR]ep[0-9]\\+_' | sed 's/_[rR]ep\\([0-9]\\+\\)_/\\1/')
                markerA_idx=${stainA_index}
                markerB_idx=${stainB_index}
                glom_idx=${glom_idx}
                nuclei_idx=${nuclei_idx}
                markerA_name=${stainA_name}
                markerB_name=${stainB_name}
                echo "Found: condition=\$condition_found, condition2=\$condition2_found, experiment=\$experiment, if_exp=\$if_exp, replicate=\$replicate, date_timestamp=\$date_timestamp,\${markerA_name}_idx=\$markerA_idx,\${markerB_name}_idx=\$markerB_idx,glom_idx=\$glom_idx,nuclei_idx=\$nuclei_idx file=\$filename threshold A = ${thresholdA}, smooth A = ${smoothA}, threshold B = ${thresholdB}, smooth B = ${smoothB}" >> \$LOG_FILE
                echo "\$condition_found,\$condition2_found,\$experiment,\$if_exp,\$replicate,\$date_timestamp,\$filename,\$filepath,\$markerA_name,\$markerB_name,\$markerA_idx,\$markerB_idx,\$glom_idx,\$nuclei_idx,${thresholdA},${smoothA},${thresholdB},${smoothB}" >> found_images.txt
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
    python3 /group/kidn4/ACTIVE/COLOCALISATION_PIPELINE/general_coloc/segmentation/find_gloms.py -i "${filepath}" \
    -o GLOM_MASK_${unique_id}.tif \
    --nphs1_index ${metadata.glom_idx} \
    -s 70 \
    -e \
    -g
    """
}

process PREPARE_MANUAL_GLOM_FILES {
    tag "${metadata.condition}_${metadata.experiment}_Rep${metadata.replicate}"
    publishDir "${params.publishDir}/INPUTS_FOR_MANUAL_MASKING/", mode: 'copy'
    container params.singularity_image

    input:
    tuple val(metadata), val(image_info), val(unique_id)
    
    output:
    tuple val(metadata), val(image_info), val(unique_id), path("${unique_id}.tif"),  emit: manual_mask_inputs    
    
    script:
    def (filename, filepath) = image_info
    """
    python3 /group/kidn4/ACTIVE/COLOCALISATION_PIPELINE/general_coloc/qc/subset_images.py -i "${filepath}" \
    -o ${unique_id}.tif \
    -c  ${metadata.nuclei_idx} ${metadata.glom_idx} \
    """
}


process LOAD_MANUAL_GLOMS {
    tag "${metadata.condition}_${metadata.experiment}_Rep${metadata.replicate}"
    publishDir "${params.publishDir}/GLOM_MASK_MANUAL/", mode: 'copy'

    input:
    tuple val(metadata), val(image_info), val(unique_id)

    output:
    tuple val(metadata), val(image_info), val(unique_id), path("GLOM_MASK_${unique_id}.tif"), optional: true, emit: glom_mask

    script:
    def manual_mask_name = "GLOM_MASK_${unique_id}.tif"
    """
    if [[ -f "${params.manual_glom_folder}/${manual_mask_name}" ]]; then
        cp "${params.manual_glom_folder}/${manual_mask_name}" "GLOM_MASK_${unique_id}.tif"
    fi
    """
}


process QC_GLOMS {
    tag "${metadata.condition}_${metadata.experiment}_Rep${metadata.replicate}"
    publishDir { "${params.publishDir}/${params.manual_glom_folder ? 'GLOM_QC_MANUAL' : 'GLOM_MASK_QC'}/" }, mode: 'copy'
    container params.singularity_image

    input:
    tuple val(metadata), val(image_info), val(unique_id), path(glom_mask)
    
    output:
    tuple val(metadata), val(image_info), val(unique_id), path("GLOM_QC_${unique_id}.gif"),  emit: glom_qc    
    
    script:
    def (filename, filepath) = image_info
    """
    python3 /group/kidn4/ACTIVE/COLOCALISATION_PIPELINE/general_coloc/qc/qc_glom_masks.py -i "${filepath}" \
    -m ${glom_mask} \
    -o GLOM_QC_${unique_id}.gif \
    -n ${metadata.glom_idx}
    """
}

process QC_MARKER_THRESHOLD {
    tag "${metadata.condition}_${metadata.experiment}_Rep${metadata.replicate}_${marker_name}"
    publishDir "${params.publishDir}/${marker_name}_QC/", mode: 'copy'
    container params.singularity_image

    input:
        tuple val(metadata), val(image_info), val(unique_id), val(marker_name),
              path(marker_mask)
    output:
        tuple val(metadata), val(image_info), val(unique_id), val(marker_name),
              path("${marker_name}_QC_${unique_id}.gif"), emit: marker_qc

    script:
    def (filename, filepath) = image_info
    def marker_idx = (marker_name == metadata.markerA_name) ? metadata.markerA_idx : metadata.markerB_idx
    """
    python3 /group/kidn4/ACTIVE/COLOCALISATION_PIPELINE/general_coloc/qc/qc_glom_masks.py \
        -i "${filepath}" \
        -m ${marker_mask} \
        -o ${marker_name}_QC_${unique_id}.gif \
        -n ${marker_idx}
    """
}

process THRESHOLD_MARKER {
    tag "${metadata.condition}_${metadata.experiment}_Rep${metadata.replicate}_${marker_name}"
    publishDir "${params.publishDir}/${marker_name}_THRESHOLD/", mode: 'copy'
    container params.singularity_image

    input:
        tuple val(metadata), val(image_info), val(unique_id),
              val(marker_name), val(marker_idx), path(glom_mask), val(threshold_method), val(smoothing_method)

    output:
        tuple val(metadata), val(image_info), val(unique_id), val(marker_name),
              path("${marker_name}_THRESHOLD_${unique_id}.tif"), emit: threshold_masks

    script:
    def (filename, filepath) = image_info
    """
    python3 /group/kidn4/ACTIVE/COLOCALISATION_PIPELINE/general_coloc/segmentation/threshold_marker.py \
        -i "${filepath}" \
        -c ${marker_idx} \
        -o ${marker_name}_THRESHOLD_${unique_id}.tif \
        --mask_apply_stage after \
        --rolling_ball_radius 51 \
        -t ${threshold_method} \
        -s ${smoothing_method}

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
    python3 /group/kidn4/ACTIVE/COLOCALISATION_PIPELINE/general_coloc/quantify/quantify_marker_in_gloms.py \
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

process CALCULATE_COLOCALISATION {
    tag "${metadata.condition}_${metadata.experiment}_Rep${metadata.replicate}"
    publishDir "${params.publishDir}/COLOCALISATION/", mode: 'copy'
    container params.singularity_image

    input:
        tuple val(metadata), val(image_info), val(unique_id),
              path(mask_a), path(mask_b), path(glom_mask)

    output:
        tuple val(metadata), val(image_info), val(unique_id),
              path("COLOC_${unique_id}.csv"),
              path("COLOC_${unique_id}_scatter.png"),
              path("COLOC_${unique_id}_snapshototsu.png"),
              emit: colocalisation

    script:
    def (filename, filepath) = image_info
    """
    python3 /group/kidn4/ACTIVE/COLOCALISATION_PIPELINE/general_coloc/colocalisation/calculate_colocalisation.py \
        -i "${filepath}" \
        --channels ${metadata.markerA_idx} ${metadata.markerB_idx} \
        --mask-a ${mask_a} \
        --mask-b ${mask_b} \
        -m   ${glom_mask} \
        --names  ${metadata.markerA_name} ${metadata.markerB_name} \
        --output-csv      COLOC_${unique_id}.csv \
        --output-plot     COLOC_${unique_id}_scatter.png \
        --output-snapshot COLOC_${unique_id}_snapshototsu.png \
        --thresh_a ${metadata.thresholdA} \
        --thresh_b ${metadata.thresholdB} \
        --smooth_a ${metadata.smoothA} \
        --smooth_b ${metadata.smoothB} \

    """
}

process WRITE_COLOC_METADATA {
    tag "${unique_id}"
    container params.singularity_image

    input:
        tuple val(metadata), val(image_info), val(unique_id),
              path(coloc_csv), path(scatter_png), path(snapshot_png)

    output:
        tuple path(coloc_csv), path("${unique_id}_metadata.json")

    script:
    def sample_name = image_info[0].replace('.ims', '').replace(' ', '_')
    """
    cat > ${unique_id}_metadata.json << EOF
    {
        "sample_name":   "${sample_name}",
        "condition":     "${metadata.condition}",
        "condition2":    "${metadata.condition2 ?: ''}",
        "experiment":    "${metadata.experiment}",
        "replicate":     "${metadata.replicate}",
        "if_exp":        "${metadata.if_exp}",
        "unique_id":     "${unique_id}",
        "markerA_name":  "${metadata.markerA_name}",
        "markerB_name":  "${metadata.markerB_name}"
    }
    EOF
    """
}


process COMBINE_COLOC_CSVS {
    tag "COMBINE_COLOC_CSVS"
    publishDir "${params.publishDir}/COLOCALISATION_CSV/", mode: 'copy'
    container params.singularity_image

    input:
        tuple path(csvs), path(meta_files)

    output:
        path "${csv_name}", emit: combined_csv

    script:
    csv_name = "${file(params.publishDir).name}.csv"
    """
    python3 /group/kidn4/ACTIVE/COLOCALISATION_PIPELINE/general_coloc/quantify/combine_csvs.py \
        -i ${csvs.join(' ')} \
        -m ${meta_files.join(' ')} \
        -o ${csv_name}
    """
}
process ORGANIZE_QC_BY_CONDITION {
    tag "ORGANIZE_QC_BY_CONDITION - ${condition}"
    
    input:
        tuple path(glom_qc_gifs), path(markerA_qc_gifs), path(markerB_qc_gifs), val(markerA_name), val(markerB_name), val(condition)
    
    output:
        tuple path("qc_by_condition"), val(markerA_name), val(markerB_name), val(condition), emit: qc_organized
    
    script:
    """
        mkdir -p qc_by_condition
        for f in ${glom_qc_gifs} ${markerA_qc_gifs} ${markerB_qc_gifs}; do
            [ -f "\$f" ] && ln -s \$(readlink -f \$f) qc_by_condition/ 2>/dev/null || true
        done
    """
}

process GENERATE_CONDITION_POWERPOINTS {
    tag "GENERATE_CONDITION_POWERPOINTS"
    publishDir "${params.publishDir}/POWERPOINT/", mode: 'copy'
    container params.singularity_image

    input:
    tuple path(qc_organized_folder),  val(markerA_name), val(markerB_name), val(condition)

    output:
    path "*.pptx", emit: powerpoints

    script:
    """
    python3 /group/kidn4/ACTIVE/COLOCALISATION_PIPELINE/general_coloc/qc/generate_qc_powerpoint.py \
        -i ${qc_organized_folder} \
        -o ${condition}_QC3.pptx \
        --marker-a ${markerA_name} \
        --marker-b ${markerB_name}
    """
}

process PLOT_COLOC_RESULTS {
    tag "PLOT_COLOC_RESULTS"
    publishDir "${params.publishDir}/COLOCALISATION_PLOTS/", mode: 'copy'
    container params.singularity_image

    input:
        path(combined_csv)
        tuple val(marker_a), val(marker_b)

    output:
        path("boxplots_${marker_a}_in_${marker_b}/*")
        path("boxplots_${marker_b}_in_${marker_a}/*")
        path("ridgeplots_${marker_a}_in_${marker_b}/*")
        path("ridgeplots_${marker_b}_in_${marker_a}/*")

    script:
    """
    python3 /group/kidn4/ACTIVE/COLOCALISATION_PIPELINE/general_coloc/plotting/plot_colocalisation.py \
        -i ${combined_csv} \
        --col-a-in-b MOC_${marker_a}in${marker_b} \
        --col-b-in-a MOC_${marker_b}in${marker_a} \
        --marker-a ${marker_a} \
        --marker-b ${marker_b} \
        --boxplot-a-in-b    boxplots_${marker_a}_in_${marker_b} \
        --boxplot-b-in-a    boxplots_${marker_b}_in_${marker_a} \
        --ridgeplot-a-in-b  ridgeplots_${marker_a}_in_${marker_b} \
        --ridgeplot-b-in-a  ridgeplots_${marker_b}_in_${marker_a} \
        --group-by condition \
        --hue-by condition2
    """
}


workflow {
    // Read samplesheet
    samplesheet_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, quote: '"')
        .map { row -> 
            // Handle optional condition2 column
            def conditions2 = row.condition2 ?: ""
            tuple(row.input_folder, row.conditions, conditions2, row.markerA_name, row.markerA_idx, row.markerB_name, row.markerB_idx, row.glom_idx, row.nuclei_idx, row.markerA_threshold, row.markerA_smooth, row.markerB_threshold, row.markerB_smooth)
        }
    
    // Process 1: Find IMS images
    FIND_IMS_IMAGES(samplesheet_ch)

    image_ch = FIND_IMS_IMAGES.out.images
        .splitCsv(header: false)
        .map { row ->
            // metadata as dictionary
            def metadata = [
                condition: row[0],
                condition2: row[1],
                experiment: row[2],
                if_exp: row[3],
                replicate: row[4],
                date_timestamp: row[5],
                markerA_name: row[8],
                markerB_name: row[9],
                markerA_idx: row[10],
                markerB_idx: row[11],
                glom_idx: row[12],
                nuclei_idx: row[13],
                thresholdA:  row[14],
                smoothA:  row[15],
                thresholdB:  row[16],
                smoothB:  row[17],
            ]
            // image info: (filename, filepath)
            def image_info = tuple(row[6], file(row[7]))
            // unique_id: concatenation including condition2 and date_timestamp
            def unique_id = "${row[0].replace('_', '-')}_${row[1].isEmpty() ? '' : row[1].replace('_', '-') + '_'}${row[2]}_${row[3]}${row[4]}_${row[5]}"
            tuple(metadata, image_info, unique_id)
        }

    PREPARE_MANUAL_GLOM_FILES(image_ch)

    // Process GLOM masks: load manual if available, otherwise auto-detect
    // Try to load manual masks first
    load_manual_ch = image_ch
    // 1. Try loading manual masks (optional emissions)
    if (params.manual_glom_folder) {
        LOAD_MANUAL_GLOMS(image_ch)
        manual_masks_ch = LOAD_MANUAL_GLOMS.out.glom_mask
    } else {
        manual_masks_ch = Channel.empty()
    }

    // 2. Create lookup maps keyed by unique_id

    manual_lookup = manual_masks_ch
        .map { metadata, image_info, unique_id, glom_mask ->
            tuple(unique_id, tuple(metadata, image_info, unique_id, glom_mask))
        }

    image_lookup = image_ch
        .map { metadata, image_info, unique_id ->
            tuple(unique_id, tuple(metadata, image_info, unique_id))
        }

    // 3. LEFT JOIN (keeps all samples)
    
    joined = image_lookup.join(manual_lookup, remainder: true)


    // 4. Split into WITH and WITHOUT manual masks

    with_manual_ch = joined
        .filter { uid, image_tuple, manual_tuple -> manual_tuple != null }
        .map { uid, image_tuple, manual_tuple -> manual_tuple }

    without_manual_ch = joined
        .filter { uid, image_tuple, manual_tuple -> manual_tuple == null }
        .map { uid, image_tuple, manual_tuple -> image_tuple }

    // 5. Run auto segmentation ONLY on missing
    FIND_GLOMS(without_manual_ch)
    auto_masks_ch = FIND_GLOMS.out.glom_mask

    // 6. Combine final channel
    glom_mask_ch = with_manual_ch.mix(auto_masks_ch)

    QC_GLOMS(glom_mask_ch)
    qc_glom_ch = QC_GLOMS.out.glom_qc

    image_with_glom_ch = image_ch
        .map { metadata, image_info, unique_id ->
            tuple(unique_id, metadata, image_info)
        }
        .join(
            glom_mask_ch.map { metadata, image_info, unique_id, glom_mask ->
                tuple(unique_id, glom_mask)
            }
        )
        .flatMap { unique_id, metadata, image_info, glom_mask ->
            [
                tuple(metadata, image_info, unique_id, metadata.markerA_name, metadata.markerA_idx, glom_mask, metadata.thresholdA, metadata.smoothA),
                tuple(metadata, image_info, unique_id, metadata.markerB_name, metadata.markerB_idx, glom_mask, metadata.thresholdB, metadata.smoothB)
            ]
        }

    THRESHOLD_MARKER(image_with_glom_ch)

    QC_MARKER_THRESHOLD(THRESHOLD_MARKER.out.threshold_masks)


    threshold_a_ch = THRESHOLD_MARKER.out.threshold_masks
        .filter { metadata, image_info, unique_id, marker_name, mask -> 
            marker_name == metadata.markerA_name 
        }
        .map { metadata, image_info, unique_id, marker_name, mask ->
            tuple(unique_id, metadata, image_info, mask)
        }

    threshold_b_ch = THRESHOLD_MARKER.out.threshold_masks
        .filter { metadata, image_info, unique_id, marker_name, mask -> 
            marker_name == metadata.markerB_name 
        }
        .map { metadata, image_info, unique_id, marker_name, mask ->
            tuple(unique_id, mask)
        }

    // Join A and B masks on unique_id, then bring in glom mask too
    coloc_input_ch = threshold_a_ch
        .join(threshold_b_ch)                           // keyed on unique_id -> (uid, meta, image_info, maskA, maskB)
        .join(                                           // bring in glom mask
            glom_mask_ch.map { metadata, image_info, unique_id, glom_mask ->
                tuple(unique_id, glom_mask)
            }
        )
        .map { unique_id, metadata, image_info, mask_a, mask_b, glom_mask ->
            tuple(metadata, image_info, unique_id, mask_a, mask_b, glom_mask)
        }

    CALCULATE_COLOCALISATION(coloc_input_ch)



    coloc_with_meta_ch = CALCULATE_COLOCALISATION.out.colocalisation
    marker_names_ch = CALCULATE_COLOCALISATION.out.colocalisation
        .first()
        .map { metadata, image_info, unique_id, csv, scatter, snapshot ->
            tuple(metadata.markerA_name, metadata.markerB_name)
        }

 

    csv_files = WRITE_COLOC_METADATA(coloc_with_meta_ch)
        .collect()
        .map { flat ->
            def pairs      = flat.collate(2)
            def csvs       = pairs.collect { it[0] }
            def meta_files = pairs.collect { it[1] }
            tuple(csvs, meta_files)
        }

    COMBINE_COLOC_CSVS(csv_files)

    // Extract marker names
    markerA_name = marker_names_ch.map { a, b -> a }
    markerB_name = marker_names_ch.map { a, b -> b }

    conditions_list = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, quote: '"')
        .map { row -> row.conditions }
        .unique()
        .collect()


    glom_qc_collected = qc_glom_ch
    .map { metadata, image_info, unique_id, gif -> [metadata.condition, gif] }
    .groupTuple()

markerA_qc_collected = QC_MARKER_THRESHOLD.out.marker_qc
    .combine(markerA_name)
    .filter { metadata, image_info, unique_id, marker_name, gif, markerA ->
        marker_name?.trim() == markerA?.trim()
    }
    .map { metadata, image_info, unique_id, marker_name, gif, markerA -> [metadata.condition, gif] }
    .groupTuple()

markerB_qc_collected = QC_MARKER_THRESHOLD.out.marker_qc
    .combine(markerB_name)
    .filter { metadata, image_info, unique_id, marker_name, gif, markerB ->
        marker_name?.trim() == markerB?.trim()
    }
    .map { metadata, image_info, unique_id, marker_name, gif, markerB -> [metadata.condition, gif] }
    .groupTuple()

    // Join by condition
    qc_inputs = glom_qc_collected
        .join(markerA_qc_collected)
        .join(markerB_qc_collected)
        .map { condition, glom_gifs, markerA_gifs, markerB_gifs ->
            [glom_gifs, markerA_gifs, markerB_gifs, markerA_name.getVal(), markerB_name.getVal(), condition]
        }
    // qc_inputs.view()
 
    // Organize QC files by condition
    ORGANIZE_QC_BY_CONDITION(
       qc_inputs
    )
    

    GENERATE_CONDITION_POWERPOINTS(
        ORGANIZE_QC_BY_CONDITION.out.qc_organized,
    )


    PLOT_COLOC_RESULTS(
        COMBINE_COLOC_CSVS.out.combined_csv,
        marker_names_ch
    )


}

// Workflow completion message
workflow.onComplete {
    log.info """
    Pipeline completed at: ${workflow.complete}
    Execution status: ${workflow.success ? 'OK' : 'failed'}
    Execution duration: ${workflow.duration}
    """
}
