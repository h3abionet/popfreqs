#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * STEP 1 - Get chip site to bed
 */
process process_chip_list {
    tag "process_chip_list_${chip_name}"
    publishDir "${params.outdir}/chip_list", mode: 'copy'

    input:
        tuple val(chip_name), file(chip_file) from chips

    output:
        tuple val(chip_name), file(chip_file), file("${out_chip_file}*") into process_chip_list

    script:
    out_chip_file = "${chip_file.baseName}"
    template "process_chip_list.py"
}