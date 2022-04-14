#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { check_files } from './modules/qc'
include { parse_manifests; generate_vcf_from_bed; snpeff_vcf } from './modules/annotate'


"""
Read a manifest file
Generate
    - VCF file
    - TSV of selected fields for annotating the VCF
"""

workflow{
    manifests = Channel.from(params.manifests)
        .map{ name, manifest, sep, coordinate_fields, annotations ->
            check_files([ manifest ])
            annots = []
            descs = []
            types = []
            annotations.each{ annot_field, annot_desc, annot_type ->
                    annots << annot_field
                    descs << annot_desc
                    types << annot_type
            }
            return [ name, file(manifest), coordinate_fields, annots.join(','), descs.join(','), types.join(','), sep ]
        }

    // //// Generate BED from manifest
    parse_manifests(manifests)

    if(generate_vcf == true){
        // //// Generate VCF from BED
        check_files([params.human_genome])
        beds = parse_manifests.out.map{ manifest, manifest_annots, manifest_hdr, manifest_bed -> [ manifest, file(manifest_bed), file(params.human_genome) ] }
        generate_vcf_from_bed( beds )
    }

    println "Cuurent folder: ${PWD}, Output folder: ${params.outdir}"
}