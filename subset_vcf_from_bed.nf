#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { check_files; check_params } from './modules/qc'
include { get_bed_vcf; concat_vcf_chrms; extract_fields } from './modules/subset'

"""
Read a a bed file and VCF file
Generate
    - Subset VCF file
"""

params.outdir = './results'

workflow{

    check_params( [ params.bed, params.vcfs] )
    vcfs = []
    vcf = file(params.vcfs)
    if (vcf instanceof List){
        vcf.each{ vcf_ ->
            check_files([ vcf_ ])
            vcfs << [ params.name, file(params.bed), file(vcf_) ]
        }
    }
    else{
        check_files([ vcf ])
        vcfs << [ params.name, file(params.bed), file(vcf) ]
    }
        
    data = Channel.from(vcfs)

    get_bed_vcf(data)
    concat_vcf_chrms(get_bed_vcf.out.groupTuple())
    extract_fields( concat_vcf_chrms.out.map{ name, vcf -> [ name, file(vcf), params.fields ] } )

    println "Output folder: ${params.outdir}"
}