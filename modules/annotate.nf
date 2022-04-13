#!/usr/bin/env nextflow
nextflow.enable.dsl=2


'''
Step: Annotate dataset using snpEff database
'''
process snpeff_vcf {
    tag "snpeff_${dataset}_${chrm}"
    label "bigmem"
    input:
        tuple val(dataset), val(chrm), file(vcf_file)
    output:
        tuple val(dataset), val(chrm), file("${vcf_out}.gz")
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_snpeff.vcf"
        """
        bcftools index --tbi -f ${vcf_file}
        snpEff \
            ${params.snpEff_human_db} \
            -lof \
            -stats ${vcf_file.baseName}.html \
            -csvStats ${vcf_file.baseName}.csv\
            ${vcf_file} > ${vcf_out} -v
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}

// -Xmx${task.memory.toGiga()}g \

process generate_vcf_from_bed {
    tag "vcf_from_bed_${manifest_name}"
    publishDir "${params.outdir}/results", overwrite: true, mode:'copy'
    label "bigmem"
    input:
        tuple val(manifest_name), file(manifest_bed), file(human_genome_ref)
    output:
        tuple val(manifest_name), file(vcf_out), file(unmapped)
    script:
        vcf_out = "${manifest_bed.baseName}.vcf.gz"
        unmapped = "${manifest_bed.baseName}.unmapped.txt"
        """
        bcftools convert -c ID,CHROM,POS,AA -f ${human_genome_ref} --tsv2vcf ${manifest_bed} --samples -- |\
        bcftools sort -T . |\
        bcftools view --drop-genotypes -Oz -o ${vcf_out}
        cp .command.err ${unmapped}
        #bcftools index --tbi -f ${vcf_out}
        """
}

process parse_manifests {
    tag "parse_manifests_${manifest_name}"
    publishDir "${params.outdir}/results", overwrite: true, mode:'copy'
    label "bigmem"
    
    input:
        tuple val(manifest_name), file(manifest_file), val(coordinate_fields), val(annots), val(annot_infos), val(annot_types), val(sep)
    
    output:
        tuple val(manifest_name), file("${manifest_out}*annots.csv"), file("${manifest_out}.hdr.csv"), file("${manifest_out}.bed")
    
    script:
        manifest_out = "${manifest_file.baseName}"
        template "parse_manifests.py"
}