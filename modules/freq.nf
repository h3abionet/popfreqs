#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * STEP - Get population frequencies
 */
process get_maf {
    tag "get_maf_${pop}_${dataset}_${chrm}_${header}"
    // publishDir "/cbio/projects/001/clients/population_allele_freqs/freqs/${dataset}", overwrite: true, mode: 'copy'
    label "bigmem"

    input:
        tuple val(pop), val(dataset), file(pop_vcf), val(chrm), val(info_field)

    output:
        tuple val(pop), val(dataset), file(pop_maf), val(chrm)

    script:
        // header = info_fields.collect{ "\t${pop}_${it.replaceAll('INFO/','')}" }.join('') // Was used in case of multiple INFO/
        header = info_field.replaceAll('INFO/','') // Was used in case of multiple INFO/
        // params = info_fields.collect{ "\\t%${it}" }.join('') // Was used in case of multiple INFO/
        pop_maf = "${pop}_${dataset}_${chrm}_${header}.csv"
        """
        ## Compute frequency
        echo "CHROM\tPOS\tPOS\tID\t${pop}_${header}" > ${pop_maf}
        bcftools query -f '%CHROM\\t%POS\\t%POS\\t%CHROM\\_%POS\\_%REF\\_%ALT\\t%${info_field}\\n' ${pop_vcf} >> ${pop_maf}
        """
}

process generate_hdr {
    tag "hdr_${pop}_${dataset}"
    publishDir "/cbio/projects/001/clients/population_allele_freqs/freqs/${dataset}", overwrite: true, mode: 'copy'
    label "bigmem"

    input:
        tuple val(pop), val(dataset), val(annots), val(annot_infos)

    output:
        tuple val(pop), val(dataset), file("${out}.hdr")

    script:
        out = "${pop}_${dataset}_annots"
        label = pop
        template "generate_hdr.py"
}

/*
 * STEP - Combine file
 */
process combine_csv {
    tag "combine_csv_${pop}_${dataset}_${chrms_}"
    publishDir "${params.outdir}/results", overwrite: true, mode:'copy'
    label "bigmem"

    input:
        tuple val(pop), val(dataset), val(pop_mafs), val(chrms), val(ext)

    output:
        tuple val(pop), val(dataset), file(pop_maf)

    script:
        if(chrms.size() > 1){
            chrms_ = "${chrms.sort()[0]}-${chrms.sort()[-1]}"
        }
        else{
            chrms_ = chrms[0]
        }
        pop_maf = "${pop}_${dataset}_${chrms_}_MAF.${ext}"
        """
        head -n1 ${pop_mafs[0]} > ${pop_maf}
        tail -q -n+2 ${pop_mafs.join(' ')} >> ${pop_maf}
        """
}

process paste_infos {
    tag "paste_info_${pop}_${dataset}_${chrm}"
    publishDir "/cbio/projects/001/clients/population_allele_freqs/freqs/${dataset}", overwrite: true, mode: 'copy'
    label "bigmem"

    input:
        tuple val(pop), val(dataset), val(annots), val(chrm)

    output:
        tuple val(pop), val(dataset), file(out), val(chrm)

    script:
        infos = annots.join(',')
        out = "${pop}_${dataset}_${chrm}_annots.csv"
        template "paste_infos.py"
}