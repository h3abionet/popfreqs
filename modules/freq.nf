#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
 * STEP - Get population frequencies
 */
process get_maf {
    tag "get_maf_${pop}_${dataset}_${chrm}"
    label "bigmem"

    input:
        tuple pop, dataset, file(pop_vcf), chrm, info_tag

    output:
        tuple pop, dataset,  file(pop_maf), chrm

    script:
        pop_maf = "${pop}_${dataset}_${chrm}_MAF.frq"
        """
        ## Compute frequency
        echo "ID\t${pop}_MAF" > ${pop_maf}
        bcftools query -f '%CHROM\\_%POS\\_%REF\\_%ALT\\t%${info_tag}\\n' ${pop_vcf} >> ${pop_maf}
        """
}

/*
 * STEP - Combine file
 */
process combine_csv {
    tag "combine_csv_${pop}_${dataset}_${chrms_}"
    publishDir "${params.outdir}", overwrite: true, mode:'copy'
    label "bigmem"

    input:
        tuple pop, dataset, pop_mafs, chrms, ext

    output:
        tuple pop, dataset, file(pop_maf)

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