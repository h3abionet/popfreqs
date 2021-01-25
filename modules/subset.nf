#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
 * STEP 1 - Get chip site to bed
 */
process get_chip_site_from_vcf {
    tag "get_chip_site_${chip_name}_${dataset}_${chrm}"
    publishDir "${params.outdir}/chip_site", mode: 'copy'
    label "bigmem"

    input:
    set chip_name, file(chip_file), dataset, file(dataset_vcf), file(dataset_sample), chrm from datasets_all

    output:
    set dataset, file(dataset_vcf), file(dataset_sample), chip_name, file(dataset_vcf_chip), chrm into get_chip_site_from_vcf

    script:
    dataset_vcf_chip = "${dataset}_${chip_name}_${chrm}.vcf.gz"
    """
    tabix -f ${dataset_vcf}
    bcftools view --regions-file ${chip_file} ${dataset_vcf} -Oz -o ${dataset}_${chip_name}_tmp.vcf.gz
    bcftools sort ${dataset}_${chip_name}_tmp.vcf.gz -Oz -o ${dataset_vcf_chip}
    rm -f ${dataset}_${chip_name}_tmp*.vcf.gz
    """
}


/*
 * STEP 2 - Get in populations in sample file
 */
process get_pops {
    tag "get_pop_${dataset}"

    input:
        tuple dataset, file(dataset_sample)

    output:
        tuple dataset, file(dataset_pops)

    script:
        dataset_pops = "${dataset_sample.baseName}_pops.csv"
        """
        tail -q -n+2 ${dataset_sample} | grep -v "^#" | awk -F' ' '{print \$2}' | sort -n | uniq > ${dataset_pops}
        """
}

process fill_tags_VCF {
    tag "fill_tags_${dataset}_${chrm}"
    label "hugemem"

    input:
        tuple dataset, file(vcf), file(sample), chrm
    output:
        tuple dataset, file(out_vcf), file(sample), chrm
    script:
        base = file(vcf.baseName).baseName
        out_vcf = "${base}_AF.vcf.gz"
        """
        tabix -f ${vcf}
        bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' ${vcf} | \
        bcftools +fill-tags -Oz -o ${out_vcf}
        """
}

/*
 * STEP 2 - Split in populations
 */
process split_population {
    tag "split_pop_${dataset}_${pop}_${chrm}"
    label "bigmem"
    // errorStrategy 'ignore'

    input:
        tuple dataset, pop, file(dataset_vcf), file(dataset_sample), chrm

    output:
        tuple pop, dataset, file(pop_vcf), file(dataset_vcf), chrm

    script:
        pop_vcf = "${pop}_${dataset}_${chrm}.vcf.gz"
        """
        awk '\$2=="${pop}" {print \$1}' ${dataset_sample} > ${pop}.samples
        ## Keep only samples for population and Recalculate AC, AN, AF
        bcftools view --samples-file ${pop}.samples --force-samples ${dataset_vcf} | \
        bcftools +fill-tags -- -t AC,AN,AF,MAF | \
        bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' | \
        bcftools view --drop-genotypes -Oz -o ${pop_vcf}
        """
}
