#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * STEP 1 - Get chip site to bed
 */
process get_chip_site_from_vcf {
    tag "get_chip_site_${chip_name}_${dataset}_${chrm}"
    publishDir "${params.outdir}/chip_site", mode: 'copy'
    label "bigmem"

    input:
        tuple val(chip_name), file(chip_file), val(dataset), file(dataset_vcf), file(dataset_sample), val(chrm) from datasets_all

    output:
        tuple val(dataset), file(dataset_vcf), file(dataset_sample), val(chip_name), file(dataset_vcf_chip), val(chrm) into get_chip_site_from_vcf

    script:
    dataset_vcf_chip = "${dataset}_${chip_name}_${chrm}.vcf.gz"
    """
    tabix -f ${dataset_vcf}
    bcftools view --regions-file ${chip_file} ${dataset_vcf} -Oz -o ${dataset}_${chip_name}_tmp.vcf.gz
    bcftools sort ${dataset}_${chip_name}_tmp.vcf.gz -T . -Oz -o ${dataset_vcf_chip}
    rm -f ${dataset}_${chip_name}_tmp*.vcf.gz
    """
}


/*
 * STEP 2 - Get in populations in sample file
 */
process get_pops {
    tag "get_pop_${dataset}"

    input:
        tuple val(dataset), file(dataset_sample)

    output:
        tuple val(dataset), file(dataset_pops)

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
        tuple val(dataset), file(vcf), file(sample), val(chrm)
    output:
        tuple val(dataset), file(out_vcf), file(sample), val(chrm)
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
        tuple val(dataset), val(pop), file(dataset_vcf), file(dataset_sample), val(chrm)

    output:
        tuple val(pop), val(dataset), file(pop_vcf), file(dataset_vcf), val(chrm)

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

process get_bed_vcf {
    tag "get_bed_${name}"
    label "bigmem"

    input:
        tuple val(name), file(bed_file), file(vcf)

    output:
        tuple val(name), file(bed_vcf)

    script:
        bed_vcf = "${name}.vcf.gz"
        """
        tabix ${vcf}
        bcftools view --regions-file ${bed_file} ${vcf} | \
        bcftools sort -T . -Oz -o ${bed_vcf}
        """
}

process concat_vcf_chrms {
   tag "concat_${name}"
   label "bigmem"
//    publishDir "${params.outdir}/${prefix}", mode: 'copy'
   
   input:
        tuple val(name), val(vcfs)
   output:
        tuple val(name), file(vcf_out)
   script:
        vcf_out = "${name}.bcf"
        if(vcfs.size() > 1){
            """
            bcftools concat ${vcfs.join(' ')} |\
            bcftools sort -T . -Ob -o ${vcf_out}
            tabix ${vcf_out}
            """
        }
        else{
            """
            bcftools sort ${vcfs.join(' ')} -T . -Ob -o ${vcf_out}
            tabix ${vcf_out}
            """
        }
}

process extract_fields1 {
   tag "extract_fields_${name}"
   label "bigmem"
   publishDir "${params.outdir}/${prefix}", mode: 'copy'
   
   input:
        tuple val(name), val(vcf)
   output:
        tuple val(name), file(csv)
   script:
        csv = "${name}.csv"
        """
        echo -e 'ID\\tCHROM\\tPOS\\tREF\\tALT\\tgnomad_b38_AF\\tgnomad_b38_AC\\tgnomad_b38_AN' > ${csv}
        bcftools query \
            -f '%ID\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/gnomad_b38_AF\\t%INFO/gnomad_b38_AC\\t%INFO/gnomad_b38_AN\\n' \
            ${vcf} >> ${csv}
        """
}

process extract_fields {
   tag "extract_fields_${name}"
   label "bigmem"
   publishDir "${params.outdir}", mode: 'copy'
   
   input:
        tuple val(name), val(vcf), val(fields)
   output:
        tuple val(name), file(csv)
   script:
        csv = "${name}.csv"
        """
        echo -e '${fields}' > ${csv}
        bcftools view ${vcf} |\
        snpsift \
            extractFields \
            - -e "." ${fields} \
            > ${csv}
        """
}

