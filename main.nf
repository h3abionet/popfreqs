#!/usr/bin/env nextflow
nextflow.preview.dsl=2


def check_files(file_list) {
    file_list.each { myfile ->
        if (!file(myfile).exists() && !file(myfile).isFile()) exit 1, "|-- ERROR: File ${myfile} not found. Please check your config file."
    }
}

process get_pops {
    tag "get_pop_${dataset}"

    input:
        tuple dataset, file(sample)

    output:
        tuple dataset, file(pops)

    script:
        pops = "${sample.baseName}_pops.csv"
        """
        cat ${sample} | grep -v "^#" | awk -F' ' '{print \$2}' | sort -n | uniq > ${pops}
        """
}

process sites_only {
    tag "sites_only_${dataset}_${chrm}"
    label "bigmem"

    input:
        tuple dataset, file(vcf), chrm

    output:
        tuple dataset, file(sites_vcf), chrm

    script:
        sites_vcf = "${dataset}_${chrm}_sites.vcf.gz"
        """
        tabix ${vcf}
        bcftools view ${vcf} --drop-genotypes -Oz -o ${sites_vcf}
        """
}

process concat_dataset {
   tag "concat_dataset_${dataset}"
   label "bigmem"
   
   input:
       tuple dataset, vcfs
   output:
       tuple dataset, file(vcf_out)
   script:
       vcf_out = "${dataset}.vcf.gz"
       """
       bcftools concat ${vcfs.join(' ')} |\
       bcftools sort -Oz -o ${vcf_out}
       tabix ${vcf_out}
       """
}

process merge_groups {
    tag "merge_groups_${group}"
    label "bigmem"
    publishDir "${params.outdir}/${group}", mode: 'copy'

    input:
        tuple dataset, datasets, vcfs
    output:
        tuple group, file(vcf_out)
    script:
        group = datasets.join('-')
        vcf_out = "${group}.vcf.gz"
        """
        bcftools merge ${vcfs.join(' ')} |\
        bcftools sort  -Oz -o ${vcf_out}
        """
}

process get_vcf_site {
    tag "get_vcf_site_${dataset}"
    label "bigmem"

    input:
        tuple dataset, file(vcf)
    
    output:
        tuple dataset, file(vcf), file(vcf_sites)
    
    script:
        base = file(vcf.baseName).baseName
        vcf_sites = "${base}.sites"
        """
        echo -e 'rsID' > ${vcf_sites}
        bcftools query -f "%CHROM\\_%POS\\_%REF\\_%ALT\\n" ${vcf} >> ${vcf_sites}
        """
}

process process_maf_dataset {
    tag "process_maf_${dataset_name}_${mafs_dataset}"
    label "bigmem"
    
    input:
        tuple dataset_name, file(dataset_vcf), file(sites), mafs_dataset, mafs_files, annot
    
    output:
        tuple dataset_name, file(dataset_vcf), mafs_dataset, file(outAnnot), file(outHdr)
    
    script:
        inTSV = mafs_files
        outAnnot = "${file(dataset_vcf.baseName).baseName}_mafs-${mafs_dataset}.tsv"
        outHdr = "${file(dataset_vcf.baseName).baseName}_mafs-${mafs_dataset}.hdr"
        template "freq_annotation.py"
}

process annotate_mafs {
    tag "mafs_${dataset_name}_${mafs_dataset}"
    label "bigmem"
    
    input:
        tuple dataset_name, file(dataset_vcf), mafs_dataset, file(maf_annot), file(header_annot)
    
    output:
        tuple dataset_name, file(outVCF)
    
    script:
        base = "${dataset_name}_${maf_annot.baseName}"
        outVCF = "${base}.vcf.gz"
        """
        columns=\$(awk 'NR==1{print \$0}' ${maf_annot} | sed -e 's/\\s\\+/,/g')
        tail -n+2 ${maf_annot} | sort -k1,1n -k2,2n -V > ${maf_annot}.annot.tsv
        bgzip ${maf_annot}.annot.tsv
        tabix -s1 -b2 -e3 ${maf_annot}.annot.tsv.gz
        bcftools sort ${dataset_vcf} -T . -Oz -o ${dataset_vcf}.sorted.vcf.gz
        tabix ${dataset_vcf}.sorted.vcf.gz
        bcftools annotate -a ${maf_annot}.annot.tsv.gz -h ${header_annot} -c \${columns} -Oz -o ${outVCF} ${dataset_vcf}.sorted.vcf.gz
        tabix ${outVCF}
        rm ${dataset_vcf}.sorted.vcf.gz
        """
}

workflow split_data_pop {
    take:
        data
    main:
        get_pops(data)

    emit:
        pops = get_pops.out
}

workflow preprocess {
    take:
        data
    main:
        sites_only(data)
        concat_dataset(sites_only.out.map{ it -> [ it[0], it[1] ]}.groupTuple(by:0))
        merge_groups(concat_dataset.out.map{ dataset, vcf -> [ 'GROUP', dataset, file(vcf) ] }.groupTuple(by:0))
        get_vcf_site(merge_groups.out)
    emit:
        vcf_sites = get_vcf_site.out
}

workflow annotate {
    take:
        data
    main:
        process_maf_dataset(data)
        annotate_mafs(process_maf_dataset.out)
    emit:
        annotated_vcfs = annotate_mafs.out
}

workflow postprocess {
    take:
        data
    main:
        merge_groups(data)
    emit:
        merged_vcfs = merge_groups.out
}


workflow{
    datasets_all = []
    datasets_samples = []
    params.datasets.each { name, vcf_, sample ->
        check_files([sample])
        datasets_samples << [name, file(sample)]
        params.chromosomes.each { chrm ->
            vcf = sprintf(vcf_, chrm)
            check_files([vcf])
            datasets_all << [name, file(vcf), file(sample), chrm]
        }
    }
    samples = Channel.from(datasets_samples)
    split_data_pop(samples)

    split_data_pop.out.pops.flatMap{ dataset, pops ->
        data_pop = []
        dataset_pops = file(pops).readLines().unique().sort()
        dataset_pops.each { pop ->
            if ('ALL' in params.POPS) {
                data_pop << [dataset, pop]
            }
            else if (pop in params.POPS) {
                data_pop << [dataset, pop]
            }
        }
        return data_pop
    }

    datasets_all_cha = Channel.from(datasets_all).map{ it -> [ it[0], it[1], it[3] ] }
    preprocess(datasets_all_cha)
    // Channel.fromPath("/scratch/users/mamana/popfreqs/**/*_1-22_MAF.frq").view()
    mafs = Channel
        .from(file(params.mafs_annotations))
        .splitCsv(strip: true, sep: ',')
        .map{ it ->
            if ( it[0][0] != '#'){
                return it
            }
        }
    annotate(preprocess.out.vcf_sites.combine(mafs)).view()
    postprocess( annotate.out.annotated_vcfs.groupTuple().map{ it -> [ it[0], [it[0]], it[1] ] }.view() )
}


// /*
//  * STEP 2 - Get in populations in sample file
//  */



// datasets_pops_all = Channel.from(datasets_all).combine(datasets_pops, by:[0])
// /*
//  * STEP 2 - Split in populations
//  */
// process split_population {
//     tag "split_population_${dataset}_${pop}_${chrm}"
//     label "bigmem"

//     input:
//     set dataset, file(dataset_vcf), file(dataset_sample), chrm, pop from datasets_pops_all

//     output:
//     set dataset, file(dataset_vcf), file(dataset_sample), chrm, pop, file(pop_vcf) into split_population

//     script:
//     pop_vcf = "${pop}_${dataset}_${chrm}.vcf.gz"
//     """
//     grep ${pop} ${dataset_sample} | cut -f1 > ${pop}.samples
//     ## Keep only samples for population and Recalculate AC, AN, AF
//     bcftools view \
//         --samples-file ${pop}.samples \
//         ${dataset_vcf} | \
//     bcftools +fill-tags -- -t MAF| \
//     bcftools annotate \
//         --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' | \
//         bgzip -c > ${pop_vcf}
//     """
// }

// /*
//  * STEP - Get population frequencies
//  */
// process pop_freq {
//     tag "pop_freq_${dataset}_${pop}_${chrm}"
//     publishDir "${params.outdir}/${pop}", mode: 'copy'
//     label "bigmem"

//     input:
//     set dataset, file(dataset_vcf), file(dataset_sample), chrm, pop, file(pop_vcf) from split_population

//     output:
//     set pop, chrm, file(pop_maf) into pop_freq

//     script:
//     pop_maf = "${pop}_${dataset}_${chrm}_MAF.frq"
//     """
//     ## Compute frequency
//     echo "rsID\tREF\tALT\t${pop}_MAF_FREQ" > ${pop_maf}
//     bcftools query \
//         -f '%ID\\t%REF\\t%ALT\\t%INFO/MAF\\n' \
//         ${pop_vcf} >> ${pop_maf}
//     """
// }

// /*
//  * STEP - Get population frequencies
//  */
// process combine_pop_freq {
//     tag "combine_pop_freq_${pop}"
//     publishDir "${params.outdir}/${pop}", mode: 'copy'
//     label "bigmem"

//     input:
//     set pop, chrm, pop_mafs from pop_freq.groupTuple(by: 0)

//     output:
//     set pop, file(pop_maf) into combine_pop_freq

//     script:
//     chrm = chrm.sort()
//     if (chrm.size() > 1){
//         chrms = "${chrm[0]}-${chrm[-1]}"
//     }
//     else{
//         chrms = chrm[0]
//     }
//     pop_maf = "${pop}_${chrms}_MAF.frq"
//     """
//     head -n1 ${pop_mafs[0]} > ${pop_maf}
//     tail -q -n +2 ${pop_mafs.join(' ')} >> ${pop_maf}
//     """
// }

// combine_pop_freq_data = combine_pop_freq.toSortedList( { a, b -> a[0] <=> b[0] } ).val

// report_data = []
// combine_pop_freq_data.each { data ->
//     report_data << data[1]
// }


// // process report_freq {
// //     tag "report_freq"
// //     publishDir "${params.outdir}/report_frqs", mode: 'copy'
// //     label "bigmem"

// //     input:
// //     val(pop_freqs) from report_data.join(' ')

// //     output:
// //     file(report_maf) into report_freq

// //     script:
// //     report_maf = "H3Africa_pops_MAF.frq"
// //     template "report_freq.py"
// // }



