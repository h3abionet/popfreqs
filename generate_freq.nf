#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { check_files } from './modules/qc'
include { get_pops; split_population; fill_tags_VCF } from './modules/subset'
include { get_maf; combine_csv as combine_freq ; generate_hdr; paste_infos } from './modules/freq'

workflow dataset_infos{
    take: data
    main:
        info_datasets = Channel.from(params.datasets)
            .flatMap{ dataset, vcf, sample, infos ->
                data = []
                fields = []
                descs = []
                infos.each{ info_field, info_desc ->
                    fields << info_field
                    descs << info_desc
                }
                data << [ dataset, fields.join(','), descs.join(',') ]
                return data }
            .combine( data, by:[0] )
            .map{ dataset, fields, infos, pop -> [ pop, dataset, fields, infos ] }
            .groupTuple( by:[0,1] )
            .map{ pop, dataset, fields, infos -> [ pop, dataset, fields[0], infos[0] ] }
    emit:
        info_datasets
}

workflow site_infos{
    main:
        info_sites = Channel.from(params.site_datasets)
            .flatMap{ dataset, vcf, infos ->
                data = []
                fields = []
                descs = []
                infos.each{ info_field, info_desc ->
                    fields << info_field
                    descs << info_desc
                }
                data << [ dataset, dataset, fields.join(','), descs.join(',') ]
                return data
            }
    emit:
        info_sites
}

workflow{
    
    //// Datasets with genotypes 
    datasets_samples = Channel.from(params.datasets)
        .map{ dataset, vcf, sample, infos -> [ dataset, file(sample) ] }
    get_pops(datasets_samples)
    datasets = Channel.from(params.datasets)
        .combine(params.chromosomes)
        .flatMap{ dataset, vcf, sample, infos, chrm ->
            dataset_vcf = sprintf(vcf, chrm)
            check_files([ dataset_vcf, sample])
            return [[dataset, file(dataset_vcf), file(sample), chrm]]
        }
    fill_tags_VCF(datasets)
    pop_data = get_pops.out
        .flatMap{ dataset1, pops_file ->
            datasets_pops = []
            pop_data = file(pops_file).readLines().unique().sort()
            pop_data.each { pop ->
                if(!(pop in ['POP'])) {
                    if ('ALL' in params.POPS) {
                        datasets_pops << [dataset1, pop]
                    }
                    else if (pop in params.POPS) {
                        datasets_pops << [dataset1, pop]
                    }
                }
            }
            return datasets_pops 
        }
        .combine(fill_tags_VCF.out, by:0)
    split_population(pop_data)
    dats = []

    //// Dataset AFs
    maf_dataset = fill_tags_VCF.out
        .combine( Channel.from(params.datasets), by:[0] )
        .flatMap{ dataset, dataset_vcf, dataset_sample, chrm, dataset_vcf_orginal, dataset_sample_original, infos ->
            data = []
            infos.each{ info_field, info_desc ->
                data << [ dataset, dataset, file(dataset_vcf), chrm, info_field ]
            }
            return data
        }
    
    //// Population AFs
    maf_dataset_pop = split_population.out
        .map{ pop, dataset, pop_vcf, dataset_vcf, chrm -> [ dataset, pop, pop_vcf, dataset_vcf, chrm ] }
        .combine( Channel.from(params.datasets), by:[0] )
        .flatMap{ dataset, pop, pop_vcf, dataset_vcf, chrm, dataset_vcf_orginal, dataset_sample, infos ->
            data = []
            infos.each{ info_field, info_desc ->
                data << [ pop, dataset, file(pop_vcf), chrm, info_field ]
            }
            return data
        }
    //// Datasets with sites only 
    site_datasets = Channel.from(params.site_datasets)
        .combine(params.chromosomes)
        .flatMap{ dataset3, vcf, infos, chrm ->
            dataset_vcf1 = sprintf(vcf, chrm)
            check_files([ dataset_vcf1 ])
            data = []
            infos.each{ info_field, info_desc ->
                data << [ dataset3, dataset3, file(dataset_vcf1), chrm, info_field ]
            }
            return data
        }
    site_all = maf_dataset.mix(maf_dataset_pop, site_datasets)
    get_maf(site_all)
    paste_infos(get_maf.out.groupTuple(by:[0,1,3]))

    //// Process info headers
    info_datasets = split_population.out
        .map{ pop, dataset, pop_vcf, dataset_vcf, chrm -> [ dataset, pop ] }
        .mix( Channel.from(params.datasets).map{ dataset, vcf, sample, infos -> [ dataset, dataset ] } )
    dataset_infos( info_datasets )
    site_infos()
    info_all = dataset_infos.out.info_datasets.mix(site_infos.out.info_sites)
    generate_hdr(info_all)

    // combine_freq(get_maf.out.groupTuple(by: [0,1]).map{ pop, dataset, pop_mafs, chrms -> [ pop, dataset, pop_mafs, chrms, 'frq' ] }).view()
    
    println "Output: ${params.outdir} ...."
}