#!/usr/bin/env nextflow
nextflow.preview.dsl=2

include { check_files } from './modules/qc'
include { get_pops; split_population; fill_tags_VCF } from './modules/subset'
include { get_maf; combine_csv as combine_freq } from './modules/freq'

workflow{
    datasets_samples = Channel.from(params.datasets).map{ dataset, vcf, sample -> [ dataset, file(sample) ] }
    datasets = Channel.from(params.datasets)
        .combine(params.chromosomes)
        .flatMap{ dataset, vcf, sample, chrm ->
            dataset_vcf = sprintf(vcf, chrm)
            check_files([ dataset_vcf, sample])
            return [[dataset, file(dataset_vcf), file(sample), chrm]]
        }
    get_pops(datasets_samples)
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
    maf_dataset = split_population.out
        .flatMap{ pop, dataset2, pop_vcf, dataset_vcf, chrm -> 
            datasets_pops = []
            datasets_pops << [ pop, dataset2, pop_vcf, chrm, 'INFO/MAF' ]
            if (!("${dataset2}:${chrm}" in dats)){
                datasets_pops << [ dataset2, dataset2, dataset_vcf, chrm, 'INFO/MAF' ]
                dats << "${dataset2}:${chrm}"
            }
            return datasets_pops
        }
    // Site only datasets
    site_datasets = Channel.from(params.site_datasets)
        .combine(params.chromosomes)
        .flatMap{ dataset3, vcf, info_tag, chrm ->
            dataset_vcf1 = sprintf(vcf, chrm)
            check_files([ dataset_vcf1 ])
            return [[dataset3, dataset3, file(dataset_vcf1), chrm, info_tag]]
        }
    site_all = maf_dataset.mix(site_datasets)
    get_maf(site_all)
    combine_freq(get_maf.out.groupTuple(by: [0,1]).map{ pop, dataset, pop_mafs, chrms -> [ pop, dataset, pop_mafs, chrms, 'frq' ] })
}