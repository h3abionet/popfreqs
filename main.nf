#!/usr/bin/env nextflow
nextflow.preview.dsl=2


def check_files(file_list) {
    file_list.each { myfile ->
        if (!file(myfile).exists() && !file(myfile).isFile()) exit 1, "|-- ERROR: File ${myfile} not found. Please check your config file."
    }
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
        if(datasets.siez() > 1){
            """
            bcftools merge ${vcfs.join(' ')} |\
            bcftools sort -Oz -o ${vcf_out}
            """
        }
        else{
            """
            bcftools sort ${vcfs.join(' ')} -Oz -o ${vcf_out}
            """
        }
        
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
    params.datasets.each { name, vcf_ ->
        params.chromosomes.each { chrm ->
            vcf = sprintf(vcf_, chrm)
            check_files([vcf])
            datasets_all << [name, file(vcf), chrm]
        }
    }

    datasets_all_cha = Channel.from(datasets_all)
    preprocess(datasets_all_cha)
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