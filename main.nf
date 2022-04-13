#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// check if files exist [name, file1, file2, ...]
def check_files(file_list) {
    file_list.each { myfile ->
        if (!file(myfile).exists() && !file(myfile).isFile()) exit 1, "|-- ERROR: File ${myfile} not found. Please check your config file."
    }
}

process sites_only {
    tag "sites_only_${dataset}_${chrm}:${start}-${end}"
    label "bigmem"

    input:
        tuple val(dataset), val(chrm), val(start), val(end), val(chip), file(vcf)

    output:
        tuple val(dataset), val(chrm), val(start), val(end), val(chip), file(sites_vcf)

    script:
        sites_vcf = "${vcf.getSimpleName()}_sites.bcf"
        """
        tabix ${vcf}
        bcftools view ${vcf} --drop-genotypes --threads ${task.cpus} -Ob -o ${sites_vcf}
        tabix ${sites_vcf}
        """
}

process merge_datasets {
   tag "merge_datasets_${dataset}_${chrm}_${start}_${end}"
   label "bigmem"
   
   input:
        tuple val(datasets), val(chrm), val(start), val(end), val(chips), val(vcfs)
   output:
        tuple val(dataset), val(chrm), val(start), val(end), val(chip), file(vcf_out)
   script:
        dataset = datasets.sort().join('-')
        chip = chips.sort().join('-')
        vcf_out = "${dataset}_${chrm}_${start}_${end}.bcf"
        if(vcfs.size() > 1){
            """
            bcftools merge ${vcfs.join(' ')} |\
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

process concat_chrms {
   tag "concat_${dataset}_${chrm}_${prefix}"
   label "bigmem"
   publishDir "${params.outdir}/${prefix}", mode: 'copy'
   //TODO: copy config file
   // Create REDME
   
   input:
        tuple val(dataset), val(chrm), val(vcfs), val(prefix)
   output:
        tuple val(dataset), val(chrm), file(vcf_out)
   script:
        vcf_out = "${prefix}_${dataset}_${chrm}.bcf"
        //println params
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

process merge_groups {
    tag "merge_groups_${dataset}_${prefix}"
    label "largemem"
    // publishDir "${params.outdir}/${dataset}", mode: 'copy'

    input:
        tuple val(dataset), val(chrms), val(vcfs), val(prefix)
    output:
        tuple val(dataset), file(vcf_out)
    script:
        vcf_out = "${prefix}.vcf.gz"
        if(vcfs.size() > 1){
            """
            bcftools merge ${vcfs.join(' ')} |\
            bcftools sort -T . -Ob -o ${vcf_out}
            """
        }
        else{
            """
            bcftools sort ${vcfs.join(' ')} -T . -Ob -o ${vcf_out}
            """
        }
        
}

process get_map {
    tag "get_map_${dataset}"
    label "bigmem"

    input:
        tuple val(dataset), file(dataset_vcf)
    output:
        tuple val(dataset), file(dataset_vcf), file(dataset_map)
    script:
        base = file(dataset_vcf.baseName).baseName
        dataset_map = "${base}.map"
        """
        bcftools query -f '%CHROM\\t%POS\\n' ${dataset_vcf} > ${dataset_map}
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
        echo -e 'ID' > ${vcf_sites}
        bcftools query -f "%CHROM\\_%POS\\_%REF\\_%ALT\\n" ${vcf} >> ${vcf_sites}
        """
}

process get_vcf_site1 {
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
        echo -e 'ID' > ${vcf_sites}
        bcftools query -f "%CHROM\\_%POS\\_%REF\\_%ALT\\n" ${vcf} >> ${vcf_sites}
        """
}

process get_vcf_site2 {
    tag "get_vcf_site_${dataset}"
    publishDir "${params.outdir}/sites", mode: 'copy', pattern: "*.sites"
    label "bigmem"

    input:
        tuple val(dataset), file(vcf)
    
    output:
        tuple val(dataset), file(vcf), file(vcf_sites)
    
    script:
        base = file(vcf.baseName).baseName
        vcf_sites = "${base}_with_gene.sites"
        """
        bcftools view ${vcf} |\
        SnpSift -Xmx${task.memory.toGiga()}g \
            extractFields \
            - -e "." CHROM POS ANN[0].GENE ANN[0].EFFECT \
            > ${vcf_sites}
        """
}

process process_annotation {
    tag "process_maf_${name}"
    label "large"
    
    input:
        tuple val(name), val(descr), path(annot_file), path(hdr_file)
    
    output:
        tuple val(name), path("${prefix_out}*.tsv"), path(hdr_file), path(annot_col)
    
    script:
        base = file(annot_file.baseName).baseName
        prefix_out = "${base}"
        annot_hdr_file = "${base}.hdr"
        annot_col = "${base}.columns"
        annotation = descr
        template "process_annotation.py"
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

process process_annot_chrm {
    tag "mafs_${annot}_${chrm}"
    label "bcftools"
    
    input:
        tuple val(annot), val(chrm), file(annot_bed), file(annot_hdr), val(columns)
    
    output:
        tuple val(annot), val(chrm), file(annot_bed_tabix), file(annot_hdr), val(columns)
    
    script:
        annot_bed_tabix = "${annot_bed}.gz"
        """
        tail -n+2 ${annot_bed} | sort -k1,1n -k2,2n -V | bgzip > ${annot_bed_tabix}
        tabix -s1 -b2 -e3 ${annot_bed_tabix}
        """
}

process annotate_annot_chunk {
    tag "annot_${dataset_name}_${annot}"
    label "bcftools"
    
    input:
        tuple val(dataset_name), val(chrm), val(start), val(end), val(chip), file(dataset_vcf), val(annot), file(annot_bed), file(annot_hdr), file(annot_col)
    
    output:
        tuple val(dataset_name), val(chrm), val(start), val(end), file(outVCF)
    
    script:
        base = "${dataset_name}_${annot}_${chrm}_${start}_${end}"
        outVCF = "${base}.vcf.gz"
        """
        columns=\$(cat ${annot_col})
        tabix -s1 -b2 -e3 ${annot_bed}
        bcftools annotate -x INFO,^FORMAT/GT,FORMAT/PL ${dataset_vcf} | \
        bcftools sort  -T . | \
        bcftools annotate -a ${annot_bed} -h ${annot_hdr} -c \${columns} | \
        bcftools annotate -x ID -Oz -o ${outVCF}
        tabix ${outVCF}
        """
}


'''
Step: Annotate dataset using snpEff database
'''
process snpeff_vcf {
    tag "snpeff_${dataset}_${chrm}"
    label "snpeff_bcftools"
    input:
        tuple val(dataset), val(chrm), file(vcf_file)
    output:
        tuple val(dataset), val(chrm), file(vcf_out)
    script:
        base = vcf_file.getSimpleName()
        vcf_out = "${base}_snpeff.vcf.gz"
        """
        snpEff \
            -Xmx${task.memory.toGiga()}g \
            ${params.snpEff_human_db} \
            -lof \
            -stats ${base}_snpeff.html \
            -csvStats ${base}_snpeff.csv \
            -dataDir ${params.snpEff_database} \
            ${vcf_file} |\
        bcftools view -Oz -o ${vcf_out}
        tabix ${vcf_out}
        #rm ${base}_snpeff.vcf
        """
}

process update_rsid_vcf {
    tag "update_rsid_vcf_${dataset}_${chrm}"
    label "bigmem"
    input:
        tuple val(dataset), val(chrm), file(vcf_file), file(dbsnp_vcf)
    output:
        tuple val(dataset), val(chrm), file(vcf_out)
    script:
        base = vcf_file.getSimpleName()
        vcf_out = "${base}_dbsnp.vcf.gz"
        """
        tabix -f ${vcf_file}
        tabix -f ${dbsnp_vcf}
        bcftools annotate ${vcf_file} -a ${dbsnp_vcf} -c ID --threads ${task.cpus} -Oz -o ${vcf_out}
        tabix ${vcf_out}
        """
}

process annot_vcf_chrm {
    tag "annot_vcf_chrm_${dataset}_${chrm}"
    label "medium"
    publishDir "/cbio/users/mamana/exome_aibst/data/AIBST/VCF/CHRS/", mode: 'copy'

    input:
        tuple val(dataset), val(chrm), file(annot_vcf), file(vcf_file)
    output:
        tuple val(dataset), val(chrm), file(vcf_chrm)
    script:
        base = file(vcf_file.baseName).baseName
        vcf_chrm = "${base}.annot.vcf.gz"
        """
        tabix ${vcf_file}
        tabix ${annot_vcf}
        bcftools view --regions ${chrm} ${vcf_file} --threads ${task.cpus} -Oz -o ${base}_${chrm}.annot.vcf.gz
        bcftools annotate -a ${annot_vcf} -c INFO ${vcf_file} --threads ${task.cpus} -Oz -o ${vcf_chrm}
        tabix ${vcf_chrm}
        """
}

process generate_chunks_vcf {
    tag "generate_chunks_${dataset}"
    label "bigmem"

    input:
        tuple val(dataset), file(vcf), file(mapFile), val(chrms), val(chunk_size)
    output:
        tuple val(dataset), file(vcf), file(chunkFile)
    script:
        chromosomes = chrms
        chunk = ''
        chunkFile = "chunks.txt"
        template "generate_chunks.py"
}

process split_target_to_chunk_sites {
    tag "split_${dataset}_${chrm}:${chunk_start}-${chunk_end}"
    label "bigmem"

    input:
        tuple val(dataset), val(chrm), val(chunk_start), val(chunk_end), val(chip), file(dataset_vcf)
    output:
        tuple val(dataset), val(chrm), val(chunk_start), val(chunk_end), val(chip), file(vcf_chunk_out)
    script:
        base = file(dataset_vcf.baseName).baseName
        vcf_chunk_out = "${base}_${chrm}_${chunk_start}-${chunk_end}_${chip}.sites.bcf"
        """
        tabix ${dataset_vcf}
        bcftools view --regions ${chrm}:${chunk_start}-${chunk_end} ${dataset_vcf} --drop-genotypes --threads ${task.cpus} -Ob -o ${vcf_chunk_out}
        tabix ${vcf_chunk_out}
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
        get_map(data)
        generate_chunks_vcf(get_map.out.map{ dataset, vcf, map_file -> [ dataset, file(vcf), file(map_file), '', params.chunk_size ] })
        chunks_datas = generate_chunks_vcf.out.flatMap{ dataset, vcf, chunk_file ->
            datas = []
            chunks = file(chunk_file).text.split()
            chunks.each{ chunk_data ->
                data = chunk_data.split(',')
                chrm = data[0]
                chunk_start = data[1]
                chunk_end = data[2]
                datas << [dataset, chrm, chunk_start, chunk_end, dataset, file(vcf)]
            }
            return datas
        }
        split_target_to_chunk_sites(chunks_datas)
        // sites_only(split_target_to_chunk.out)
        merge_ch = split_target_to_chunk_sites.out.map{ dataset, chrm, start, end, chip, sites_vcf -> [ dataset, chrm, start, end, chip, file(sites_vcf) ]}.groupTuple(by:[1,2,3])
        merge_datasets(merge_ch)

    emit:
        vcf_sites = merge_datasets.out
}

workflow annotate {
    take:
        data
    main:
        //// Preprocess freq files, split into chromosomes
        check_files([ params.mafs_annotations, params.dbsnp_vcf ])
        freqs = Channel
            .from(file(params.mafs_annotations))
            .splitCsv(strip: true, sep: ',')
            .flatMap{ name, descr, annots, hdr ->
                annot_files = file(annots)
                datas = []
                // TODO check chromosomes in datasets
                if ( name[0] != '#'){
                    if (annot_files instanceof List){
                        annot_files.each{ annot ->
                            check_files([ annot, hdr ])
                            datas << [ name, descr, file(annot), file(hdr) ]
                        }
                    }
                    else{
                        check_files([ annots, hdr ])
                        datas << [ name, descr, file(annots), file(hdr) ]
                    }
                }
                return datas
            }
        process_annotation(freqs)
        process_annotation_chrm = process_annotation.out
            .map{ annot_name, annot_file, annot_hdr, annot_col ->
            // annot_files.each{ annot -> // Not needed now as only have one chrm 
                chrm = annot_file =~ /__(.+)\.tsv/
                datas = [ annot_name, chrm[0][1], file(annot_file), file(annot_hdr), annot_col ]
            // }
            return datas
        }
        process_annot_chrm(process_annotation_chrm)
        annotate_pop_cha =  data.combine( process_annot_chrm.out, by:1 )
            .map{ chrm, dataset, start, end, chip, dataset_vcf, annot, annot_bed, annot_hdr, annot_col ->[ dataset, chrm, start, end, chip, file(dataset_vcf), annot, file(annot_bed), file(annot_hdr), file(annot_col) ] }
        annotate_annot_chunk(annotate_pop_cha)
        if(params.update_rsid == true){
            // annotate_annot_chunk.out.flatMap{ dataset, chrm, start, end, vcf ->
            //     println file(vcf)
            //     println file(vcf).countLines()
            // }
            update_rsid_vcf( annotate_annot_chunk.out.map{ dataset, chrm, start, end, vcf -> [ dataset, chrm, file(vcf), file(params.dbsnp_vcf) ] } )
            snpeff_vcf( update_rsid_vcf.out )
            annotated_vcfs = snpeff_vcf.out
        }
        else{
            annotated_vcfs = annotate_annot_chunk.out.map{ dataset, chrm, start, end, vcf -> [ dataset, chrm, file(vcf) ] }
        }

    emit:
        data
        // annotated_vcfs = update_rsid_vcf.out
        annotated_vcfs
}

workflow annotate_pop {
    take:
        data
    main:
        snpeff_vcf(data.map{ dataset, vcf, sites -> [ dataset, '', file(vcf)] })
        get_vcf_site2(snpeff_vcf.out.map{ dataset, chrm, vcf -> [ dataset, file(vcf) ] })
    emit:
        annotated_vcfs = data
}

workflow postprocess {
    take:
        data
    main:
        concat_chrms( data.map{ dataset, chrm, vcfs -> [ dataset, chrm, vcfs, params.out_prefix ]} )
    emit:
        concat_vcfs = concat_chrms.out
}

workflow{

    datasets = []
    params.datasets.each { dataset, dataset_vcf ->
        datas = []
        dataset_vcfs = file(dataset_vcf)
        if (dataset_vcfs instanceof List){
            dataset_vcfs.each{ vcf ->
                datas = [ dataset ]
                check_files([ vcf ])
                datas << file(vcf)
                datasets << datas 
            }
        }
        else{
            datas = [ dataset ]
            check_files([ dataset_vcf ])
            datas << file(dataset_vcf)
            datasets << datas 
        }
        
    }
    datasets_ch = Channel.from(datasets)
    preprocess(datasets_ch)
    annotate( preprocess.out.vcf_sites )
    annotated_data = annotate.out.annotated_vcfs.groupTuple(by:[0,1])
    postprocess( annotated_data )

    println "Output folder: ${params.outdir}"
}