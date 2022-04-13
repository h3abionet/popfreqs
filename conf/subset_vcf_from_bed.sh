#!/usr/bin/env bash

cd /cbio/users/mamana/popfreq/OJ_gnomad
nextflow /users/mamana/popfreqs/subset_vcf_from_bed.nf \
-with-singularity docker://quay.io/mamanambiya/snpeff_bcftools:latest \
--vcfs='/cbio/projects/001/clients/population_allele_freqs/results/h3africa_annot_v0.0.1/h3africa_annot_v0.0.1_gnomad_3_0_b38-h3africa_v6_chr*.bcf' \
--bed='/cbio/users/mamana/popfreq/OJ_gnomad/snps_chrm_pos.csv' \
--fields='ID CHROM POS REF ALT ANN[0].EFFECT ANN[0].IMPACT ANN[0].GENE ANN[0].FEATURE ANN[0].BIOTYPE ANN[0].FEATURE gnomad_b38_AF gnomad_b38_AC gnomad_b38_AN gnomad_b38_AF_afr gnomad_b38_AC_afr gnomad_b38_AN_afr gnomad_b38_AF_amr gnomad_b38_AC_amr gnomad_b38_AN_amr gnomad_b38_AF_eas gnomad_b38_AC_eas gnomad_b38_AN_eas gnomad_b38_AF_fin gnomad_b38_AC_fin gnomad_b38_AN_fin gnomad_b38_AF_nfe gnomad_b38_AC_nfe gnomad_b38_AN_nfe' \
--name='oj_gnomad_af' \
--outdir='' \
-profile slurm \
-resume

