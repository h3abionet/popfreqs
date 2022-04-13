#!/usr/bin/env bash

#SBATCH --partition=Main
#SBATCH --nodes=1 --ntasks=5 --mem=7000
#SBATCH --time=48:00:00
#SBATCH --job-name="popfreq"
#SBATCH --mail-user=mbymam001@myuct.ac.za
#SBATCH --mail-type=BEGIN,END,FAIL

OUTDIR="./annotate_vcf"
mkdir -p ${OUTDIR}
cd ${OUTDIR}

nextflow -log nextflow.log \
run h3abionet/popfreqs/main.nf \
-c conf/test.1.config \
-profile slurm,singularity \
-resume
