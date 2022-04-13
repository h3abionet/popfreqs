#!/usr/bin/env bash

#SBATCH --partition=Main
#SBATCH --nodes=1 --ntasks=5 --mem=7000
#SBATCH --time=48:00:00
#SBATCH --job-name="popfreq"
#SBATCH --mail-user=mbymam001@myuct.ac.za
#SBATCH --mail-type=BEGIN,END,FAIL

OUTDIR="/scratch3/users/mamana/popfreq_test"
mkdir -p ${OUTDIR}
cd ${OUTDIR}

nextflow run /users/mamana/popfreqs/parse_manifests.nf \
    -c /users/mamana/popfreqs/conf/parse_manifests.config \
    -profile slurm,singularity \
    -resume 
