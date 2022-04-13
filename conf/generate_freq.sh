#!/usr/bin/env bash

#SBATCH --partition=Main
#SBATCH --nodes=1 --ntasks=5 --mem=7000
#SBATCH --time=48:00:00
#SBATCH --job-name="popfreq"
#SBATCH --mail-user=mbymam001@myuct.ac.za
#SBATCH --mail-type=BEGIN,END,FAIL

OUTDIR="/scratch3/users/mamana/popfreq/generate_af"
mkdir -p ${OUTDIR}
cd ${OUTDIR}

nextflow run /users/mamana/popfreqs/generate_freq.nf \
    -c /users/mamana/popfreqs/conf/test1.config \
    -profile slurm,singularity \
    -resume 

# nextflow run h3abionet/popfreqs/generate_freq.nf -resume -profile slurm,singularity,test1
