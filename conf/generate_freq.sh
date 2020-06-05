#!/usr/bin/env bash

nextflow run h3abionet/popfreqs/generate_freq.nf -resume -profile slurm,singularity,test1
