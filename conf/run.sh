#!/usr/bin/env bash

nextflow run h3abionet/popfreqs/main.nf -resume -profile slurm,singularity,test
