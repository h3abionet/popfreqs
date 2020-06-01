#!/usr/bin/env bash

nextflow run h3abionet/main.nf -c config/test.config -resume -profile slurm,singularity,test