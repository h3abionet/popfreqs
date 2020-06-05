#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// check if files exist [name, file1, file2, ...]
def check_files(file_list) {
    file_list.each { myfile ->
        if (!file(myfile).exists() && !file(myfile).isFile()) exit 1, "|-- ERROR: File ${myfile} not found. Please check your config file."
    }
}