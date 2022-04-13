#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// check if files exist [name, file1, file2, ...]
def check_files(file_list) {
    file_list.each { myfile ->
        if (!file(myfile).exists() && !file(myfile).isFile()) exit 1, "|-- ERROR: File ${myfile} not found. Please check your config file."
    }
}

def check_params(params_list) {
    params_list.each { param ->
        if (param == '') exit 1, "|-- ERROR: Empty params --bed. Please check your parameters."
    }
}