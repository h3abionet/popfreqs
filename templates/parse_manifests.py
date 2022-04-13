#!/usr/bin/env python3.7
'''

'''
import argparse,sys,time
import pandas as pd
import numpy as np
import csv

parser = argparse.ArgumentParser()
parser.add_argument("--manifest_file",
                    default="${manifest_file}", help="One TSV annotation files")
parser.add_argument(
    "--manifest_name", default="${manifest_name}", help="TSV annotation file")
parser.add_argument("--coordinate_fields",
                    default="${coordinate_fields}", help="")
parser.add_argument("--sep", default="${sep}", help="")
parser.add_argument("--annots", default="${annots}", help="")
parser.add_argument("--annot_infos", default="${annot_infos}", help="")
parser.add_argument("--annot_types", default="${annot_types}", help="")
parser.add_argument("--manifest_out",
                    default="${manifest_out}", help="")
args = parser.parse_args()


def parse_manifest(manifest_name, manifest_file, coordinate_fields, manifest_out, annots, annot_infos, annot_types, sep=','):
    """_summary_

    Args:
        manifest_name (_type_): _description_
        manifest_file (_type_): _description_
        coordinate_fields (_type_): _description_
        manifest_out (_type_): _description_
        annots (_type_): _description_
        annot_infos (_type_): _description_
        sep (str, optional): _description_. Defaults to ','.
    """
    # manifest_file = "/users/mamana/H3Africa_2019_20037295_B1_100.csv"
    # manifest_name = "Illumina HumanOmni5Exome"
    # coordinate_fields = "Chr,MapInfo,SNP"
    # sep = ","
    # manifest_out = "HumanOmni5Exome-4-v1-1-B"
    # annots = "IlmnID,IlmnStrand"
    # annot_infos = "Illumina ID, Illumina Strand"

    data = [it.replace('"', '').strip().split(sep) for it in open(
        manifest_file).readlines() if not it.startswith('#')]
    out_bed = open(f'{manifest_out}.bed', 'w')
    out_bed.writelines('\\t'.join(['ID', 'CHROM', 'POS', 'AA'])+'\\n')
    # out_annot = open(f'{manifest_out}.annot.csv', 'w')
    out_hdr = open(f'{manifest_out}.hdr.csv', 'w')
    out_annot = {}

    # These are required
    coordinate_fields = coordinate_fields.split(',')
    chrm_string = coordinate_fields[0]
    pos_string = coordinate_fields[1]
    alleles = coordinate_fields[2]

    annotations = annots.split(',')
    annot_infos = [ it.strip() for it in annot_infos.split(',') ]
    annot_types = [ it.strip() for it in annot_types.split(',') ]
    header = False
    section = False
    annot_idx = {}
    datas = []
    n = 0
    myHeader = []
    manifest_name = manifest_name.strip().replace(' ', '-')
    for it in data:
        if '[Assay]' in it:
            section = True
        if chrm_string in it and pos_string in it:
            header = True
            chr_idx = it.index(chrm_string)
            pos_idx = it.index(pos_string)
            alleles_idx = it.index(alleles)
            myHeader = ['CHROM', 'POS', 'POS', f'IS_{manifest_name}']
            info = f'##INFO=<ID=IS_{manifest_name},Number=1,Type=Integer,Description=\"1 if site on {manifest_name} chip\">'
            out_hdr.writelines(info+'\\n')
            for annot in annotations:
                if annot in it:
                    if annot not in annot_idx:
                        annot_idx[annot] = it.index(annot)
                        myHeader.append(f'{manifest_name}_{annot}')
                info = f'##INFO=<ID={manifest_name}_{annot},Number=1,Type={annot_types[annotations.index(annot)]},Description=\"{annot_infos[annotations.index(annot)]}\">'
                out_hdr.writelines(info+'\\n')
        if '[Controls]' in it:
            header = False
        if header and section:
            try:
                chr = it[chr_idx]
                if chrm_string not in chr:
                    if chr not in out_annot:
                        out_annot[chr] = open(f'{manifest_out}_{chr}_annots.csv', 'w')
                        out_annot[chr].writelines('\\t'.join(myHeader)+'\\n')
                pos = it[pos_idx]
                alleles = it[alleles_idx]
                id = f'{chr}_{pos}'
                for char in ['[', ']', '/']:
                    alleles = alleles.replace(char, '')
                if 'Chr' not in id: 
                    out_bed.writelines('\\t'.join([id, chr, pos, alleles])+'\\n')
                dat = [id, chr, pos, alleles]
                if chr != chrm_string:
                    # Get annotation fields
                    myData = [chr, pos, pos, '1']
                    for annot in annotations:
                        myData.append(it[annot_idx[annot]])
                    out_annot[chr].writelines('\\t'.join(myData)+'\\n')
            except:
                print("Wrong", it)
    out_bed.close()
    out_hdr.close()
    for chr in out_annot:
        out_annot[chr].close()

if __name__ == '__main__':
    parse_manifest(args.manifest_name, args.manifest_file,
                   args.coordinate_fields, args.manifest_out, args.annots, args.annot_infos, args.annot_types, args.sep)

