#!/usr/bin/env python3.7
'''

'''
import argparse,sys,time
import pandas as pd
import numpy as np
import csv

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="${label}", help="")
parser.add_argument("--annots", default="${annots}", help="")
parser.add_argument("--annot_infos", default="${annot_infos}", help="")
parser.add_argument("--out",
                    default="${out}", help="")
args = parser.parse_args()


def parse_manifest(dataset, annots, annot_infos, out ):
    
    out_hdr = open(f'{out}.hdr', 'w')
    annots = annots.split(',')
    annot_infos = [ it.strip() for it in annot_infos.split(',') ]

    for annot in annots:
        info = f'##INFO=<ID={dataset}_{annot.replace("INFO/","")},Number=1,Type=Float,Description=\"{annot_infos[annots.index(annot)]}\">'
        print(info)
        out_hdr.writelines(info+'\\n')

    out_hdr.close()

if __name__ == '__main__':
    parse_manifest(args.dataset, args.annots, args.annot_infos, args.out)


