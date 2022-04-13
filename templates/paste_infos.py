#!/usr/bin/env python3.7
'''

'''
import argparse,sys,time
import pandas as pd
import numpy as np
import csv

parser = argparse.ArgumentParser()
parser.add_argument("--infos", default="${infos}", help="")
parser.add_argument("--out",
                    default="${out}", help="")
args = parser.parse_args()


def paste_infos(infos, out ):
    
    out = open(f'{out}', 'w')
    infos = [it.strip() for it in infos.split(',')]

    datas = pd.read_csv(infos[0], delimiter='\\\\s+',
                    quotechar='\\"', engine='python')
    
    for info in infos[1:]:
        info1 = pd.read_csv(info, delimiter='\\\\s+',
                        quotechar='\\"', engine='python')
        datas = pd.merge(datas, info1, how='outer', on=['CHROM', 'POS', 'POS.1', 'ID'])

    datas.to_csv(out, sep='\\t', index=False)

if __name__ == '__main__':
    paste_infos(args.infos, args.out)


