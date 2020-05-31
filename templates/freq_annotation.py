#!/usr/bin/env python3.7
'''

'''
import argparse,sys,time
import pandas as pd
import numpy as np
import csv

parser = argparse.ArgumentParser()
parser.add_argument("--inTSV", default="${inTSV}", help="One TSV annotation files")
parser.add_argument("--outAnnot", default="${outAnnot}", help="TSV annotation file")
parser.add_argument("--outHdr", default="${outHdr}", help="TSV annotation header file")
parser.add_argument("--sites", default="${sites}", help="Base annotation file")
parser.add_argument("--annot", default="${annot}", help="INFO annotation")
args = parser.parse_args()


def readTSV(inTSV, sites='', outAnnot='', outHdr='', annot=''):
    '''
    :param inBed:
    :param outBed:
    :return:
    '''
    
    base = pd.read_csv(sites, sep='\\\\s+', engine='python', error_bad_lines=False)
    data = pd.read_csv(inTSV, sep='\\\\s+', engine='python', error_bad_lines=False)
    base = pd.merge(base, data, on='rsID', how='left')
    base = base.replace(np.nan, '0.0', regex=True)
    base = base.replace('.', '0.0')

    datas = base.to_dict()  # Transform to dict
    new_datas = {}
    for col in datas:
        if '_MAF_FREQ' in col or 'rsID' in col:
            new_datas[col] = {}
            for idx in datas[col]:
                if '_MAF_FREQ' in col:
                    frq = str(datas[col][idx]).strip().split(',')[0]
                    if frq == '.':
                        frq = '0.0'
                    try:
                        frq = float(frq)
                    except:
                        print(frq)
                    new_datas[col][idx] = format(frq, '.5f')
                if 'rsID' in col:
                    new_datas[col][idx] = datas[col][idx]
                    
    base = pd.DataFrame.from_dict(new_datas)

    chrom = base['rsID'].str.split('_').str[0]
    pos = base['rsID'].str.split('_').str[1]
    base = base.drop(columns=['rsID'])
    base.insert(0, 'POS', pos, allow_duplicates=True)
    base.insert(0, 'POS', pos, allow_duplicates=True)
    base.insert(0, 'CHROM', chrom)

    # Writing to file
    base.to_csv(outAnnot, sep='\\t', index=False)

    info = pd.DataFrame([], columns=['info'])
    for col in base.columns:
        if '_MAF_FREQ' in col:
            info = info.append({'info': "##INFO=<ID={},Number=1,Type=Float,Description=\"{}\">".format(col, annot)}, ignore_index=True)

    #  Writing header
    info.to_csv(outHdr, sep=';', index=False, header=False, quoting=csv.QUOTE_NONE, quotechar='\\0')
    
if __name__ == '__main__':
    readTSV(args.inTSV, args.sites, args.outAnnot, args.outHdr, args.annot)


