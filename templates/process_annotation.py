#!/usr/bin/env python3.7
'''

'''
import argparse,sys,time
import pandas as pd
import numpy as np
import csv

parser = argparse.ArgumentParser()
parser.add_argument("--annot_file", default="${annot_file}", help="One TSV annotation files")
parser.add_argument(
    "--prefix_out", default="${prefix_out}", help="TSV annotation file")
parser.add_argument("--annotation", default="${annotation}", help="INFO annotation")
args = parser.parse_args()


def readTSV(annot_file, prefix_out='', annotation=''):
    '''
    :param inBed:
    :param outBed:
    :return:
    '''
        
    # Split into chromosome
    out_datas = {}
    header = []
    col_out = open(f'{prefix_out}.columns', 'w')
    for df in pd.read_csv(annot_file, sep='\\\\s+', engine='python', error_bad_lines=False, iterator=True, chunksize=1000000):
        chrm_datas = df.groupby('CHROM')
        for idx in range(len(list(chrm_datas))):
            chrm = list(chrm_datas)[idx][0]
            data = list(chrm_datas)[idx][1]
            data = data.replace(np.nan, '0.0', regex=True)
            data = data.replace('.', '0.0')
            if len(header) == 0:
                header = [ it.replace('POS.1', 'POS') for it in data.columns.tolist() ]
                col_out.writelines(','.join(header)+'\\n')
            data = data.values.tolist()
            if chrm not in out_datas:
                # Writing to file
                out_datas[chrm] = open(f'{prefix_out}__{chrm}.tsv', 'w')
                out_datas[chrm].writelines('\\t'.join(header)+'\\n')
            for record in data:
                # Writing to file
                out_datas[chrm].writelines('\\t'.join( [ str(it) for it in record ] )+'\\n')    
    for chrm in out_datas:
        out_datas[chrm].close()
    col_out.close()    
    
    #  Writing header
    # info = pd.DataFrame([], columns=['info'])
    # for col in header:
    #     if col not in ['ID', 'POS', 'CHROM', 'FROM', 'TO']:
    #         info = info.append({'info': f'##INFO=<ID={col},Number=1,Type=Float,Description=\"{annotation}\">'}, ignore_index=True)
    # info.to_csv(annot_hdr_file, sep=';', index=False, header=False,
    #             quoting=csv.QUOTE_NONE, quotechar='\\0')
    
if __name__ == '__main__':
    readTSV(args.annot_file, args.prefix_out, args.annotation)


