#!/usr/bin/env python

import pandas as pd

"""
Reads map file, chunk size
Returns file with chromosome chunk_start chunk_end
"""


def chunk_split(map_file, output, chunk_size, chrms='', chunk=''):
    '''
    Return: chunk files in the output folder
    '''
    # data = [ [dat.split('\\t')[0], dat.split('\\t')[1]] for dat in open(map_file).readlines() ]
    data = []
    chunk_size = int(chunk_size)
    for df in pd.read_csv(map_file, sep='\\s+', engine='python', error_bad_lines=False, iterator=True, chunksize=chunk_size/100, names=["CHROM", "POS"]):
        chrm_datas = df.groupby('CHROM')
        for idx in range(len(list(chrm_datas))):
            chrm = list(chrm_datas)[idx][0]
            d = list(chrm_datas)[idx][1]
            data.append([chrm, d["POS"].max()])
    datas = {}
    out = open(output, 'w')
    for dat in data:
        chrm = dat[0]
        myPos = int(dat[1])
        if chrm not in datas:
            datas[chrm] = []
        datas[chrm].append(myPos)
    data = {}
    if chrms != '':
        chrms = sorted([it for it in set(chrms.split(','))])
    else:
        try:
            chrms = sorted([it for it in datas])
        except TypeError as e:
            print("TypeError: '<' not supported between instances of 'str' and 'int'")
            chrms = [it for it in datas]
            

    max_ = {}
    min_ = {}
    myPos = {}
    for chrm in chrms:
        max_[chrm] = max(datas[chrm]) + (max(datas[chrm]) % 10)
        # min_[chrm] = min(datas[chrm]) - (min(datas[chrm]) % 10) + 1
        myPos[chrm] = list(range(1, max_[chrm], chunk_size))
    for chrm in myPos:
        for pos in myPos[chrm]:
            start_ = pos
            end_ = start_ + chunk_size - 1
            out.writelines(','.join([str(chrm), str(start_), str(end_)]) + '\\n')
    out.close()


if __name__ == '__main__':
    mapFile = "${mapFile}"
    outputFile = "${chunkFile}"
    chunk_size = "${chunk_size}"
    chromosomes = "${chromosomes}"
    chunk = "${chunk}"
    chunk_split(mapFile, outputFile, chunk_size, chromosomes, chunk)
