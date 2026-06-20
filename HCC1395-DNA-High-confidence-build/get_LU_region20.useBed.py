import pandas as pd
import numpy as np
from scipy import stats
import os
import sys
import subprocess
import io
import argparse
import re
import math
import warnings
import multiprocessing

parser = argparse.ArgumentParser(description='Ffilter VAF')
parser.add_argument('-i','--bed_input', help='Input Merge.bed')
parser.add_argument('-v','--variant_type',help='The type of variant, SNV or Indel')
parser.add_argument('-o','--output_dir',help='The path of output dir')

args = parser.parse_args()

bed_input = args.bed_input
variant_type = args.variant_type
output_dir = args.output_dir

if not os.path.exists(output_dir):
    os.makedirs(output_dir,exist_ok=True)

#输出bed文件：
bed_file = pd.read_csv(bed_input,sep='\t')
bed_file['start'] = bed_file['POS'].apply(lambda x:int(x)-21)
bed_file['end'] = bed_file['POS'].apply(lambda x:int(x)+20)
bed_out= bed_file[bed_file['Confidence_new'].isin(['LowConf','Unclassified'])][['#CHROM','start','end']]

print(bed_out.shape)
row1 = ['chr6',0,59800000]
row2 = ['chr16',36800000,90338345]
bed_out.loc[len(bed_out)] = row1
bed_out.loc[len(bed_out)] = row2
bed_out.drop_duplicates().to_csv(output_dir+'/'+'UnConf.{0}.bed'.format(variant_type),sep='\t',index=None,header=None)
