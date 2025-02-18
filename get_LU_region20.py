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
parser.add_argument('-i','--vcf_input', help='Input VCFs,separate by ,. such as "Merge.strelka.Indels.filter.vcf,Merge.strelka.Indels.filter.vcf" ')
parser.add_argument('-t','--vcf_type',help='SEQC2 or PGx')
parser.add_argument('-c','--caller_num',help='The number of variant callers 2 / 3...')
parser.add_argument('-v','--variant_type',help='The type of variant, SNV or Indel')
parser.add_argument('-o','--output_dir',help='The path of output dir')
parser.add_argument('-n','--num',type=int,default=9,help='Level1 cutoff')

args = parser.parse_args()

vcf_input = args.vcf_input
vcf_type = args.vcf_type
caller_num = int(args.caller_num)
variant_type = args.variant_type
output_dir = args.output_dir
num = int(args.num)

if not os.path.exists(output_dir):
    os.makedirs(output_dir,exist_ok=True)

def get_vcf_info(vcf_input):
    if vcf_input.endswith('.HCR.vcf'):
        vcf_file = pd.read_csv(vcf_input,sep='\t')
    else:
        #取出VCF标题所在的行的位置
        A = subprocess.check_output("cat %s | grep -n '#CHROM' | awk -F ':' '{print $1}'" % (vcf_input),shell=True)

        vcf_header = subprocess.check_output("cat %s | grep -n '##' " % (vcf_input),shell=True)

        #跳过前几行读取vcf文件
        vcf_file = pd.read_table(vcf_input,skiprows=int(A)-1)
    
    caller = vcf_input.split('.')[1]

    pattern = r'(\.\/)?\.:+'
    aim_header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    for i in vcf_file.columns:
        if i.endswith('.T'):
            aim_header.append(i)
    aim_header = aim_header+['Total sum','Confidence']
    #仅保留1-22号染色体的信息
    chrs = ['chr'+str(i) for i in range(1,23)]
    vcf_file_filter = vcf_file[vcf_file['#CHROM'].isin(chrs)]
    vcf_sub = vcf_file_filter[aim_header]
    vcf_sub['tag']=vcf_sub.apply(lambda x:'&'.join([str(x['#CHROM']),str(x['POS']),str(x['REF']),str(x['ALT'])]),axis=1)    
    vcf_sub['caller'] = caller

    return vcf_sub

vcf_lst = []


for vcf1 in vcf_input.split(','):
    print('Reading '+vcf1+'......')
    vcf_sub1 = get_vcf_info(vcf1)
    vcf_lst.append(vcf_sub1)


vcf_sub = pd.concat(vcf_lst)
score_dict = {'Level1':3,'Level2':2,'Level3':1,'Level4':0}
vcf_sub['Score'] = vcf_sub['Confidence'].replace(score_dict)
print(vcf_sub.head())
vcf_sub_filter = vcf_sub.groupby('tag').filter(lambda x:len(x)>1)
vcf_sub_score =  vcf_sub_filter.groupby('tag')['Score'].sum().reset_index()
vcf_sub_result = pd.merge(vcf_sub_filter, vcf_sub_score, on='tag', suffixes=('', '_Sum'))

print('Scoring....')
def classify_2callers(total_score):
    if total_score > 5:
        return 'HighConf'
    elif 4 < total_score <= 5:
        return 'MedConf'
    elif 2 < total_score <= 4:
        return 'LowConf'
    else:
        return 'Unclassified'


# def classify_3callers(total_score):
#     if total_score > 8:
#         return 'HighConf'
#     elif 6 < total_score <= 8:
#         return 'MedConf'
#     elif 3 < total_score <= 6:
#         return 'LowConf'
#     else:
#         return 'Unclassified'
def classify_3callers(total_score):
    if total_score > num:
        return 'HighConf'
    elif num-2 < total_score <= num:
        return 'MedConf'
    elif num-4 < total_score <= num-2:
        return 'LowConf'
    else:
        return 'Unclassified'

if caller_num == 2:
    vcf_sub_result['Confidence_new'] = vcf_sub_result['Score_Sum'].apply(lambda x:classify_2callers(x))
elif caller_num == 3:
    vcf_sub_result['Confidence_new'] = vcf_sub_result['Score_Sum'].apply(lambda x:classify_3callers(x))
aim_header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
for i in vcf_sub_result.columns:
            if i.endswith('.T'):
                        aim_header.append(i)
aim_header = aim_header+['Confidence_new','tag','Score_Sum','caller','Total sum']

#输出bed文件：
vcf_sub_result[['#CHROM', 'POS','POS', 'ID', 'REF', 'ALT','Confidence_new','Score_Sum','caller','Total sum']].drop_duplicates().to_csv(output_dir+'/'+'Merge.{0}.{1}.bed'.format(variant_type,num),sep='\t',index=None)

#输出bed文件：
bed_file = vcf_sub_result[['#CHROM', 'POS', 'ID', 'REF', 'ALT','Confidence_new','Score_Sum','caller','Total sum']].drop_duplicates()
bed_file['start'] = bed_file['POS'].apply(lambda x:int(x)-21)
bed_file['end'] = bed_file['POS'].apply(lambda x:int(x)+20)
bed_out= bed_file[bed_file['Confidence_new'].isin(['LowConf','Unclassified'])][['#CHROM','start','end']]
row1 = ['chr6',0,59800000]
row2 = ['chr16',36800000,90338345]
bed_out.loc[len(bed_out)] = row1
bed_out.loc[len(bed_out)] = row2
bed_out.drop_duplicates().to_csv(output_dir+'/'+'UnConf.{0}.bed'.format(variant_type),sep='\t',index=None,header=None)
