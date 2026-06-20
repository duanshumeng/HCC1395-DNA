#Calculate VAF of Strelka
import pandas as pd
import subprocess
import argparse
import os

parser = argparse.ArgumentParser(description='FIlter RNA editing')
parser.add_argument('-i','--vcf_input', help='Input VCFs, *.vcf or *.vcf.gz')

#parser.add_argument('-o','--output_dir',help='The path of output dir')

args = parser.parse_args()

vcf_input = args.vcf_input

vcf_out = vcf_input.replace('.vcf.gz','.csv').replace('.vcf','.csv')

def VAF_calculate(alt,v_info):
    info_lst = v_info.split(':')
    print(alt)
    print(v_info)
    total_count=info_lst[0]
    if alt=='A':
        altCount = info_lst[4].split(',')[0]
    elif alt=='C':
        altCount = info_lst[5].split(',')[0]
    elif alt=='G':
        altCount = info_lst[6].split(',')[0]
    elif alt=='T':
        altCount = info_lst[7].split(',')[0]
    else:
        altCount = info_lst[2].split(',')[0]
    print(altCount)
    print(total_count)
    vaf=0
    if float(total_count) != 0:
       vaf = float(altCount)/float(total_count)
    return vaf


def get_vcf_info(vcf_input):
    if vcf_input.endswith('.vcf'):
        #取出VCF标题所在的行的位置
        A = subprocess.check_output("cat %s | grep -n '#CHROM' | awk -F ':' '{print $1}'" % (vcf_input),shell=True)

        #vcf_header = subprocess.check_output("cat %s | grep -n '#' " % (vcf_input),shell=True)
        

        #跳过前几行读取vcf文件
        vcf_file = pd.read_table(vcf_input,skiprows=int(A)-1)
    elif vcf_input.endswith('.vcf.gz'):
        A = subprocess.check_output("zcat %s | grep -n '#CHROM' | awk -F ':' '{print $1}'" % (vcf_input),shell=True)

        #vcf_header = subprocess.check_output("zcat %s | grep -n '##'" % (vcf_input),shell=True)

        vcf_file = pd.read_table(vcf_input,skiprows=int(A)-1,compression='gzip')

    vcf_file['VAF'] = vcf_file.apply(lambda x:VAF_calculate(x['ALT'],x['TUMOR']),axis=1)
    vcf_file['tag']=vcf_file.apply(lambda x:'&'.join([str(x['#CHROM']),str(x['POS']),str(x['REF']),str(x['ALT'])]),axis=1)    

    return vcf_file


vcf_sub = get_vcf_info(vcf_input)

vcf_sub.to_csv(vcf_out,index=None)

