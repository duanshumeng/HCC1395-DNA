# -*- coding: UTF-8 -*
#!/opt/local/cobweb/envs/hrd/bin/python3.
import pandas as pd
import os
import sys
import subprocess
import io
import argparse

parser = argparse.ArgumentParser(description='Calculate VAF of TNseq or Haplotyper')
parser.add_argument('-i', '--vcf_input', help='Path to the VCF files, includes *.vcf.gz.')
parser.add_argument('-c', '--caller', help='TNseq/Haplotyper')
parser.add_argument('-R','--region',required=False,help='指定区间的bed文件')
args = parser.parse_args()



vcf_input = args.vcf_input
caller = args.caller
region = args.region

if vcf_input.endswith('.vcf'):
    out_file = vcf_input.replace('.vcf','.VAF.bed')

    #取出VCF标题所在的行的位置
    A = subprocess.check_output("cat %s | grep -n '#CHROM' | awk -F ':' '{print $1}'" % (vcf_input),shell=True)
    print(A)
    #跳过前几行读取vcf文件
    vcf_file = pd.read_table(vcf_input,skiprows=int(A)-1)
elif vcf_input.endswith('.vcf.gz'):
    out_file = vcf_input.replace('.vcf.gz','.VAF.bed')

    A = subprocess.check_output("zcat %s | grep -n '#CHROM' | awk -F ':' '{print $1}'" % (vcf_input),shell=True)

    vcf_file = pd.read_table(vcf_input,skiprows=int(A)-1,compression='gzip')

# if len([i for i in vcf_file.columns if i.endswith(('_T','_T_1'))]) > 0:
#     Tumor = [i for i in vcf_file.columns if i.endswith(('_T','_T_1'))][0]
#     if caller == 'TNseq':
#         vcf_file['VAF'] = vcf_file[Tumor].apply(lambda x:x.split(':')[2])
#     elif caller == 'Haplotyper':
#         vcf_file['VAF'] = vcf_file[Tumor].apply(lambda x:int(x.split(':')[1].split(',')[1])/int(x.split(':')[1].split(',')[0]) if int(x.split(':')[1].split(',')[0]) != 0 else 0)
# else:
#     Normal = [i for i in vcf_file.columns if i.endswith(('_N','_N_1'))][0]
#     if caller == 'TNseq':
#         vcf_file['VAF'] = vcf_file[Normal].apply(lambda x:x.split(':')[2])
#     elif caller == 'Haplotyper':
#         vcf_file['VAF'] = vcf_file[Normal].apply(lambda x:int(x.split(':')[1].split(',')[1])/int(x.split(':')[1].split(',')[0]) if int(x.split(':')[1].split(',')[0]) != 0 else 0)

print([i for i in vcf_file.columns if ('HCC1395' in i and 'HCC1395BL' not in i)])

if len([i for i in vcf_file.columns if ('HCC1395' in i and 'HCC1395BL' not in i)]) > 0:
    Tumor = [i for i in vcf_file.columns if ('HCC1395' in i and 'HCC1395BL' not in i)][0]
    if caller == 'TNseq':
        vcf_file['VAF'] = vcf_file[Tumor].apply(lambda x:x.split(':')[2])
    elif caller == 'Haplotyper':
        vcf_file['VAF'] = vcf_file[Tumor].apply(lambda x:int(x.split(':')[1].split(',')[1])/int(x.split(':')[1].split(',')[0]) if int(x.split(':')[1].split(',')[0]) != 0 else 0)
else:
    Normal = [i for i in vcf_file.columns if ('HCC1395BL' in i)][0]
    if caller == 'TNseq':
        vcf_file['VAF'] = vcf_file[Normal].apply(lambda x:x.split(':')[2])
    elif caller == 'Haplotyper':
        vcf_file['VAF'] = vcf_file[Normal].apply(lambda x:int(x.split(':')[1].split(',')[1])/int(x.split(':')[1].split(',')[0]) if int(x.split(':')[1].split(',')[0]) != 0 else 0)

vcf_file[['#CHROM', 'POS','POS','REF', 'ALT','VAF']].to_csv(out_file,index=None,sep='\t',header=False)
vcf_df = vcf_file[['#CHROM', 'POS','REF', 'ALT','VAF']]

hrc_VAF = subprocess.check_output("bedtools intersect -b {0} -a {1}".format(region,out_file),
                                  shell=True)

hrc_VAF = hrc_VAF.decode('utf-8')

hrc_df = pd.read_table(io.StringIO(hrc_VAF),lineterminator='\n',header=None)
hrc_df.columns = ['#CHROM', 'POS','POS_1','REF', 'ALT','VAF']
hrc_df = hrc_df[['#CHROM', 'POS','REF', 'ALT','VAF']].drop_duplicates()
hrc_df['POS'] = hrc_df['POS'].apply(lambda x:str(x))
hrc_df['Type'] = 'HRC'
hrc_df['tag'] = hrc_df.apply(lambda x:'_'.join([x['#CHROM'],x['POS'],x['REF'],x['ALT']]),axis=1)

vcf_df['POS'] = vcf_df['POS'].apply(lambda x:str(x))

vcf_df['tag'] = vcf_df.apply(lambda x:'_'.join([x['#CHROM'],x['POS'],x['REF'],x['ALT']]),axis=1)

#vcf_df.apply(lambda x:'HRC' if x['tag'] in list(hrc_df['tag']) else 'Other',axis=1).drop_duplicates()
vcf_df['Type'] = vcf_df.apply(lambda x:'HRC' if x['tag'] in list(hrc_df['tag']) else 'Other',axis=1)

os.system('rm {}'.format(out_file))
vcf_df.to_csv(out_file.replace('.VAF.bed','.VAF.csv'),index=None)
