import pandas as pd
import os
import sys
import subprocess
import io
import argparse
import re


parser = argparse.ArgumentParser(description='Filter VAF')
parser.add_argument('-i', '--vcf_input', help='Path to the VCF files, includes *.vcf.gz.')
parser.add_argument('-t','--vcf_type',help='SEQC2 or PGx')
parser.add_argument('-n','--num',type=int,default=9,help='Level1 cutoff')

args = parser.parse_args()

vcf_input = args.vcf_input
vcf_type = args.vcf_type
num = int(args.num)

if vcf_input.endswith('.vcf'):
    out_file = vcf_input.replace('.vcf','.{}.filter.vcf'.format(num))

    #取出VCF标题所在的行的位置
    A = subprocess.check_output("cat %s | grep -n '#CHROM' | awk -F ':' '{print $1}'" % (vcf_input),shell=True)

    vcf_header = subprocess.check_output("cat %s | grep -n '##' " % (vcf_input),shell=True)
    

    #跳过前几行读取vcf文件
    vcf_file = pd.read_table(vcf_input,skiprows=int(A)-1)
elif vcf_input.endswith('.vcf.gz'):
    out_file = vcf_input.replace('.vcf.gz','.{}.filter.vcf'.format(num))

    A = subprocess.check_output("zcat %s | grep -n '#CHROM' | awk -F ':' '{print $1}'" % (vcf_input),shell=True)

    vcf_header = subprocess.check_output("zcat %s | grep -n '##'" % (vcf_input),shell=True)

    vcf_file = pd.read_table(vcf_input,skiprows=int(A)-1,compression='gzip')

vcf_header = vcf_header.decode('utf-8')
lines = vcf_header.split('\n')
with open(out_file, 'w') as w:
        for i, line in enumerate(lines, 1):
                    line_split = line.split(':',maxsplit=1)
                    if len(line_split) > 1:
                        w.write(line_split[1]+'\n')
pattern = r'(\.\/)*?\.:+'
aim_header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
for i in vcf_file.columns:
    #if i.endswith('.T') and ('2021' not in i):
    if i.endswith('.T'):
        aim_header.append(i)

vcf_sub = vcf_file[aim_header]
for i in vcf_sub.columns:
    if i.endswith('.T'):
        print(i)
        vcf_sub[i+' sum'] = vcf_sub[i].apply(lambda x: 0 if re.match(pattern,x) else 1)

if vcf_type == 'SEQC2':
    #对样本重复进行打分
    samples = [i for i in vcf_file.columns if i.endswith('.T')]
    groups = list(set([i.rsplit('_',1)[0] for i in samples]))
    group_new = []
    for g in groups:
        if g not in ['WGS_EA','WGS_LL','WGS_NC']:
            group_new.append(g)
            vcf_sub[g+' sum'] = vcf_sub[[i for i in vcf_sub.columns if 'sum' in i and re.search(g,i)]].sum(axis=1)
        else:
            group_new.append('WGS_EA_LL_NC')
            vcf_sub['WGS_EA_LL_NC sum'] = vcf_sub[[i+' sum' for i in samples if re.search(r'WGS_EA|WGS_LL|WGS_NC',i)]].sum(axis=1)

    #对分组进行打分
    group_new = list(set(group_new))
    vcf_sub['Total sum'] = vcf_sub[[i+' sum' for i in group_new]].sum(axis=1)

    #根据打分进行等级评定
    print("Level1 cutoff: num-5")
    def classify(total_score):
        #15~8
        if total_score > num:
            return 'Level1'
        #5~7
        elif num-3 < total_score <= num:
            return 'Level2'
        #3~4
        elif num-5 < total_score <= num-3: 
            return 'Level3'
        #1~2
        else:
            return 'Level4'
    #print("Level1 cutoff: num-7")
    # def classify(total_score):
    #     #15~8
    #     if total_score > num:
    #         return 'Level1'
    #     #5~7
    #     elif num-3 < total_score <= num:
    #         return 'Level2'
    #     #3~4
    #     elif num-7 < total_score <= num-3: 
    #         return 'Level3'
    #     #1~2
    #     else:
    #         return 'Level4'
    # print("Level1 cutoff: num-6")
    # def classify(total_score):
    #     #15~8
    #     if total_score > num:
    #         return 'Level1'
    #     #5~7
    #     elif num-3 < total_score <= num:
    #         return 'Level2'
    #     #3~4
    #     elif num-6 < total_score <= num-3: 
    #         return 'Level3'
    #     #1~2
    #     else:
    #         return 'Level4'
    vcf_sub['Confidence'] = vcf_sub['Total sum'].apply(lambda x:classify(x))
elif vcf_type == 'PGx':
    #对样本重复进行打分
    samples = [i for i in vcf_file.columns if i.endswith('.T')]
    groups = list(set([i.rsplit('_',1)[0] for i in samples]))
    group_new = []
    for g in groups:
        if g not in ['BGI_T10_WGS','BGI_T20_WGS','WGS_WUX']:
            group_new.append(g)
            vcf_sub[g+' sum'] = vcf_sub[[i for i in vcf_sub.columns if 'sum' in i and re.search(g,i)]].sum(axis=1)
        else:
            group_new.append('WGS_ARD_NVG_WGE_T7')
            vcf_sub['WGS_ARD_NVG_WGE_T7 sum'] = vcf_sub[[i+' sum' for i in samples if re.search(r'WGS_EA|WGS_LL|WGS_NC',i)]].sum(axis=1)

    #对分组进行打分
    group_new = list(set(group_new))
    vcf_sub['Total sum'] = vcf_sub[[i+' sum' for i in group_new]].sum(axis=1)

    #根据打分进行等级评定
    def classify(total_score):
        if total_score > num:
            return 'Level1'
        elif num-3 < total_score <= num:
            return 'Level2'
        elif num-5 < total_score <= num-3:
            return 'Level3'
        else:
            return 'Level4'

    # def classify(total_score):
    #     if total_score > num:
    #         return 'Level1'
    #     elif num-3 < total_score <= num:
    #         return 'Level2'
    #     elif num-7 < total_score <= num-3:
    #         return 'Level3'
    #     else:
    #         return 'Level4'

    # def classify(total_score):
    #     if total_score > num:
    #         return 'Level1'
    #     elif num-3 < total_score <= num:
    #         return 'Level2'
    #     elif num-6 < total_score <= num-3:
    #         return 'Level3'
    #     else:
    #         return 'Level4'

    vcf_sub['Confidence'] = vcf_sub['Total sum'].apply(lambda x:classify(x))


vcf_sub[aim_header+['Total sum','Confidence']].to_csv(out_file,mode='a', index=False,sep='\t')
#docker run -it --rm -v /data1/HCC1395/tumor_normal_WGS/SEQC2_VCF/bwa_VCF/Filter_by_callers/Filter_by_repeats/Filter_by_batch/BED_files/venn_HCR/Check_loss:/data/ hap.py:v1 bcftools view /data/Filter_merge.vcf -Oz -o /data/Filter_merge.vcf.gz
#docker run -it --rm -v /data1/HCC1395/tumor_normal_WGS/SEQC2_VCF/bwa_VCF/Filter_by_callers/Filter_by_repeats/Filter_by_batch/BED_files/venn_HCR/Check_loss:/data/ hap.py:v1 bcftools index /data/Filter_merge.vcf.gz
