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
parser.add_argument('-n','--num',type=int,default=9,help='Level1 cutoff')
parser.add_argument('-o','--output_dir',help='The path of output dir')

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
    caller = vcf_input.split('.')[1]

    pattern = r'(\.\/)?\.:+'
    aim_header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    for i in vcf_file.columns:
        if i.endswith('.T'):
            aim_header.append(i)
    aim_header = aim_header+['Total sum','Confidence']
    #仅保留1-22号染色体的信息
    chrs = ['chr'+str(i) for i in range(1,23)] + ['chrX','chrY']
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

# print('Scoring....')
# def classify_2callers(total_score):
#     if total_score > 5:
#         return 'HighConf'
#     elif 4 < total_score <= 5:
#         return 'MedConf'
#     elif 2 < total_score <= 4 :
#         return 'LowConf'
#     else:
#         return 'Unclassified'

# def classify_3callers(total_score):
#     if total_score > num:
#         return 'HighConf'
#     elif num-2 < total_score <= num:
#         return 'MedConf'
#     elif num-5 < total_score <= num-2:
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


# def classify_3callers(total_score):
#     if total_score > 8:
#         return 'HighConf'
#     elif 6 < total_score <= 8:
#         return 'MedConf'
#     elif 3 < total_score <= 6:
#         return 'LowConf'
#     else:
#         return 'Unclassified'



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
#输出HighConf/MedConf
vcf_sub_result[vcf_sub_result['Confidence_new'].isin(['HighConf','MedConf'])][['#CHROM', 'POS','POS', 'ID', 'REF', 'ALT','Confidence_new']].drop_duplicates().to_csv(output_dir+'/'+'Merge.HC.{0}.{1}.bed'.format(variant_type,num),sep='\t',index=None)


#输出中间文件:
print('Output intermediate file Merge_Indels.txt......')
vcf_sub_result[aim_header].to_csv(output_dir+'/'+'Merge_{}.tmp.txt'.format(variant_type),sep='\t')

#仅保留HighConf和MedConf的位点进行信息融合
vcf_sub_result = vcf_sub_result[vcf_sub_result['Confidence_new'].isin(['HighConf','MedConf'])]

print('Combining INFO.......')
#合并目标信息并重新组合信息
def combine_info(df,tag):
    try:
        print('This is  ',tag)
        df_sub = df[df['tag'] == tag][['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL']].head(1)
        print(df_sub)
        #print(df[(df['tag'] == tag)]['INFO'])
        #print(df[((df['tag'] == tag) & (df['INFO'].str.contains('SOMATIC')))])
        info = df[((df['tag'] == tag) & (df['INFO'].str.contains('SOMATIC')))]['INFO'].values[0]
        MQ=info.split(';')[7]
        MQ0=info.split(';')[8]

        smps = [i for i in df.columns if i.endswith('.T')] 
        DPs = df[((df['tag'] == tag) & (df['INFO'].str.contains('SOMATIC')))][smps].apply(lambda x:x.str.split(':').str[0].str.split(',').str[0]).apply(lambda x: ','.join(x.replace('.',0).astype(str)), axis=1).values[0]
        #print(DPs)
        df['caller sum'] = df.apply(lambda x:str(x['caller']+'_n')+'='+str(int(x['Total sum']))+',13',axis=1)
        call_score = ';'.join(list(df[df['tag'] == tag]['caller sum']))
        if len(DPs) > 0:
            #DPs_list = [int(i) for i in DPs.str.split(',').sum() if int(i) != 0]
            DPs_list = [int(i) for i in DPs.split(',') if int(i) != 0]
            #print(DPs_list)
            if len(DPs_list) > 0:
                DP_mean = np.mean(DPs_list)
                DP_std = np.std(DPs_list)
                try:
                    DP95 = stats.norm.interval(0.95, loc=DP_mean, scale=DP_std)
                    DP95_rod = tuple(round(x) for x in DP95)
                except Exception as e:
                    DP95 = round(DP_mean)
                    DP95_rod = (round(DP_mean),round(DP_mean))
            else:
                DP_mean = 0
                DP95_rod = ('.','.')
        else:
                DP_mean = 0
                DP95_rod = ('.','.')


        VAFs = df[((df['tag'] == tag) & (df['INFO'].str.contains('ECNT')))][smps].apply(lambda x:x.str.split(':').str[2]).apply(lambda x: ','.join(x.replace('.',0).astype(str)), axis=1)
        #print(VAFs)
        if len(VAFs) > 0:
            VAFs_list = [float(i) for i in VAFs.str.split(',').sum() if float(i) != 0]
            #print(VAFs_list)
            if len(DPs_list) > 0:
                VAF_mean = np.mean(VAFs_list)
                VAF_std = np.std(VAFs_list)
                try:
                    VAF95 = stats.norm.interval(0.95, loc=VAF_mean, scale=VAF_std)
                    VAF95_rod = tuple(round(x,3) for x in VAF95)
                except Exception as e:
                    VAF95 = round(VAF_mean,3)
                    VAF95_rod = (round(VAF_mean,3),round(VAF_mean,3))
            else:
                VAF_mean = 0
                VAF95_rod = ('.','.')
        else:
                VAF_mean = 0
                VAF95_rod = ('.','.')

    
        all_info = str(MQ)+';'+str(MQ0)+';'+'DP='+str(int(np.ceil(DP_mean)) if not math.isnan(DP_mean) else 0)+';'+ \
            'DP95='+str(DP95_rod[0])+','+str(DP95_rod[1])+';'+'VAF='+str(round(VAF_mean,3))+';'+ \
                'VAF95='+str(VAF95_rod[0])+','+str(VAF95_rod[1])+';'+call_score
        df_sub['FILTER'] = df[df['tag'] == tag]['Confidence_new'].values[0]
        df_sub['INFO'] = all_info
        return df_sub
    except Exception as e:
        column_names = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL','FILTER','INFO']
        df_sub = pd.DataFrame(columns=column_names)
        return df_sub


#使用多进程处理，每条染色体一个进程
chrs = ['chr'+str(i) for i in range(1,23)]+ ['chrX','chrY']
#用来测试的染色体
#chrs = ['chr1']
#"""
for i in chrs:
    print(i)
    globals()[i] = vcf_sub_result[vcf_sub_result['#CHROM']==i]['tag'].drop_duplicates().values
#tags = vcf_sub_result['tag'].drop_duplicates().values

def get_work(df_lst,vcf_sub_result,tags):
    # 忽略 RuntimeWarning 类型的警告
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    for t in tags:
        #print(t)
        try:
            df_res = combine_info(vcf_sub_result,t)
            df_lst.append(df_res)
        except Exception as e:
            print(e)
    # 取消忽略警告设置
    warnings.filterwarnings("default", category=RuntimeWarning)

manager = multiprocessing.Manager()
df_lst = manager.list()

num_processes = len(chrs)
# 创建进程池
pool = multiprocessing.Pool(processes=num_processes)

# 在多个进程中并行执行操作
pool_args = [(df_lst, vcf_sub_result,globals()[chr_name]) for chr_name in chrs]
print(pool_args)
pool.starmap(get_work, pool_args)

# 关闭进程池
pool.close()
pool.join()

#"""

# tags = vcf_sub_result['tag'].drop_duplicates().values
# df_lst = []
# def get_work(df_lst,vcf_sub_result,tags):
#     # 忽略 RuntimeWarning 类型的警告
#     warnings.filterwarnings("ignore", category=RuntimeWarning)
#     for t in tags:
#         print(t)
#         df_res = combine_info(vcf_sub_result,t)
#         df_lst.append(df_res)

#     return df_lst
#     # 取消忽略警告设置
#     warnings.filterwarnings("default", category=RuntimeWarning)

# df_lst = get_work(df_lst,vcf_sub_result,tags)

print('Combining '+'dataframe......')
df_final = pd.concat(df_lst)


#重新制作VCF的header信息
vcf_header = '''##fileformat=VCFv4.2
##fileDate=20230802
##reference=GRCh38.d1.vd1
##FILTER=<ID=HighConf,Description="highly confident that it is a real somatic mutation">
##FILTER=<ID=MedConf,Description="confident that it is a real somatic mutation">
##FILTER=<ID=LowConf,Description="not very confident that it is a real somatic mutation">
##FILTER=<ID=Unclassified,Description="likely not a real somatic mutation">
##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Average depth of tumor samples">
##INFO=<ID=DP95,Number=2,Type=Integer,Description="Estimated 95% confidence interval for DP in tumor samples">
##INFO=<ID=VAF,Number=1,Type=Float,Description="Average VAF of tumor samples">
##INFO=<ID=VAF95,Number=2,Type=Float,Description="Estimated 95% confidence interval for VAF in tumor samples">
##INFO=<ID=strelka_n,Number=2,Type=Integer,Description="The number of files support the variant in strelka, total 13">
##INFO=<ID=TNseq_n,Number=2,Type=Integer,Description="The number of files support the variant in TNseq, total 13">
##INFO=<ID=TNscope_n,Number=2,Type=Integer,Description="The number of files support the variant in TNscope, total 13">
##INFO=<ID=muTect2_n,Number=2,Type=Integer,Description="The number of files support the variant in TNseq, total 13">
##INFO=<ID=somaticSniper_n,Number=2,Type=Integer,Description="The number of files support the variant in TNscope, total 13">
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr4,length=190214555>
##contig=<ID=chr5,length=181538259>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr12,length=133275309>
##contig=<ID=chr13,length=114364328>
##contig=<ID=chr14,length=107043718>
##contig=<ID=chr15,length=101991189>
##contig=<ID=chr16,length=90338345>
##contig=<ID=chr17,length=83257441>
##contig=<ID=chr18,length=80373285>
##contig=<ID=chr19,length=58617616>
##contig=<ID=chr20,length=64444167>
##contig=<ID=chr21,length=46709983>
##contig=<ID=chr22,length=50818468>
##contig=<ID=chrX,length=156040895>
##contig=<ID=chrY,length=57227415>
'''

print('Writing to file......')
out_file = output_dir+'/'+'high-confidence_s{}_v1.vcf'.format(variant_type)
with open(out_file, 'w') as w:
    w.writelines(vcf_header)

df_final.to_csv(output_dir+'/'+'All_s{}_v1.vcf'.format(variant_type),index=False,sep='\t')

df_out = df_final[df_final['FILTER'].isin(['HighConf','MedConf'])]
df_out['FILTER'] = 'PASS;'+df_out['FILTER']
df_out.to_csv(out_file,mode='a', index=False,sep='\t')

print('Process successed!!!')