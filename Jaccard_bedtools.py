#!/opt/local/cobweb/envs/hrd/bin/python3.9
import pandas as pd
import os
import itertools
import sys
import argparse

parser = argparse.ArgumentParser(description='Calculate the reproducibility with Jaccard Index')
parser.add_argument('-i', '--input_dir', help='Path to the input files, includes *.vcf.gz.')
parser.add_argument('-o', '--out_prefix', help='Prefix of output path')
parser.add_argument('-c', '--HCR',default='', help='如果只针对高置信区间内进行重复性比较，则添加对应的bed文件,否则就在全基因组范围内比较')
args = parser.parse_args()

input_dir = os.path.abspath(args.input_dir)
#year = sys.argv[2]
#输出文件夹的前缀
out_prefix = args.out_prefix
#caller为‘All'的时候，只针对年份进行比较
#caller = sys.argv[3]

#seqc2为yes的时候，和SEQC2的标准数据集进行比较，为no的时候两两进行比较

HCR = args.HCR


def unique_combinations(arr):
    combos = itertools.combinations(arr, 2)
    unique_combos = set(combos)
    return unique_combos

tail=''
aim_arr = []
print(input_dir)
for i in os.listdir(input_dir):
    if i.endswith('vcf.gz'):
        aim_arr.append(i)
        tail='vcf.gz'
    elif i.endswith('vcf'):
        aim_arr.append(i)
        tail = 'vcf'
        
    # if caller != 'All':
    #   if year in i and caller in i and i.endswith('.gz'):
    #       aim_arr.append(i)
    # else:
    #     if year in i and i.endswith('.gz'):
    #       aim_arr.append(i)

print(aim_arr)
out_dir = input_dir+'/'+out_prefix
os.makedirs(out_dir,exist_ok=True)
print(out_dir)
JI_dict = {'type':[],
           'JI':[],
           'source':[]}


for out in unique_combinations(aim_arr):
        if HCR != '':
            print('Compared in '+HCR)
            a1,b1 = out[0],out[1]
            print('bedtools intersect -header -a {0} -b {1} > {2}'.format(
                input_dir+'/'+a1,
                HCR,
                input_dir+'/'+a1.replace(tail,'HCR'+'.'+tail)
            ))
            os.system('bedtools intersect -header -a {0} -b {1} > {2}'.format(
                input_dir+'/'+a1,
                HCR,
                input_dir+'/'+a1.replace(tail,'HCR'+'.'+tail)
            ))
            os.system('bedtools intersect -header -a {0} -b {1} > {2}'.format(
                input_dir+'/'+b1,
                HCR,
                input_dir+'/'+b1.replace(tail,'HCR'+'.'+tail)
            ))
            a = a1.replace(tail,'HCR'+'.'+tail)
            b = b1.replace(tail,'HCR'+'.'+tail)
            name = a.replace('.'+tail,'')+'VS'+b.replace('.'+tail,'')

        else:
            a,b = out[0],out[1]
            print(a,b)
            name = a.replace('.'+tail,'')+'VS'+b.replace('.'+tail,'')
            #拆分Indels和SNVs
        for i in ['snps','indels']:
            print('bcftools view -v {0} {1} -O z -o {2}'.format(
                i,input_dir+'/'+a,
                out_dir+'/'+a.replace(tail,i+'.'+tail)))
            os.system('bcftools view -v {0} {1} -O z -o {2}'.format(
                            i,input_dir+'/'+a,
                            out_dir+'/'+a.replace(tail,i+'.'+tail)))
            os.system('bcftools view -v {0} {1} -O z -o {2}'.format(
                            i,input_dir+'/'+b,
                            out_dir+'/'+b.replace(tail,i+'.'+tail)))

            #计算Jaccard Index
            os.system('bedtools jaccard -a {0} -b {1} > {2}'.format(out_dir+'/'+a.replace(tail,i+'.'+tail),out_dir+'/'+b.replace(tail,i+'.'+tail),out_dir+'/'+name+'_'+i+'.txt'))
            print('bedtools jaccard -a {0} -b {1} > {2}'.format(out_dir+'/'+a.replace(tail,i+'.'+tail),out_dir+'/'+b.replace(tail,i+'.'+tail),out_dir+'/'+name+'_'+i+'.txt'))
            print(out_dir+'/'+name+'_'+i+'.txt')
            JI_result = pd.read_csv(out_dir+'/'+name+'_'+i+'.txt',sep='\t')
            print(JI_result.shape[0])
            if (JI_result.shape[0]!=0):
                JI_dict['type'].append(i)
                JI_dict['JI'].append(JI_result['jaccard'].values[0])
                JI_dict['source'].append(name)
        

JI_df = pd.DataFrame.from_dict(JI_dict)
JI_df.to_csv(out_dir+'/'+'All.Jacard_index_bedtools.csv')
#os.system('rm {}/*.txt'.format(out_dir))
#os.system('rm {}/*.vcf.gz'.format(out_dir))

