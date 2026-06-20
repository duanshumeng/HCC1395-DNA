import pandas as pd
import os
import itertools
import sys
import argparse

parser = argparse.ArgumentParser(description='Get intersection or union of Callers(TNseq/TNscope/strelka)')
parser.add_argument('-i', '--input_dir', help='Path to the input files, includes *.vcf.gz.')
parser.add_argument('-o', '--out_prefix', help='Prefix of output path')
parser.add_argument('-n', '--n_caller', help='需要合并的Caller类型，2:任意两个caller取交集；3:三个caller取交集；    \
    1:取TNseq和strelka的交集; 0: 取TNseq和strelka的并集')
args = parser.parse_args()

input_dir = args.input_dir
out_prefix = args.out_prefix
#n_caller=2(caller任意组合两两取交集)
#n_caller=3(三个中任意两个取交集)
#n_caller=1(取TNseq和strelka的交集)
n_caller = int(args.n_caller)

out_dir = input_dir+'/'+out_prefix
os.makedirs(out_dir,exist_ok=True)

def unique_combinations(arr):
    combos = itertools.combinations(arr, 2)
    unique_combos = set(combos)
    return unique_combos

#统计每个caller两两组合后得到的vcf数量
def vcf_stats(isec_txt,name):
    vcf_dict = {'type':[],'count':[],'source':[]}
    with open(isec_txt,'r') as f:
        for s in f.readlines():
            if 'SNPs' in s:
                    a_snv = s.split(':')[1]
                    vcf_dict['type'].append('SNVs')
                    vcf_dict['count'].append(a_snv.replace(' ',''))
                    vcf_dict['source'].append(name)
            elif 'indels' in s:
                    a_indel = s.split(':')[1]
                    vcf_dict['type'].append('Indels')
                    vcf_dict['count'].append(a_indel.replace(' ',''))
                    vcf_dict['source'].append(name)
    return vcf_dict 

vcf_lst = list()
#计算每个样本任意两个caller之间的交集
for year in list(set([i.split('.')[0] for i in os.listdir(input_dir) if i.endswith('.gz')])):
#for year in ['WGS_EA','WGS_FD','WGS_NS']:
      print(year)
      if n_caller==2:
        for caller in unique_combinations(['TNseq','TNscope','strelka']):
            c1,c2 = caller[0],caller[1]
            name = year + '_'+ caller[0] +'_' + caller[1]
            print(c1,'-----',c2)
            c1_vcf = ''
            c2_vcf = ''
            for vcf in os.listdir(input_dir):
                if vcf.endswith('.gz'):
                  if year in vcf and c1 in vcf:
                      c1_vcf = vcf
                  elif year in vcf and c2 in vcf:
                      c2_vcf = vcf
            print(c1_vcf,'-----',c2_vcf)
            os.system('bcftools isec {0} {1} -n +2 -w1 -O z -o {2} && bcftools stats -s - {2} > {3} && grep "^SN" {3} | cut -f 2- > {4} '.format(input_dir+'/'+c1_vcf,input_dir+'/'+c2_vcf,out_dir+'/'+name+'.isec.vcf.gz',out_dir+'/'+name+'.isec.stats',out_dir+'/'+name+'.isec.SNVs.txt'))
            vcf_dict_out = vcf_stats(out_dir+'/'+name+'.isec.SNVs.txt')
      
      elif n_caller==3:
            name = year + '_'+ 'TNseq' +'_' + 'TNscope'+'_'+'strelka'
            c1_vcf = ''
            c2_vcf = ''
            c3_vcf = ''
            for vcf in os.listdir(input_dir):
                if vcf.endswith('.gz'):
                     print(vcf)
                     if year in vcf and 'TNseq' in vcf:
                          c1_vcf = vcf
                     elif year in vcf and 'TNscope' in vcf:
                          c2_vcf = vcf
                     elif year in vcf and ('stralka' in vcf or 'strelka' in vcf):
                          c3_vcf = vcf
            os.system('bcftools isec {0} {1} {2} -n +2 -w1 -O z -o {3} && bcftools stats -s - {3} > {4} && grep "^SN" {4} | cut -f 2- > {5} '.format(
                 input_dir+'/'+c1_vcf,
                 input_dir+'/'+c2_vcf,
                 input_dir+'/'+c3_vcf,
                 out_dir+'/'+name+'.isec.vcf.gz',
                 out_dir+'/'+name+'.isec.stats',
                 out_dir+'/'+name+'.isec.SNVs.txt'))
            print('bcftools isec {0} {1} {2} -n +2 -w1 -O z -o {3} && bcftools stats -s - {3} > {4} && grep "^SN" {4} | cut -f 2- > {5} '.format(
                 input_dir+'/'+c1_vcf,
                 input_dir+'/'+c2_vcf,
                 input_dir+'/'+c3_vcf,
                 out_dir+'/'+name+'.isec.vcf.gz',
                 out_dir+'/'+name+'.isec.stats',
                 out_dir+'/'+name+'.isec.SNVs.txt'))
            
            vcf_dict_out = vcf_stats(out_dir+'/'+name+'.isec.SNVs.txt')

      elif n_caller ==1:
            c1,c2 = 'TNseq','strelka'
            name = year + '_'+ c1 +'_' + c2
            print(c1,'-----',c2)
            c1_vcf = ''
            c2_vcf = ''
            for vcf in os.listdir(input_dir):
                if vcf.endswith('.gz'):
                  if year in vcf and c1 in vcf:
                      c1_vcf = vcf
                  elif year in vcf and c2 in vcf:
                      c2_vcf = vcf
            print(c1_vcf,'-----',c2_vcf)
            os.system('docker run --rm -v {0}:/data hap.py:v1 bcftools isec /data/{1} /data/{2} -n +2 -w1 -O z -o /data/{3} && \
                docker run --rm -v {0}:/data hap.py:v1 bcftools stats -s - /data/{3} > {4} && \
                    grep "^SN" {4} | cut -f 2- > {5} '.format(
                        os.path.abspath(input_dir),
                        c1_vcf,
                        c2_vcf,
                        out_prefix+'/'+name+'.isec.vcf.gz',
                        out_dir+'/'+name+'.isec.stats',
                        out_dir+'/'+name+'.isec.SNVs.txt'))
            vcf_dict_out = vcf_stats(out_dir+'/'+name+'.isec.SNVs.txt',name)
            vcf_df = pd.DataFrame.from_dict(vcf_dict_out)
            vcf_lst.append(vcf_df)
      
      elif n_caller ==0:
            c1,c2 = 'TNseq','strelka'
            name = year + '_'+ c1 +'_' + c2
            print(c1,'-----',c2)
            c1_vcf = ''
            c2_vcf = ''
            for vcf in os.listdir(input_dir):
                if vcf.endswith('.gz'):
                  if year in vcf and c1 in vcf:
                      c1_vcf = vcf
                  elif year in vcf and c2 in vcf:
                      c2_vcf = vcf
            print(c1_vcf,'-----',c2_vcf)
            os.system('docker run --rm -v {0}:/data hap.py:v1 bcftools merge /data/{1} /data/{2} -Oz -o /data/{3} && \
                docker run --rm -v {0}:/data hap.py:v1 bcftools stats -s - /data/{3} > {4} && \
                    grep "^SN" {4} | cut -f 2- > {5} '.format(
                        os.path.abspath(input_dir),
                        c1_vcf,
                        c2_vcf,
                        out_prefix+'/'+name+'.union.vcf.gz',
                        out_dir+'/'+name+'.union.stats',
                        out_dir+'/'+name+'.union.SNVs.txt'))
            vcf_dict_out = vcf_stats(out_dir+'/'+name+'.union.SNVs.txt',name)
            vcf_df = pd.DataFrame.from_dict(vcf_dict_out)
            vcf_lst.append(vcf_df)



vcf_df_all = pd.concat(vcf_lst)
vcf_df_all.to_csv(out_dir+'/'+out_prefix+'.vcf_stats.csv')
os.system('rm {}/*.txt'.format(out_dir))
os.system('rm {}/*.stats'.format(out_dir))
          



              
        

    

