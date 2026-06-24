#!/opt/local/cobweb/envs/hrd/bin/python3.9
# -*- coding: UTF-8 -*-
import pandas as pd
import os
import sys

input_dir = os.path.abspath(sys.argv[1])
out_dir = os.path.abspath(sys.argv[2])
#是否只统计HCR区间的结果，yes为统计HCR区间的结果，no为统计全基因组范围内的结果
HCR = sys.argv[3]

os.makedirs(out_dir,exist_ok=True)

vcf_dict = {'type':[],
           'count':[],
           'source':[]}

if HCR == 'yes':
    for i in os.listdir(input_dir):
        if i.endswith('.gz'):
            name = i.replace('.vcf.gz','')
            #取出在高置信区间内的vcf
            os.system('bedtools intersect -header -a {0} -b {1} > {2}'.format(
                input_dir+'/'+i,
                '/data1/HCC1395/wux_WGS/somatic_SNVs/High-Confidence_Regions_v1.2.bed',
                out_dir+'/'+i.replace('vcf.gz','HCR'+'.vcf.gz')
            ))

            os.system('bcftools stats -s - {0}/{1} > {2}'.format(out_dir,i.replace('vcf.gz','HCR'+'.vcf.gz'),
                                                            out_dir+'/'+name+'.stats'))
            os.system('grep "^SN" {0} | cut -f 2- > {1}'.format(out_dir+'/'+name+'.stats',
                                                                out_dir+'/'+name+'.SNVs.txt'))
            with open(out_dir+'/'+name+'.SNVs.txt','r') as f:
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


    vcf_df = pd.DataFrame.from_dict(vcf_dict)
    vcf_df.to_csv(out_dir+'/'+os.path.basename(out_dir)+'.HCR.vcf_stats.csv',index=None)
    os.system('rm {}/*.txt'.format(out_dir))
    os.system('rm {}/*.stats'.format(out_dir))
    os.system('rm {}/*.vcf.gz'.format(out_dir))
else:
    for i in os.listdir(input_dir):
        if i.endswith('.gz') or i.endswith('.vcf'):
            name = i.replace('.vcf.gz','').replace('.vcf','')
            os.system('bcftools stats -s - {0}/{1} > {2}'.format(input_dir,i,
                                                            out_dir+'/'+name+'.stats'))
            print('bcftools stats -s - {0}/{1} > {2}'.format(input_dir,i,
                                                            out_dir+'/'+name+'.stats'))
            os.system('grep "^SN" {0} | cut -f 2- > {1}'.format(out_dir+'/'+name+'.stats',
                                                                out_dir+'/'+name+'.SNVs.txt'))
            with open(out_dir+'/'+name+'.SNVs.txt','r') as f:
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


    vcf_df = pd.DataFrame.from_dict(vcf_dict)
    vcf_df.to_csv(out_dir+'/'+os.path.basename(out_dir)+'.vcf_stats.csv',index=None)
    os.system('rm {}/*.txt'.format(out_dir))
    os.system('rm {}/*.stats'.format(out_dir))

