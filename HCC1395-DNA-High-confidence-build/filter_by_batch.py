import os
import pandas as pd
import argparse
### input arguments
parser = argparse.ArgumentParser(description="This script is to build HCC1395 High Confidence Region")

parser.add_argument('-d', '--FileDirectoryPath', type=str, help='Path of input VCF files',  required=True)
parser.add_argument('-o','--out_prefix',type=str,help='Output prefix of direct',required=True)

args = parser.parse_args()
input_dir = args.FileDirectoryPath
out_prefix = args.out_prefix

out_dir = input_dir+'/'+out_prefix
os.makedirs(out_dir,exist_ok=True)

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

vcf_sub = []
for i in os.listdir(input_dir):
    if i.endswith('.gz'):
        s= i.split('.')[0].rsplit('_',1)[0]
        vcf_sub.append(i)
print(vcf_sub)
name = 'High_confidence'
os.system('docker run --rm -v {0}:/data hap.py:v1 bcftools isec \
    /data/{1} /data/{2} /data/{3} /data/{4} /data/{5} /data/{6} /data/{7} \
    -n +4 -w1 -O z -o /data/{8} && \
            docker run --rm -v {0}:/data hap.py:v1 bcftools stats -s - /data/{8} > {9} && \
            grep "^SN" {9} | cut -f 2- > {10} '.format(
                        os.path.abspath(input_dir),
                        vcf_sub[0],
                        vcf_sub[1],
                        vcf_sub[2],
                        vcf_sub[3],
                        vcf_sub[4],
                        vcf_sub[5],
                        vcf_sub[6],
                        out_prefix+'/'+name+'.vcf.gz',
                        out_dir+'/'+name+'.stats',
                        out_dir+'/'+name+'.SNVs.txt'))
os.system('docker run --rm -v {0}:/data hap.py:v1 bcftools index /data/{1}'.format(
        os.path.abspath(out_dir),
        name+'.vcf.gz'
    ))
print('docker run --rm -v {0}:/data hap.py:v1 bcftools isec \
    /data/{1} /data/{2} /data/{3} /data/{4} /data/{5} /data/{6} /data/{7} \
    -n +2 -w1 -O z -o /data/{8} && \
            docker run --rm -v {0}:/data hap.py:v1 bcftools stats -s - /data/{8} > {9} && \
            grep "^SN" {9} | cut -f 2- > {10} '.format(
                        os.path.abspath(input_dir),
                        vcf_sub[0],
                        vcf_sub[1],
                        vcf_sub[2],
                        vcf_sub[3],
                        vcf_sub[4],
                        vcf_sub[5],
                        vcf_sub[6],
                        out_prefix+'/'+name+'.vcf.gz',
                        out_dir+'/'+name+'.stats',
                        out_dir+'/'+name+'.SNVs.txt'))


vcf_dict_out = vcf_stats(out_dir+'/'+name+'.SNVs.txt',name)
vcf_df = pd.DataFrame.from_dict(vcf_dict_out)
vcf_df.to_csv(out_dir+'/'+out_prefix+'.vcf_stats.csv')