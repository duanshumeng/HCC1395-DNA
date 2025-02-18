import os
import pandas as pd
import argparse
### input arguments
parser = argparse.ArgumentParser(description="This script is to build HCC1395 High Confidence Region")

parser.add_argument('-d', '--FileDirectoryPath', type=str, help='Path of input VCF files',  required=True)
parser.add_argument('-t','--Type',type=str,help = 'Type of WGS source, SEQC2 or PGx')
parser.add_argument('-o','--out_prefix',type=str,help='Output prefix of direct',required=True)

args = parser.parse_args()
input_dir = args.FileDirectoryPath
out_prefix = args.out_prefix
Type = args.Type

print(Type)
out_dir = input_dir+'/'+out_prefix
os.makedirs(out_dir,exist_ok=True)

#统计vcf数量
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
if Type == 'PGx':
    smps = set([i.split('.')[0] for i in os.listdir(input_dir) if '2023' in i and i.endswith('.vcf.gz')])

elif Type == 'SEQC2':
    smps = set([i.split('.')[0] for i in os.listdir(input_dir) if i.endswith('.vcf.gz')])
else:
    print('Type must be PGx or SEQC2!!!')
print(smps)
for s in smps:
    vcf_sub = [i for i in os.listdir(input_dir) if (i.endswith('.gz') & i.startswith(s)) ]
    smp = s
    print(smp,vcf_sub)
    if Type == 'PGx':
        name = smp+'.'+'TNseq' +'_' + 'TNscope'+'_'+'strelka'
    
    elif Type == 'SEQC2':
        name = smp+'.'+'muTect2' +'_' + 'somaticSniper'+'_'+'strelka'
        
    else:
        print('Type must be PGx or SEQC2!!!')

    os.system('docker run --rm -v {0}:/data hap.py:v1 bcftools isec /data/{1} /data/{2} /data/{3} -n +2 -w1 -O z -o /data/{4} && \
            docker run --rm -v {0}:/data hap.py:v1 bcftools stats -s - /data/{4} > {5} && \
            grep "^SN" {5} | cut -f 2- > {6} '.format(
                        os.path.abspath(input_dir),
                        vcf_sub[0],
                        vcf_sub[1],
                        vcf_sub[2],
                        out_prefix+'/'+name+'.isec.vcf.gz',
                        out_dir+'/'+name+'.isec.stats',
                        out_dir+'/'+name+'.isec.SNVs.txt'))
    os.system('docker run --rm -v {0}:/data hap.py:v1 bcftools index /data/{1}'.format(
        os.path.abspath(out_dir),
        name+'.isec.vcf.gz'
    ))
    print('docker run --rm -v {0}:/data hap.py:v1 bcftools isec /data/{1} /data/{2} /data/{3} -n +2 -w1 -O z -o /data/{4} && \
            docker run --rm -v {0}:/data hap.py:v1 bcftools stats -s - /data/{4} > {5} && \
            grep "^SN" {5} | cut -f 2- > {6} '.format(
                        os.path.abspath(input_dir),
                        vcf_sub[0],
                        vcf_sub[1],
                        vcf_sub[2],
                        out_prefix+'/'+name+'.isec.vcf.gz',
                        out_dir+'/'+name+'.isec.stats',
                        out_dir+'/'+name+'.isec.SNVs.txt'))

    vcf_dict_out = vcf_stats(out_dir+'/'+name+'.isec.SNVs.txt',name)
    vcf_df = pd.DataFrame.from_dict(vcf_dict_out)
    vcf_lst.append(vcf_df)

vcf_df_all = pd.concat(vcf_lst)
vcf_df_all.to_csv(out_dir+'/'+out_prefix+'.vcf_stats.csv')
os.system('rm {}/*.txt'.format(out_dir))
os.system('rm {}/*.stats'.format(out_dir))