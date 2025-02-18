#conda环境 hcc1395
import pandas as pd
import subprocess
import argparse
import os
from functools import reduce

parser = argparse.ArgumentParser(description='FIlter RNA editing')
parser.add_argument('-i','--vcf_input', help='Input VCFs, *.vcf or *.vcf.gz')
parser.add_argument('-e','--edit_ed', help='Input RNA editing database,such as TABLE1_hg38.txt')
#parser.add_argument('-o','--output_dir',help='The path of output dir')

args = parser.parse_args()

vcf_input = args.vcf_input
edit_ed = args.edit_ed
#output_dir = args.output_dir


ed_df = pd.read_csv(edit_ed,sep='\t')

ed_df['tag']=ed_df.apply(lambda x:'&'.join([str(x['Region']),str(x['Position']),str(x['Ref']),str(x['Ed'])]),axis=1)
ed_df['type'] = 'RNA editing'


#获取fp.vcf文件：
print('Get fp.vcf.gz.....')
if 'D5' in vcf_input:
    hc = "/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Quartet_high_confidence/LCL5.high.confidence.calls.vcf.gz"
elif 'D6' in vcf_input:
    hc = "/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Quartet_high_confidence/LCL6.high.confidence.calls.vcf.gz"
elif 'F7' in vcf_input:
     hc = "/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Quartet_high_confidence/LCL7.high.confidence.calls.vcf.gz"
elif 'M8' in vcf_input:
     hc = "/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Quartet_high_confidence/LCL8.high.confidence.calls.vcf.gz"
        
sdf = "/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/Reference/GRCh38.sdf"
bed = "/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Quartet_high_confidence/Quartet_benchmark_regions_CDS_snv_indel_202306.bed"


rtg_out = os.path.basename(vcf_input).replace('.vcf.gz','').replace('.vcf','')
if os.path.exists(rtg_out):
    os.rmdir(rtg_out)
    
os.system('rtg vcfeval -b {0} -c {1} -o {2} -t {3} --bed-regions {4}'.format(
            hc,
            vcf_input,
            rtg_out,
            sdf,
            bed))

print('rtg finished....')

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

    pattern = r'(\.\/)?\.:+'
    aim_header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    for i in vcf_file.columns:
        if i.endswith('.T'):
            aim_header.append(i)
    #仅保留1-22号染色体的信息
    chrs = ['chr'+str(i) for i in range(1,23)] + ['chrX','chrY']
    vcf_file_filter = vcf_file[vcf_file['#CHROM'].isin(chrs)]
    vcf_sub = vcf_file_filter[aim_header]
    vcf_sub['tag']=vcf_sub.apply(lambda x:'&'.join([str(x['#CHROM']),str(x['POS']),str(x['REF']),str(x['ALT'])]),axis=1)    

    return vcf_sub


rna_df = get_vcf_info(rtg_out+'/fp.vcf.gz')

rna_df_edit = pd.merge(rna_df,ed_df[['tag','type']],on='tag',how='left')
#rna_df_edit.to_csv('rna_df_edit.csv')

tp = int(subprocess.check_output("cat %s | grep None | awk -F ' ' '{print $3}'" % (rtg_out+'/summary.txt'),shell=True))


sum_dict = {'RNA editing':len(rna_df_edit['type'].dropna()),
         'Potential artifacts':rna_df_edit.shape[0]-len(rna_df_edit['type'].dropna()),
         'True positive':tp}

sum_df = pd.DataFrame.from_dict(sum_dict,orient='index')

sum_df.to_csv(rtg_out+'.summary.csv')


#合并summary.csv
df_lst = []
for i in os.listdir():
    if i.endswith('.summary.csv'):
        df = pd.read_csv(input_dir+'/'+i)
        df.columns = [['Variant type',i.split('.')[0]]]
        df_lst.append(df)


all_df = reduce(lambda x, y: pd.merge(x, y), df_lst)
all_df.to_csv('RNA_variants.csv',index=None)
