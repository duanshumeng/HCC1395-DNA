import pandas as pd
import subprocess
import argparse


parser = argparse.ArgumentParser(description='Check and clean VAF')
parser.add_argument('-i', '--vcf_input', help='Path to the VCF files, includes *.vcf')
parser.add_argument('-o','--out_file',help='*output.txt')

args = parser.parse_args()

vcf_input = args.vcf_input
out_file = args.out_file

##清理VCF文件

def uniq_lst(lst):
    lst_new = []
    for i in lst:
        if i not in lst_new:
            lst_new.append(i)
    return ';'.join(lst_new)


def get_credibility(info):
    lst_info = info.split(';')
    if len(lst_info) == 8:
        count = int(lst_info[6].split('=')[1].split(',')[0])+int(lst_info[7].split('=')[1].split(',')[0])
    elif len(lst_info) == 9:
        count = int(lst_info[6].split('=')[1].split(',')[0])+int(lst_info[7].split('=')[1].split(',')[0])+int(lst_info[8].split('=')[1].split(',')[0])
    else:
        print(len(lst_info))
        print(info)
    return (count/39)*100
        
def get_vcf_info(vcf_input):
    #取出VCF标题所在的行的位置
    A = subprocess.check_output("cat %s | grep -n '#CHROM' | awk -F ':' '{print $1}'" % (vcf_input),shell=True)

    vcf_header = subprocess.check_output("cat %s | grep -n '##' " % (vcf_input),shell=True)

    #跳过前几行读取vcf文件
    vcf_file = pd.read_table(vcf_input,skiprows=int(A)-1)

    print(vcf_file.shape)
    #去重
    vcf_file_new = vcf_file.drop_duplicates()
    print(vcf_file_new.shape)

    #去除VAF=0的位点
    vcf_file_new1 = vcf_file_new[~vcf_file_new['INFO'].str.contains('VAF=0;')]

    #去除INFO中的重复信息
    vcf_file_new1['INFO'] = vcf_file_new1['INFO'].apply(lambda x:uniq_lst(x.split(';')))

    #检查vcf_file_new1['INFO']信息
    vcf_file_new1['INFO'].apply(lambda x:len(x.split(';'))).unique()
    
    

    vcf_header = vcf_header.decode('utf-8')
    lines = vcf_header.split('\n')
    with open(out_file, 'w') as w:
            for i, line in enumerate(lines, 1):
                        line_split = line.split(':',maxsplit=1)
                        if len(line_split) > 1:
                            w.write(line_split[1]+'\n')
    vcf_file_new1['Credibility (%)'] = vcf_file_new1['INFO'].apply(lambda x:get_credibility(x))
    vcf_file_new1.to_csv(out_file,mode='a', index=False,sep='\t')
    #提取每个位点被检测到的次数：TNseq_n + strelka_n + TNscope_n

    vcf_file_new1.to_csv(out_file.replace('vcf','txt'),index=False,sep='\t')


get_vcf_info(vcf_input)