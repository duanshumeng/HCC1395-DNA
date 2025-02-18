import subprocess
import pandas as pd
import argparse
parser = argparse.ArgumentParser(description='Ffilter VAF')
parser.add_argument('-i','--vcf_input', help='Input VCF output by TNseq')
parser.add_argument('-o','--output',help='The output file csv')
parser.add_argument('-n','--num',type=int,default=9,help='Level1 cutoff')

args = parser.parse_args()
vcf_input = args.vcf_input
output = args.output

A = subprocess.check_output("cat %s | grep -n '#CHROM' | awk -F ':' '{print $1}'" % (vcf_input),shell=True)
vcf_header = subprocess.check_output("cat %s | grep -n '##' " % (vcf_input),shell=True)
vcf_file = pd.read_table(vcf_input,skiprows=int(A)-1)
column_mapping = {vcf_file.columns[9]:'Normal',vcf_file.columns[10]:'Tumor'}
vcf_file.rename(columns=column_mapping, inplace=True)

#取出BQ
vcf_file['BQ.Numor'] = vcf_file['INFO'].apply(lambda x:x.split(';')[5].split('=')[1].split(',')[0])
vcf_file['BQ.Tumor'] = vcf_file['INFO'].apply(lambda x:x.split(';')[5].split('=')[1].split(',')[1])

#取出MQ
vcf_file['MQ.Numor'] = vcf_file['INFO'].apply(lambda x:x.split(';')[7].split('=')[1].split(',')[0])
vcf_file['MQ.Tumor'] = vcf_file['INFO'].apply(lambda x:x.split(';')[7].split('=')[1].split(',')[1])

#取出支持变异位点的reads depth
vcf_file['AD.alt.Numor'] = vcf_file['Normal'].apply(lambda x:x.split(':')[1].split(',')[1])
vcf_file['AD.alt.Tumor'] = vcf_file['Tumor'].apply(lambda x:x.split(':')[1].split(',')[1])

#取出突变频率VAF
vcf_file['VAF.Normal'] = vcf_file['Normal'].apply(lambda x:x.split(':')[2].split(',')[0])
vcf_file['VAF.Tumor'] = vcf_file['Tumor'].apply(lambda x:x.split(':')[2].split(',')[0])

#取出total depth
vcf_file['DP.Numor'] = vcf_file['Normal'].apply(lambda x:x.split(':')[3])
vcf_file['DP.Tumor'] = vcf_file['Tumor'].apply(lambda x:x.split(':')[3])

#位点筛选
#1. MQ >= 40 (Tumor)
vcf_file['MQ.gteq40'] =  vcf_file['MQ.Tumor'].apply(lambda x:1 if int(x)>=40 else 0)
#2. BQ >= 30 (Tumor)
vcf_file['BQ.gteq30'] =  vcf_file['BQ.Tumor'].apply(lambda x:1 if int(x)>=30 else 0)
#3. Tumor variant depth >2
vcf_file["AD.alt.Tumor.gt2"] = vcf_file['AD.alt.Tumor'].apply(lambda x:1 if int(x)>2 else 0)
#4. tumor VAF > 10 times the normal VAF
vcf_file['TVAF.gt.10NVAF'] = vcf_file.apply(
    lambda x:1 if float(x['VAF.Tumor'])/float(x['VAF.Normal']) >10 else 0,axis=1
)
#5. Normal depth <= 20
vcf_file['DP.Numor.lt20'] = vcf_file['DP.Numor'].apply(lambda x:1 if int(x)<=20 else 0)

#6. Normal VAF <= 0.1
vcf_file['VAF.Normal.lt0.1'] = vcf_file['VAF.Normal'].apply(lambda x:1 if float(x)<=0.1 else 0)
vcf_file['Validated'] = vcf_file[['MQ.gteq40','BQ.gteq30', 'AD.alt.Tumor.gt2', 'TVAF.gt.10NVAF', 'DP.Numor.lt20','VAF.Normal.lt0.1']].apply(lambda x:'yes' if x.sum()==6 else 'no',axis=1)
vcf_file.to_csv(output,sep=',',index=None)

