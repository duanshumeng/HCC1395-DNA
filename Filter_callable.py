import os
import pandas as pd
import argparse
### input arguments
parser = argparse.ArgumentParser(description="This script is to build HCC1395 High Confidence Region")

parser.add_argument('-d', '--FileDirectoryPath', type=str, help='Path of input bed files',  required=True)
parser.add_argument('-p', '--Prefix', type=str, help='Prefix of three repeat sampleIDs,eg."WGS_WUX_2023,BGI_T10_WGS_2023,BGI_T20_WGS_2023"',  required=True)
parser.add_argument('-o','--out_prefix',type=str,help='Output prefix of direct',required=True)

args = parser.parse_args()
input_dir = args.FileDirectoryPath
out_prefix = args.out_prefix
prefix = args.Prefix

out_dir = input_dir+'/'+out_prefix
os.makedirs(out_dir,exist_ok=True)

bed_files = [i for i in os.listdir() if i.endswith('_callable.bed')]
smps = list(set([i.split('_T_')[0] for i in bed_files if '_T_' in i]))

#Step1: Intersection of tumor and normal
for i in smps:
    sub_bed = [a for a in bed_files if a.startswith(i)]
    print(sub_bed)
    os.system('bedtools multiinter -i {0}/{1} {0}/{2} > {3}/{4}'.format(
        input_dir,
        sub_bed[0],
        sub_bed[1],
        out_dir,
        i+'.tmp.callable.bed'))

    call_bed = pd.read_csv(out_dir+'/'+i+'.tmp.callable.bed',header=None,sep='\t')
    call_bed[call_bed[3]>1][[0,1,2]].to_csv(out_dir+'/'+i+'.callable.bed',sep='\t',header=None,index=None)
    os.system('rm {0}/{1}'.format(out_dir,i+'.tmp.callable.bed'))

#Step2: 2/3 of 3 repeats
out_dir_2 = input_dir+'/'+out_prefix+'/'+'Filter_by_repeat'
os.makedirs(out_dir_2,exist_ok=True)

rep_lst = prefix.split(',')
print(rep_lst)
print(out_dir)
for s in rep_lst:
    bed_sub = [i for i in os.listdir(out_dir) if (i.endswith('.callable.bed') & i.startswith(s))]
    print(bed_sub)
    os.system('bedtools multiinter -i {0}/{1} {0}/{2} {0}/{3} > {4}/{5}'.format(
        out_dir,
        bed_sub[0],
        bed_sub[1],
        bed_sub[2],
        out_dir_2,
        s+'.tmp.callable.bed'))
    call_bed = pd.read_csv(out_dir_2+'/'+s+'.tmp.callable.bed',header=None,sep='\t')
    call_bed[call_bed[3]>1][[0,1,2]].to_csv(out_dir_2+'/'+s+'.callable.bed',sep='\t',header=None,index=None)
    os.system('rm {0}/{1}'.format(out_dir_2,s+'.tmp.callable.bed'))

#no repeat files
for i in os.listdir(out_dir):
    if i.endswith('.bed'):
        s= i.split('.')[0].rsplit('_',1)[0]
        out = s+'.callable.bed'
        if s not in rep_lst:
            print(i)
            os.system('cp {0}/{1} {2}/{3}'.format(
               out_dir,
               i, 
               out_dir_2,
               out
            ))


    
