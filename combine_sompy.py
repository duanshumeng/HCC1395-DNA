#!/opt/local/cobweb/envs/hrd/bin/python3.9
import pandas as pd
import os
import itertools
import sys
import argparse

parser = argparse.ArgumentParser(description='Compare with SEQC2 reference datasets and calculate F1 score')
parser.add_argument('-i', '--input_dir', help='Path to the input files, includes *.vcf.gz.')
parser.add_argument('-o', '--out_dir', help='Path to the output dir')
parser.add_argument('-s', '--seqc2', help='如果是和seqc2的高置信数据集比较，则设置为yes，否则为no')
args = parser.parse_args()

input_dir = args.input_dir
#year = sys.argv[2]
#caller = sys.argv[3]
seqc2 = args.seqc2
out_dir = args.out_dir
#如果是和seqc2的高置信数据集比较，则设置为yes，否则为no
#command demo

 

def unique_combinations(arr):
    combos = itertools.combinations(arr, 2)
    unique_combos = set(combos)
    return unique_combos

aim_arr = []

for i in os.listdir(input_dir):
    #if year in i and caller in i and i.endswith('.gz'):
    if i.endswith('.gz'):
        aim_arr.append(i)

print(aim_arr)
#out_dir = caller+'_'+year
os.makedirs(out_dir,exist_ok=True)

if seqc2 == 'yes':
    for a in aim_arr:
        #仅在高置信区间内比较
        name = a.replace('.vcf.gz','')+'VS_SEQC2'
        os.system('docker run --rm -v /data1/HCC1395/wux_WGS/somatic_SNVs/:/data hap.py:v1 \
            /opt/hap.py/bin/som.py /data/SEQC2_high_confidence/high-confidence_in_HC_regions_v1.2.1_sort_SNV_Indels.vcf.gz \
                /data/{0} -r /data/GRCh38.d1.vd1.fa -R /data/High-Confidence_Regions_v1.2.bed -o /data/{1}'.
                  format(input_dir+'/'+a,
                         out_dir+'/'+name+'_HCR'))

else:
    for out in unique_combinations(aim_arr):
        a,b = out[0],out[1]
        name = a.replace('.vcf.gz','')+'VS'+b.replace('.vcf.gz','')
        #全基因组范围内比较
        print('-------Processing '+name+'------')
        os.system('docker run --rm -v {3}:/data hap.py:v1 /opt/hap.py/bin/som.py /data/{0} /data/{1} -r /data/GRCh38.d1.vd1.fa -o /data/{2}'.format(
            input_dir+'/'+a,
            input_dir+'/'+b,
            out_dir+'/'+name,
            os.path.dirname(os.path.abspath(input_dir))))

        #仅在高置信区间内比较
        os.system('docker run --rm -v {3}:/data hap.py:v1 /opt/hap.py/bin/som.py /data/{0} /data/{1} -r /data/GRCh38.d1.vd1.fa -R /data/High-Confidence_Regions_v1.2.bed -o /data/{2}'.format(
            input_dir+'/'+a,
            input_dir+'/'+b,
            out_dir+'/'+name+'_HCR',
            os.path.dirname(os.path.abspath(input_dir))))


file_list=[]
for i in os.listdir(out_dir):
    if i.endswith('.stats.csv'):
        stats = pd.read_csv(out_dir+'/'+i,sep=',')
        name = i.split('.stats')[0]
        stats_sub = stats[['type','total.truth','total.query','tp','fp','fn','recall','precision']]
        stats_sub['source']=name
        stats_sub['F1 score'] = stats_sub.apply(lambda df:2 * (df['precision'] * df['recall']) / (df['precision'] + df['recall']) if (df['precision'] + df['recall']) != 0 else 0,axis=1)
        file_list.append(stats_sub)

stats_df = pd.concat(file_list)
stats_df.to_csv(out_dir+'/'+'F1score_stats_total.csv',index=None)
os.system('rm {}/*.txt'.format(out_dir))
os.system('rm {}/*.vcf'.format(out_dir))
os.system('rm {}/*.stats'.format(out_dir))