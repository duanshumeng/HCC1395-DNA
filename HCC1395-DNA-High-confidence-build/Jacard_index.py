import pandas as pd
import os
import itertools
import sys

input_dir = sys.argv[1]
year = sys.argv[2]
caller = sys.argv[3]

def unique_combinations(arr):
    combos = itertools.combinations(arr, 2)
    unique_combos = set(combos)
    return unique_combos

aim_arr = []

for i in os.listdir(input_dir):
    if year in i and caller in i and i.endswith('.gz'):
        aim_arr.append(i)

print(aim_arr)
out_dir = caller+'_'+year
os.makedirs(out_dir,exist_ok=True)

JI_dict = {'type':[],
           'JI':[],
           'source':[]}
for out in unique_combinations(aim_arr):
    a,b = out[0],out[1]
    name = a.split('.')[0]+'VS'+b.split('.')[0]+'_'+caller
    #取交集
    os.system('bcftools isec {0} {1} -n +2 -w1 -Ov -o {2}'.format(input_dir+'/'+a,
                                                                  input_dir+'/'+b,
                                                                  out_dir+'/'+name+'.isec.vcf'))
    os.system('bcftools stats -s - {0} > {1}'.format(out_dir+'/'+name+'.isec.vcf',
                                                     out_dir+'/'+name+'.isec.stats'))
    os.system('grep "^SN" {0} | cut -f 2- > {1}'.format(out_dir+'/'+name+'.isec.stats',
                                                        out_dir+'/'+name+'.isec.SNVs.txt'))
    
    #取并集
    if caller in ['TNseq','TNscope']:
          os.system('bcftools merge {0} {1} -o {2}'.format(input_dir+'/'+a,
                                                      input_dir+'/'+b,
                                                      out_dir+'/'+name+'.concat.vcf'))
    else:
          os.system('bcftools concat {0} {1} -o {2}'.format(input_dir+'/'+a,
                                                      input_dir+'/'+b,
                                                      out_dir+'/'+name+'.concat.vcf'))
    
    os.system('bcftools stats -s - {0} > {1}'.format(out_dir+'/'+name+'.concat.vcf',
                                                     out_dir+'/'+name+'.concat.stats'))
    os.system('grep "^SN" {0} | cut -f 2- > {1}'.format(out_dir+'/'+name+'.concat.stats',
                                                        out_dir+'/'+name+'.concat.SNVs.txt'))


    def get_snv_number(f):
        a_snv=int()
        a_indel=int()
        for s in f.readlines():
            if 'SNPs' in s:
                    a_snv = s.split(':')[1]
            elif 'indels' in s:
                    a_indel = s.split(':')[1]
        return(a_snv,a_indel)
    
    with open(out_dir+'/'+name+'.isec.SNVs.txt') as f1, open(out_dir+'/'+name+'.concat.SNVs.txt') as f2:
            a1_snv,a1_indel = get_snv_number(f1)
            a2_snv,a2_indel = get_snv_number(f2)
            print(name+' SNVs Jaccard index: {:.3f}'.format(int(a1_snv)/int(a2_snv)))
            print(name+' Indels Jaccard index: {:.3f}'.format(int(a1_indel)/int(a2_indel)))
            JI_dict['type'].append('SNVs')
            JI_dict['source'].append(name)
            JI_dict['JI'].append('{:.3f}'.format(int(a1_snv)/int(a2_snv)))

            JI_dict['type'].append('Indels')
            JI_dict['source'].append(name)
            JI_dict['JI'].append('{:.3f}'.format(int(a1_indel)/int(a2_indel)))


JI_df = pd.DataFrame.from_dict(JI_dict)
JI_df.to_csv(out_dir+'/'+out_dir+'.Jacard_index.csv')
os.system('rm {}/*.txt'.format(out_dir))
os.system('rm {}/*.vcf'.format(out_dir))
os.system('rm {}/*.stats'.format(out_dir))

