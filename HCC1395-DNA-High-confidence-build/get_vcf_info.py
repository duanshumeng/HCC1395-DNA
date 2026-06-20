#!/usr/bin/env python3
import pandas as pd
import numpy as np
from scipy import stats
import os
import subprocess
import argparse
import re
import warnings
import multiprocessing
from glob import glob

# 忽略计算过程中的 RuntimeWarning
warnings.filterwarnings("ignore", category=RuntimeWarning)

def get_args():
    parser = argparse.ArgumentParser(description='Filter and Merge VCF VAF')
    parser.add_argument('-i','--vcf_input', required=True, help='Input VCFs, separated by comma')
    parser.add_argument('-t','--vcf_type', help='SEQC2 or PGx')
    parser.add_argument('-c','--caller_num', type=int, required=True, help='Number of callers')
    parser.add_argument('-v','--variant_type', required=True, help='SNV or Indel')
    parser.add_argument('-n','--num', type=int, default=9, help='Level1 cutoff score')
    parser.add_argument('-o','--output_dir', required=True, help='Output directory')
    parser.add_argument('-a','--all_count', type=int, default=13, help='Total sample count')
    return parser.parse_args()

all_count=99
def safe_int(val):
    """
    非递归的标准安全整型转换函数。
    防止空值、单点号 '.'、或带有不可见字符的占位符强转 int 时导致流程中断。
    """
    if val is None:
        return None
    
    # 如果是 subprocess 返回的 bytes 类型，先 decode 成字符串
    if isinstance(val, bytes):
        try:
            val = val.decode('utf-8')
        except:
            pass

    clean_val = str(val).strip()
    
    # 过滤掉常见的生信缺失占位符
    if clean_val == '.' or clean_val == '' or clean_val.lower() == 'none':
        return None
        
    try:
        return int(clean_val)
    except ValueError:
        return None


def calc_95ci(data, decimals=3):
    """通用统计函数：计算均值与95%置信区间"""
    if not data or len(data) == 0: 
        return ".", (".", ".")
    mu = np.mean(data)
    sd = np.std(data, ddof=1)  # 使用样本标准差（ddof=1）
    if sd > 0 and len(data) >= 2:
        ci = stats.norm.interval(0.95, loc=mu, scale=sd / np.sqrt(len(data)))
        ci_clean = (max(0, round(ci[0], decimals)), round(ci[1], decimals))
        return round(mu, decimals), ci_clean
    return round(mu, decimals), (round(mu, decimals), round(mu, decimals))


def format_ci(val, decimals=0):
    """处理空值的字符串拼接逻辑"""
    if val == ".":
        return "."
    if isinstance(val, tuple):
        return f"{val[0]},{val[1]}"
    return f"{val:.{decimals}f}" if decimals > 0 else f"{val}"


def normalize_sample_name(name):
    if name in ["Total sum", "Confidence"]:
        return name
    # 规则1：_T → .T
    if name.endswith("_T"):
        return name[:-2] + ".T"
    # 规则2：已经是 .T
    if name.endswith(".T"):
        return name
    # 规则3：没有后缀
    return name + ".T"


def fix_header(header):
    cols = header
    print(cols)
    fixed_cols = cols[:9]
    sample_cols = cols[9:]
    new_samples = [normalize_sample_name(x) for x in sample_cols]
    return fixed_cols + new_samples


def get_vcf_info(vcf_path):
    """读取VCF并提取核心列及样本列，压缩内存"""
    if vcf_path.endswith('.gz'):
        cmd = f"zcat {vcf_path} | grep -n '#CHROM' | head -n 1 | cut -d: -f1"
        compression = 'gzip'
    else:
        cmd = f"grep -n '#CHROM' {vcf_path} | head -n 1 | cut -d: -f1"
        compression = None

    # 配合 safe_int 健壮解析系统返回的行号字节流
    line_num = safe_int(subprocess.check_output(cmd, shell=True))
    if line_num is None:
        raise ValueError(f"无法从文件找到 #CHROM 表头行，请检查文件格式或命令: {cmd}")
    
    df = pd.read_table(
        vcf_path, skiprows=line_num-1, compression=compression, low_memory=True)
    caller_name = os.path.basename(vcf_path).split('.')[0]

    valid_chrs = [f'chr{i}' for i in range(1,23)] + ['chrX','chrY']
    df = df[df['#CHROM'].isin(valid_chrs)]
    aim_header = [i for i in df.columns]
    df.columns = fix_header(aim_header)
    
    sample_cols = [c for c in df.columns if c.endswith('.T')]
    # 修复点1：移除不存在的列（Total sum/Confidence），避免样本列丢失
    keep_cols = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + sample_cols
    df_sub = df[keep_cols]

    # 压缩 dtype
    for col in ['#CHROM','REF','ALT']:
        df_sub[col] = df_sub[col].astype('category')
    df_sub['caller'] = caller_name
    df_sub['tag'] = df_sub['#CHROM'].astype(str) + "&" + df_sub['POS'].astype(str) + "&" + df_sub['REF'].astype(str) + "&" + df_sub['ALT'].astype(str)

    # 修复点2：补充Score和Confidence列（后续逻辑需要）
    df_sub['Score'] = 0
    df_sub['Confidence'] = 'Level1'  # 默认值，后续会重新计算
    print(df_sub.head())  # 防止内存爆炸
    return df_sub


def combine_info_logic(df_site, tag, smps):
    """
    终极修正版核心整合逻辑（已整合代码1的突变联合解析与边界防御体系）
    """
    try:
        #drop_duplicates
        df_site = df_site.drop_duplicates()
        res_row = df_site[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL']].head(1).copy()
        ref_base = str(res_row['REF'].iloc[0]).strip()
        
        # 核心修复：支持多等位基因联合解析，将其切分为列表 (例如: ['G', 'A'])
        raw_alt = df_site['ALT'].iloc[0]
        alt_bases = [b.strip() for b in str(raw_alt).split(',')]

        # 1. 计算全局 Mapping Quality (MQ)
        mq_vals, mq0_vals = [], []
        for info_str in df_site['INFO'].dropna():
            m_mq = re.search(r'MQ=([\d\.]+)', str(info_str))
            m_mq0 = re.search(r'MQ0=(\d+)', str(info_str))
            if m_mq: mq_vals.append(float(m_mq.group(1)))
            if m_mq0: mq0_vals.append(float(m_mq0.group(1)))
        
        MQ = round(np.mean(mq_vals), 2) if mq_vals else "."
        MQ0 = int(np.mean(mq0_vals)) if mq0_vals else "."

        # 2. 样本层级流式迭代
        all_dps, all_vafs = [], []
        caller_tags = []
        support_num = 0
        for caller_name, df_caller in df_site.groupby('caller'):
            caller_valid_sample_count = 0
            for _, row in df_caller.iterrows():
                fmt = str(row['FORMAT']).split(':')
                
                for s in smps:
                    sample_val = str(row[s]).strip()
                    if sample_val == '.' or ':' not in sample_val or sample_val.startswith('./.'):
                        continue
                    
                    parts = sample_val.split(':')
                    is_sample_informative = False
                    current_dp = None
                    current_vaf = None

                    # A. 提取测序深度 DP
                    if 'DP' in fmt:
                        dp_idx = fmt.index('DP')
                        if len(parts) > dp_idx:
                            current_dp = safe_int(parts[dp_idx])
                            if current_dp is not None:
                                all_dps.append(current_dp)

                    # B. 动态提取/计算 VAF
                    try:
                        if 'VAF' in fmt:
                            vaf_idx = fmt.index('VAF')
                            current_vaf = float(parts[vaf_idx].split(',')[0])
                        elif 'AF' in fmt:
                            af_idx = fmt.index('AF')
                            current_vaf = float(parts[af_idx].split(',')[0])
                        
                        # C. 兼容 Strelka2 包含 AU/CU/GU/TU 字段的多等位基因 SNV 情况
                        elif any(f'{b}U' in fmt for b in alt_bases) and current_dp and current_dp > 0:
                            total_alt_tier1_reads = 0
                            has_valid_reads = False
                            for b in alt_bases:
                                field_name = f'{b}U'
                                if field_name in fmt:
                                    b_idx = fmt.index(field_name)
                                    reads = safe_int(parts[b_idx].split(',')[0])
                                    if reads is not None:
                                        total_alt_tier1_reads += reads
                                        has_valid_reads = True
                            
                            if has_valid_reads:
                                current_vaf = total_alt_tier1_reads / current_dp
                            else:
                                current_vaf = None
                            
                        # D. 兼容 Strelka2 标准 InDel 模式 (TAR/TIR) —— 完美边界防御版
                        elif 'TIR' in fmt and 'TAR' in fmt:
                            tir_idx = fmt.index('TIR')
                            tar_idx = fmt.index('TAR')
                            
                            tir_reads = safe_int(parts[tir_idx].split(',')[0])
                            tar_reads = safe_int(parts[tar_idx].split(',')[0])
                            
                            if tir_reads is not None and tar_reads is not None:
                                total_reads = tir_reads + tar_reads
                                if total_reads > 0:
                                    current_vaf = tir_reads / total_reads
                                else:
                                    current_vaf = None
                            else:
                                current_vaf = None
                                
                    except Exception as e:
                        print(f"[-] 警告: 位点 {tag} 样本 {s} 解析失败. 原因: {e}")
                        current_vaf = None

                    if current_vaf is not None:
                        all_vafs.append(current_vaf)
                        is_sample_informative = True

                    if is_sample_informative or current_dp is not None:
                        caller_valid_sample_count += 1
            
            caller_tags.append(f"{caller_name}_n={caller_valid_sample_count}")
            support_num += caller_valid_sample_count
            print("support_num============",support_num)

        # 3. 计算均值和置信区间
        dp_avg, dp_ci = calc_95ci(all_dps, 0)
        vaf_avg, vaf_ci = calc_95ci(all_vafs, 3)

        # 4. 结果包装
        info_field_parts = [
            f"MQ={MQ}", f"MQ0={MQ0}",
            f"DP={dp_avg}", f"DP95={format_ci(dp_ci)}",
            f"VAF={vaf_avg}", f"VAF95={format_ci(vaf_ci, 3)}"
        ] + sorted(caller_tags)+[f"Confidence={str(round(support_num/all_count,3))}"]
        
        # 兼容原有流程中的表头映射
        res_row['FILTER'] = df_site['Confidence_new'].iloc[0] if 'Confidence_new' in df_site.columns else (df_site['Confidence'].iloc[0] if 'Confidence' in df_site.columns else df_site['FILTER'].iloc[0])
        res_row['INFO'] = ";".join(info_field_parts)
        #res_row['Confidence'] = support_num/all_count

        return res_row

    except Exception as e:
        print(f"[-] 严重致命错误 - 位点 {tag}: {e}")
        return pd.DataFrame()


def worker(chrom_df, tags, smps, out_file):
    """每个进程处理一个染色体"""
    # 修复点7：进程内写入时添加列名（仅首行）
    first_write = True
    with open(out_file,'w') as w:
        for t in tags:
            site_data = chrom_df[chrom_df['tag']==t]
            res = combine_info_logic(site_data, t, smps)
            if not res.empty:
                # 首行写入列名，后续行不写
                res.to_csv(w, sep='\t', header=first_write, index=False)
                first_write = False


def main():
    args = get_args()
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir, exist_ok=True)
    output_dir = args.output_dir
    variant_type = args.variant_type

    # 1. 读取VCF
    vcf_lst = []
    for f in args.vcf_input.split(','):
        print(f"Reading {f}...")
        vcf_lst.append(get_vcf_info(f))
    vcf_all = pd.concat(vcf_lst, ignore_index=True)
    vcf_all.to_csv(output_dir+'/'+'vcf_all.tsv',sep='\t',index=None)

    # 2. Score
    score_dict = {'Level1':3,'Level2':2,'Level3':1,'Level4':0}
    vcf_all['Score'] = vcf_all['Confidence'].replace(score_dict)
    print("=======vcf_all=========",vcf_all)

    # 3. 多Caller过滤
    caller_counts = vcf_all.groupby(['tag','caller']).size().reset_index()
    print("====caller_counts===",caller_counts)
    multi_tags = caller_counts.groupby('tag').size()
    multi_tags = multi_tags[multi_tags>1].index
    vcf_sub = vcf_all[vcf_all['tag'].isin(multi_tags)]
    vcf_sub_score = vcf_sub.groupby('tag')['Score'].sum().reset_index()
    vcf_res = pd.merge(vcf_sub, vcf_sub_score, on='tag', suffixes=('','_Sum'))

    # Confidence_new
    if args.caller_num==2:
        vcf_res['Confidence_new'] = vcf_res['Score_Sum'].apply(lambda x: "HighConf" if x >= 5 else ("MedConf" if x>=4 else "LowConf"))
    else:
        n = args.num
        vcf_res['Confidence_new'] = vcf_res['Score_Sum'].apply(lambda x:'HighConf' if x >= n else ('MedConf' if x >= n-2 else ('LowConf' if x >= n-4 else 'Unclassified')))

    tmp_file = os.path.join(output_dir,'Merge_{}.tmp.txt'.format(variant_type))
    # 修复点8：确保临时文件保留所有样本列
    vcf_res.to_csv(tmp_file, sep='\t', index=False, header=True)

    # 4. 只处理 High/Med
    vcf_to_combine = vcf_res[vcf_res['Confidence_new'].isin(['HighConf','MedConf'])]
    smps = [c for c in vcf_to_combine.columns if c.endswith('.T')]

    # 5. 多进程
    print("Starting Multi-processing combination...")
    chrs = vcf_to_combine['#CHROM'].unique()
    pool = multiprocessing.Pool(processes=min(len(chrs), multiprocessing.cpu_count()))
    tmp_files = []

    for chrom in chrs:
        chr_df = vcf_to_combine[vcf_to_combine['#CHROM']==chrom]
        chr_tags = chr_df['tag'].unique()
        out_file = os.path.join(output_dir,f"worker_{chrom}.tsv")
        tmp_files.append(out_file)
        pool.apply_async(worker, args=(chr_df, chr_tags, smps, out_file))

    pool.close()
    pool.join()

    # 6. 合并
    print("Finalizing results...")
    # 修复点9：合并时指定列名，确保df_final有列名
    df_list = []
    for f in tmp_files:
        if os.path.exists(f) and os.path.getsize(f) > 0:
            df = pd.read_table(f)
            df_list.append(df)
    df_final = pd.concat(df_list, ignore_index=True) if df_list else pd.DataFrame()
    
    # 确保列名存在（兜底）
    if not df_final.empty and 'FILTER' not in df_final.columns:
        df_final.columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
    
    print(df_final)

    # 7. 输出VCF
    vcf_header = '''##fileformat=VCFv4.2
##fileDate=20230802
##reference=GRCh38.d1.vd1
##FILTER=<ID=HighConf,Description="highly confident that it is a real somatic mutation">
##FILTER=<ID=MedConf,Description="confident that it is a real somatic mutation">
##FILTER=<ID=LowConf,Description="not very confident that it is a real somatic mutation">
##FILTER=<ID=Unclassified,Description="likely not a real somatic mutation">
##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Average depth of tumor samples">
##INFO=<ID=DP95,Number=2,Type=Integer,Description="Estimated 95% confidence interval for DP in tumor samples">
##INFO=<ID=VAF,Number=1,Type=Float,Description="Average VAF of tumor samples">
##INFO=<ID=VAF95,Number=2,Type=Float,Description="Estimated 95% confidence interval for VAF in tumor samples">
##INFO=<ID=strelka_n,Number=2,Type=Integer,Description="The number of files support the variant in strelka, total Smaples">
##INFO=<ID=Deepsomatic_n,Number=2,Type=Integer,Description="The number of files support the variant in Deepsomatic, total Samples">
##INFO=<ID=Mutect2_n,Number=2,Type=Integer,Description="The number of files support the variant in Mutect2, total Samples">
##INFO=<ID=Strelka_n,Number=2,Type=Integer,Description="The number of files support the variant in Strelka, total Samples">
##INFO=<ID=somaticSniper_n,Number=2,Type=Integer,Description="The number of files support the variant in TNscope, total Samples">
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr4,length=190214555>
##contig=<ID=chr5,length=181538259>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr12,length=133275309>
##contig=<ID=chr13,length=114364328>
##contig=<ID=chr14,length=107043718>
##contig=<ID=chr15,length=101991189>
##contig=<ID=chr16,length=90338345>
##contig=<ID=chr17,length=83257441>
##contig=<ID=chr18,length=80373285>
##contig=<ID=chr19,length=58617616>
##contig=<ID=chr20,length=64444167>
##contig=<ID=chr21,length=46709983>
##contig=<ID=chr22,length=50818468>
##contig=<ID=chrX,length=156040895>
##contig=<ID=chrY,length=57227415>
'''
    final_vcf_path = os.path.join(output_dir,'All_s{}_v1.vcf'.format(variant_type))
    with open(final_vcf_path, 'w') as f:
        f.write(vcf_header)
        if not df_final.empty:
            header_line = '\t'.join(df_final.columns) + '\n'
            f.write(header_line)
            df_final.to_csv(f, index=False, sep='\t', header=False)
    
    out_file = os.path.join(output_dir,'high-confidence_s{}_v1.vcf'.format(variant_type))
    with open(out_file,'w') as w:
        w.writelines(vcf_header)
        if not df_final.empty and 'FILTER' in df_final.columns:
            df_out = df_final[df_final['FILTER'].isin(['HighConf','MedConf'])].copy()
            df_out['FILTER'] = 'PASS;' + df_out['FILTER']
            header_line = '\t'.join(df_out.columns) + '\n'
            w.write(header_line)
            df_out.to_csv(w, mode='a', index=False, sep='\t', header=False)

if __name__=='__main__':
    main()