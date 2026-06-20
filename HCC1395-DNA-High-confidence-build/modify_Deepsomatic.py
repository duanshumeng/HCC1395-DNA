## -*- coding: UTF-8 -*-
#source activate mapping

# Post-process GATK4's MuTect2 output. The main purpose is to split multi-allelic records into one variant record per line.

import sys, os, argparse, gzip
import regex as re

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Variant Call Type, i.e., snp or indel

class VcfLine:
    def __init__(self, line):
        parts = line.split('\t')
        self.chrom = parts[0]
        self.pos = parts[1]
        self.id = parts[2]
        self.ref = parts[3]
        self.alt = parts[4]
        self.qual = parts[5]
        self.filter = parts[6]
        self.info = parts[7]
        self.format = parts[8].split(':')
        self.sample = parts[9].split(':')

    def get_sample_dict(self):
        return dict(zip(self.format, self.sample))

def split_deepsomatic_vcf(infile, snv_out_path, indel_out_path):
    with genome.open_textfile(infile) as vcf_in, open(snv_out_path, 'w') as snv_out, open(indel_out_path, 'w') as indel_out:
        for line in vcf_in:
            if line.startswith('#'):
                snv_out.write(line)
                indel_out.write(line)
                continue

            v = VcfLine(line.rstrip())
            alts = v.alt.split(',')
            
            # 如果是单等位基因，直接分类输出
            if len(alts) == 1:
                target = snv_out if len(v.ref) == 1 and len(v.alt) == 1 else indel_out
                target.write(line)
                continue

            # 处理多等位基因
            s_dict = v.get_sample_dict()
            ad = s_dict['AD'].split(',')
            vaf = s_dict['VAF'].split(',')
            pl = s_dict['PL'].split(',')

            for i, alt_base in enumerate(alts):
                # 1. 拆分 AD: 保留 Ref 和 当前的 Alt
                new_ad = f"{ad[0]},{ad[i+1]}"
                
                # 2. 拆分 VAF
                new_vaf = vaf[i]

                # 3. 拆分 PL 并根据评分重新赋值 GT
                # 原始 PL 顺序 (二倍体, 2个ALT): 0/0, 0/1, 1/1, 0/2, 1/2, 2/2
                # 索引: 0=0/0, 1=0/1, 2=1/1, 3=0/2, 4=1/2, 5=2/2
                # 拆分给第一个 ALT (i=0): 0/0, 0/1, 1/1 -> 索引 0, 1, 2
                # 拆分给第二个 ALT (i=1): 0/0, 0/2, 2/2 -> 索引 0, 3, 5
                if i == 0:
                    sub_pl = [int(pl[0]), int(pl[1]), int(pl[2])]
                else:
                    sub_pl = [int(pl[0]), int(pl[3]), int(pl[5])]
                
                # 归一化 PL (最小值减为0)
                min_pl = min(sub_pl)
                norm_pl = [p - min_pl for p in sub_pl]
                new_pl_str = ",".join(map(str, norm_pl))

                # 根据 PL 最小值确定 GT
                best_gt_idx = norm_pl.index(0)
                if best_gt_idx == 0: new_gt = "0/0"
                elif best_gt_idx == 1: new_gt = "0/1"
                else: new_gt = "1/1"

                # 4. 组装新的 Sample 列
                # 注意：这里保留了原始的 GQ 和 DP
                new_sample = f"{new_gt}:{s_dict['GQ']}:{s_dict['DP']}:{new_ad}:{new_vaf}:{new_pl_str}"
                
                new_line = "\t".join([
                    v.chrom, v.pos, v.id, v.ref, alt_base, 
                    v.qual, v.filter, v.info, ":".join(v.format), new_sample
                ]) + "\n"

                # 5. 分流输出
                if len(v.ref) == 1 and len(alt_base) == 1:
                    snv_out.write(new_line)
                else:
                    indel_out.write(new_line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-infile',  '--input_vcf', type=str, help='Input Deepsomatic VCF file', required=True)
    parser.add_argument('-snv',     '--snv_out',   type=str, help='Output VCF file', required=True)
    parser.add_argument('-indel',   '--indel_out', type=str, help='Output VCF file', required=True)
    args = parser.parse_args()
    split_deepsomatic_vcf(args.input_vcf, args.snv_out, args.indel_out)