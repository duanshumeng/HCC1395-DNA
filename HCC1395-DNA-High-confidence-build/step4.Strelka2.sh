#!/bin/bash
sample=$1
normal_bam=$2
tumor_bam=$3
ref_bed=$4
out_dir=$5
ref_fasta="/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/Reference/GRCh38.d1.vd1.fa" # 你可以根据实际情况修改这个路径

# 创建输出目录，如果目录已存在则先删除
sample_out_dir="${out_dir}/${sample}"
if [ -d "$sample_out_dir" ]; then
echo "Output directory $sample_out_dir already exists. Removing it."
rm -rf $sample_out_dir
fi

# 创建新的输出目录
mkdir -p $sample_out_dir

# 运行Strelka2的配置脚本
configureStrelkaSomaticWorkflow.py --normalBam $normal_bam --tumorBam $tumor_bam \
--callRegions $ref_bed \
--referenceFasta $ref_fasta --runDir $sample_out_dir

nt=24
# 运行Strelka2工作流
python2.7 ${sample_out_dir}/runWorkflow.py -m local -j $nt

# 合并和标准化VCF文件
bcftools concat ${sample_out_dir}/results/variants/somatic.indels.vcf.gz ${sample_out_dir}/results/variants/somatic.snvs.vcf.gz -a -Oz -o ${sample_out_dir}/$sample.strelka.vcf.gz
bcftools norm -m -both ${sample_out_dir}/$sample.strelka.vcf.gz -O z -o ${sample_out_dir}/$sample.strelka.norm.vcf.gz

