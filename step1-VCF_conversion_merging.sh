

#TNscope-----------------------------------------
cd /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/VCF_raw
for i in `ls VCF_raw/*.TNscope.vcf.gz | awk -F '/' '{print $2}'`; do echo $i;  python /data2/HCC1395/wux_WGS/somatic_SNVs/python_scripts/modify_TNseq.py -infile VCF_raw/$i -snv Modify_TNscope/${i/vcf.gz/SNVs.vcf} -indel Modify_TNscope/${i/vcf.gz/Indels.vcf}; done

/data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/Modify_TNscope
# Compress the VCF file and build an index
for i in `ls *.vcf`; do echo $i; docker run -it --rm -v /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/Modify_TNscope:/data/ hap.py:v1 bcftools view /data/$i -Oz -o /data/${i/vcf/vcf.gz}; docker run -it --rm -v /data2/HCC1395/tumor_normal_WGS/SEQC2_VCF/bwa_VCF/Modify_muTect2:/data/ hap.py:v1 bcftools index /data/${i/vcf/vcf.gz}; done
for i in `ls *.gz`; do echo $i; docker run -it --rm -v /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/Modify_TNscope:/data/ hap.py:v1 bcftools index /data/$i; done

# Merge SNVs / Merge Indels
sh /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/Modify_TNscope/merge_SNV.sh
sh /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/Modify_TNscope/merge_Indel.sh

# Score the variant sites
#Indels
nohup python /data2/HCC1395/wux_WGS/somatic_SNVs/python_scripts/vcf_vote.py -i Merge.TNscope.Indels.vcf.gz -t PGx &

#SNVs
nohup python /data2/HCC1395/wux_WGS/somatic_SNVs/python_scripts/vcf_vote.py -i Merge.TNscope.SNVs.vcf.gz -t PGx &

#TNseq-----------------------------------------
cd /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/VCF_raw
#先对vcf进行norm
for i in `ls *.TNseq.vcf.gz`; do echo $i; docker run -it --rm -v /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/VCF_raw:/data/ hap.py:v1 bcftools norm -m -both /data/$i -O z -o /data/${i/vcf.gz/norm.vcf.gz}; done

for i in `ls VCF_raw/*.TNseq.norm.vcf.gz | awk -F '/' '{print $2}'`; do echo $i;  python /data2/HCC1395/wux_WGS/somatic_SNVs/python_scripts/modify_TNseq.py -infile VCF_raw/$i -snv Modify_TNseq/${i/norm.vcf.gz/SNVs.vcf} -indel Modify_TNseq/${i/norm.vcf.gz/Indels.vcf}; done

cd /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/Modify_TNseq
#压缩vcf文件并构建索引
for i in `ls *.vcf`; do echo $i; docker run -it --rm -v /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/Modify_TNseq:/data/ hap.py:v1 bcftools view /data/$i -Oz -o /data/${i/vcf/vcf.gz}; docker run -it --rm -v /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/Modify_TNseq:/data/ hap.py:v1 bcftools index /data/${i/vcf/vcf.gz}; done

#合并SNVs/#合并Indels
sh /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/Modify_TNseq/merge_SNV.sh

sh /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/Modify_TNseq/merge_Indel.sh


#对变异位点进行打分
#Indels
python /data2/HCC1395/wux_WGS/somatic_SNVs/python_scripts/vcf_vote.py -i Merge.TNseq.Indels.vcf.gz -t PGx

#SNVs
python /data2/HCC1395/wux_WGS/somatic_SNVs/python_scripts/vcf_vote.py -i Merge.TNseq.SNVs.vcf.gz -t PGx


#Strelka2-----------------------------------------
cd /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/VCF_raw
#修改文件格式并区分SNVs和Indels
for i in `ls VCF_raw/*.strelka.vcf.gz | awk -F '/' '{print $2}'`; do echo $i;  python /data2/HCC1395/wux_WGS/somatic_SNVs/python_scripts/modify_Strelka.py -infile VCF_raw/$i -snv Modify_Strelka/${i/vcf.gz/SNVs.vcf} -indel Modify_Strelka/${i/vcf.gz/Indels.vcf}; done

#压缩文件并构建索引
for i in `ls *.vcf`; do echo $i; docker run -it --rm -v /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/Modify_Strelka:/data/ hap.py:v1 bcftools view /data/$i -Oz -o /data/${i/vcf/vcf.gz}; docker run -it --rm -v /data2/HCC1395/tumor_normal_WGS/SEQC2_VCF/bwa_VCF/Modify_muTect2:/data/ hap.py:v1 bcftools index /data/${i/vcf/vcf.gz}; done

for i in `ls *.gz`; do echo $i; docker run -it --rm -v /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/Modify_Strelka:/data/ hap.py:v1 bcftools index /data/$i; done


#合并SNVs/#合并Indels
docker run -it --rm -v /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/Modify_Strelka:/data/ hap.py:v1 bcftools merge --force-samples /data/ARD_WGS_2023_1.strelka.Indels.vcf.gz /data/BGI_T10_WGS_2023_3.strelka.Indels.vcf.gz /data/BGI_T20_WGS_2023_3.strelka.Indels.vcf.gz /data/WGS_WUX_2021_3.strelka.Indels.vcf.gz /data/BGI_T10_WGS_2021_1.strelka.Indels.vcf.gz /data/BGI_T10_WGS_2021_2.strelka.Indels.vcf.gz /data/BGI_T10_WGS_2021_3.strelka.Indels.vcf.gz /data/BGI_T10_WGS_2023_1.strelka.Indels.vcf.gz /data/BGI_T10_WGS_2023_2.strelka.Indels.vcf.gz /data/BGI_T20_WGS_2021_1.strelka.Indels.vcf.gz /data/BGI_T20_WGS_2021_2.strelka.Indels.vcf.gz /data/BGI_T20_WGS_2021_3.strelka.Indels.vcf.gz /data/BGI_T20_WGS_2023_1.strelka.Indels.vcf.gz /data/BGI_T20_WGS_2023_2.strelka.Indels.vcf.gz /data/BGI_T7_WGS_2023_1.strelka.Indels.vcf.gz /data/NVG_WGS_2023_1.strelka.Indels.vcf.gz /data/WGE_WGS_2023_1.strelka.Indels.vcf.gz /data/WGS_WUX_2021_1.strelka.Indels.vcf.gz /data/WGS_WUX_2021_2.strelka.Indels.vcf.gz /data/WGS_WUX_2023_1.strelka.Indels.vcf.gz /data/WGS_WUX_2023_2.strelka.Indels.vcf.gz /data/WGS_WUX_2023_3.strelka.Indels.vcf.gz -Oz -o /data/Merge.strelka.Indels.vcf.gz
sh /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/Modify_Strelka/merge_SNV.sh
sh /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/Modify_Strelka/merge_Indel.sh

#对变异位点进行打分
#Indels
python /data2/HCC1395/wux_WGS/somatic_SNVs/python_scripts/vcf_vote.py -i Merge.strelka.Indels.vcf.gz -t PGx

#SNVs
python /data2/HCC1395/wux_WGS/somatic_SNVs/python_scripts/vcf_vote.py -i Merge.strelka.SNVs.vcf.gz -t PGx