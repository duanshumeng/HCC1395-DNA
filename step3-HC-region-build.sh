cd /data2/HCC1395/wux_WGS/somatic_SNVs/Callable_region
#First, extract the regions marked as CALLABLE.
for i in `ls *.bed`; do echo $i; cat $i | grep 'CALLABLE' > ${i/_status/} ; done

nohup python /data2/HCC1395/wux_WGS/somatic_SNVs/python_scripts/Filter_callable.py -d ./ -p WGS_WUX_2023,BGI_T10_WGS_2023,BGI_T20_WGS_2023 -o Callable_isec &

#Take the intersection of high - confidence intervals among different platforms and score them.
bedtools multiinter -header -i *.bed > PGx.callable.bed

#Extract the regions that are considered callable regions in at least 4 platforms.
awk '$4 > 3' PGx.callable.HCR.bed > PGx.HCR.bed
awk '$4 == 7' PGx.callable.HCR.bed > PGx.HCR.7.bed

#Merge the bed intervals
bedtools merge -d 1 -i PGx.HCR.7.bed > PGx.HCR.7.merge.bed

#Only retain chr1 - 22, X, and Y.
cat PGx.HCR.7.merge.bed | grep -v '_' > PGx.HCR.7.merge.filter.bed


#Compare with the high - confidence intervals of SEQC2.
Rename PGx.HCR.7.merge.filter.bed as PGx_High - Confidence_Regions_v1.3.bed.
bedtools jaccard -a PGx.HCR.7.merge.filter.bed -b /data2/HCC1395/tumor_normal_WGS/SEQC2_high_confidence/High-Confidence_Regions_v1.2.bed
#Jaccard index
0.884677

#RenamePGx.HCR.7.merge.filter.bed as PGx_High-Confidence_Regions_v1.3.bed
cp PGx.HCR.7.merge.filter.bed /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/High_confidence_datasets/PGx_High-Confidence_Regions_v1.3.bed

#Remove the regions 20 bp upstream and downstream of the sites judged as LowConf and Unclassified, 'chr6', 0, 59800000, and 'chr16', 36800000, 90338345.
#Output the results of high - confidence Indels
python /data2/HCC1395/wux_WGS/somatic_SNVs/python_scripts/get_LU_region20.py -i Merge.strelka.Indels.12.filter.vcf,Merge.TNscope.Indels.12.filter.vcf,Merge.TNseq.Indels.12.filter.vcf -t PGx -c 3 -v Indel -o Test_Indel

#Remove the intervals where low - confidence Indels are located.
cat PGx_High-Confidence_Regions_v1.3.bed | grep -v 'chrX\|chrY' >PGx_High-Confidence_Regions_tmp.bed
bedtools subtract -a PGx_High-Confidence_Regions_tmp.bed -b /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/High_confidence_datasets_12/UnConf.Indel.bed > PGx_High-Confidence_Regions_v1.4.bed

#和SEQC2的参考数据集比较一致性
bedtools subtract -a PGx_High-Confidence_Regions_tmp.bed -b /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/High_confidence_datasets_12/UnConf.Indel.bed > PGx_High-Confidence_Regions_v1.4.bed
jaccard -a SEQC2_High-Confidence_Regions_v1.2.bed -b PGx_High-Confidence_Regions_v1.4.bed
#0.940401

#Compare the consistency with the reference dataset of SEQC2.
cat high-confidence_sIndel_v1.sort.HighConf.vcf | grep '#' > high-confidence_sIndel_v1.sort.HighConf.HCRnew.vcf
bedtools intersect -a high-confidence_sIndel_v1.sort.HighConf.vcf -b /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/High_confidence_datasets/PGx_High-Confidence_Regions_v1.4.bed >> high-confidence_sIndel_v1.sort.HighConf.HCRnew.vcf

#Calculate F1 score：mean value is 0.76
for i in `ls *.vcf.gz`; do echo $i; docker run -i --rm -v /data2/HCC1395/wux_WGS/somatic_SNVs/:/data/ hap.py:v1 /opt/hap.py/bin/som.py /data/VCF_noPASS/High_confidence_datasets_12/High_confidence_sIndel_v1/high-confidence_sIndel_v1.sort.HighConf.HCRnew.vcf /data/VCF_noPASS/VCF_raw/${i} -r /data/GRCh38.d1.vd1.fa -R /data/VCF_noPASS/High_confidence_datasets/PGx_High-Confidence_Regions_v1.4.bed -o /data/VCF_noPASS/F1score_Results/${i/vcf.gz/Indel.F1}; done

#Output the results of high-confidenced SNV
python /data2/HCC1395/wux_WGS/somatic_SNVs/python_scripts/get_LU_region20.py -i Merge.strelka.SNVs.12.filter.vcf,Merge.TNscope.SNVs.12.filter.vcf,Merge.TNseq.SNVs.12.filter.vcf -t PGx -c 3 -v SNV -o Test_SNV
#Remobe SNVs within low-confidenced region
bedtools subtract -a PGx_High-Confidence_Regions_v1.4.bed -b /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/High_confidence_datasets_12/UnConf.SNV.bed > PGx_High-Confidence_Regions_v1.5.bed

#The consistency with SEQC2
jaccard -a SEQC2_High-Confidence_Regions_v1.2.bed -b PGx_High-Confidence_Regions_v1.5.bed
#0.921974

#Re - extract the variant sites within the high - confidence intervals
#Indel
cat high-confidence_sIndel_v1.sort.HighConf.HCRnew.vcf | grep '#' > high-confidence_sIndel_v1.sort.HighConf.HCRnew2.vcf
bedtools intersect -a high-confidence_sIndel_v1.sort.HighConf.HCRnew.vcf -b /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/High_confidence_datasets/PGx_High-Confidence_Regions_v1.5.bed >> high-confidence_sIndel_v1.sort.HighConf.HCRnew2.vcf
#Calcuate F1 score：
for i in `ls *.vcf.gz`; do echo $i; docker run -i --rm -v /data2/HCC1395/wux_WGS/somatic_SNVs/:/data/ hap.py:v1 /opt/hap.py/bin/som.py /data/VCF_noPASS/High_confidence_datasets_12/High_confidence_sIndel_v1/high-confidence_sIndel_v1.sort.HighConf.HCRnew2.vcf /data/VCF_noPASS/VCF_raw/${i} -r /data/GRCh38.d1.vd1.fa -R /data/VCF_noPASS/High_confidence_datasets/PGx_High-Confidence_Regions_v1.5.bed -o /data/VCF_noPASS/F1score_Results/${i/vcf.gz/Indel.F1}; done

#SNV
cd /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/High_confidence_datasets_12/High_confidence_sSNV_raw
cat high-confidence_sSNV_v1.sort.HighConf.vcf | grep '#' > high-confidence_sSNV_v1.sort.HighConf.HCR.vcf
bedtools intersect -a high-confidence_sSNV_v1.sort.HighConf.vcf -b /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/High_confidence_datasets/PGx_High-Confidence_Regions_v1.5.bed >> high-confidence_sSNV_v1.sort.HighConf.HCR.vcf
#calculate F1 score：
for i in `ls *.vcf.gz`; do echo $i; docker run -i --rm -v /data2/HCC1395/wux_WGS/somatic_SNVs/:/data/ hap.py:v1 /opt/hap.py/bin/som.py /data/VCF_noPASS/High_confidence_datasets_12/High_confidence_sSNV_raw/high-confidence_sSNV_v1.sort.HighConf.HCR.vcf /data/VCF_noPASS/VCF_raw/${i} -r /data/GRCh38.d1.vd1.fa -R /data/VCF_noPASS/High_confidence_datasets/PGx_High-Confidence_Regions_v1.5.bed -o /data/VCF_noPASS/F1score_Results/${i/vcf.gz/SNV.F1}; done

#Determine the final version of high-confidence
cd /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/High_confidence_datasets/High_confidence_datasets
cp /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/High_confidence_datasets_12/High_confidence_sIndel_v1/high-confidence_sIndel_v1.sort.vcf ./
cp /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/High_confidence_datasets_12/High_confidence_sSNV_raw/high-confidence_sSNV_v1.sort.vcf ./
#Indel
cat high-confidence_sIndel_v1.sort.vcf | grep '#' > high-confidence_sIndel_v1.sort.HCR.vcf
bedtools intersect -a high-confidence_sIndel_v1.sort.vcf -b /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/High_confidence_datasets/PGx_High-Confidence_Regions_v1.5.bed >> high-confidence_sIndel_v1.sort.HCR.vcf
#SNVs
cat high-confidence_sSNV_v1.sort.vcf | grep '#' > high-confidence_sSNV_v1.sort.HCR.vcf
bedtools intersect -a high-confidence_sSNV_v1.sort.vcf -b  /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/High_confidence_datasets/PGx_High-Confidence_Regions_v1.5.bed >> high-confidence_sSNV_v1.sort.HCR.vcf


#Screen out the SEQC2 sites that are missing using the voting method
#INDEL
docker run -it --rm -v /data2/HCC1395/tumor_normal_WGS/:/data/ hap.py:v1 bcftools isec -C /data/SEQC2_high_confidence/high-confidence_sINDEL_in_HC_regions_v1.2.1._sort.vcf.gz /data/SEQC2_high_confidence_PGx/high-confidence_sIndel_v1.sort.HighConf.vcf.gz -w 1 -o /data/SEQC2_high_confidence_PGx/sIndel_missing.vcf

#Count each confidence level
class=("HighConf" "MedConf" "LowConf" "Unclassified")
for i in "${class[@]}"; do echo $i; cat Merge.SNV.5.bed | grep $i | awk -F ' ' '{print $1,$2,$3,$5,$6,$7}' | uniq | wc -l; done
