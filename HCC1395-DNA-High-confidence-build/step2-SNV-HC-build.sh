cd /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS
for num in {12..12}
do
        cd /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS
        echo $num
        echo 'Vote...'
        python /data2/HCC1395/wux_WGS/somatic_SNVs/python_scripts/vcf_vote.py -i Modify_TNseq/Merge.TNseq.SNVs.vcf.gz -t PGx -n $num
        python /data2/HCC1395/wux_WGS/somatic_SNVs/python_scripts/vcf_vote.py -i Modify_Strelka/Merge.strelka.SNVs.vcf.gz -t PGx -n $num
        python /data2/HCC1395/wux_WGS/somatic_SNVs/python_scripts/vcf_vote.py -i Modify_TNscope/Merge.TNscope.SNVs.vcf.gz -t PGx -n $num

        mkdir High_confidence_datasets_${num} 

        echo 'Move file...'
        mv Modify_TNseq/*.$num.filter.vcf High_confidence_datasets_${num} 
        mv Modify_Strelka/*.$num.filter.vcf High_confidence_datasets_${num} 
        mv Modify_TNscope/*.$num.filter.vcf High_confidence_datasets_${num} 

        cd High_confidence_datasets_${num} 
        echo 'Filter HCR...' 
        for i in `ls *filter.vcf`; do echo $i; bedtools intersect -a $i -b /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/High_confidence_datasets/PGx_High-Confidence_Regions_v1.3.bed -header > ${i/vcf/HCR.vcf}; done 
        echo 'Filter Level1/2...' 
        for i in `ls *.HCR.vcf`; do echo $i; cat $i | grep -v '##' | grep -v 'Level3\|Level4' > ${i/vcf/Level12.vcf}; done 

        echo 'Get HC SNV...' 
        python /data2/HCC1395/wux_WGS/somatic_SNVs/python_scripts/get_vcf_info.py -i Merge.TNseq.SNVs.$num.filter.HCR.Level12.vcf,Merge.strelka.SNVs.$num.filter.HCR.Level12.vcf,Merge.TNscope.SNVs.$num.filter.HCR.Level12.vcf -t PGx -v SNV -c 3 -o High_confidence_sSNV 

        echo 'Jaccard index...' 
        cd 'High_confidence_sSNV' 
        docker run -i --rm -v /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/High_confidence_datasets_${num}/High_confidence_sSNV:/data/ hap.py:v1 bcftools sort /data/high-confidence_sSNV_v1.vcf -Oz -o /data/high-confidence_sSNV_v1.sort.vcf 
        bedtools jaccard -a /data2/HCC1395/tumor_normal_WGS/SEQC2_VCF/bwa_VCF/High_confidence_datasets/SNVS_candidates/High_confidence_sSNV/high-confidence_sSNV_in_HC_regions_v1.2.1_sort.vcf.gz -b high-confidence_sSNV_v1.sort.vcf > Jaccard_SNV_PGxvsSEQC2 
        cat high-confidence_sSNV_v1.sort.vcf | grep '#' > high-confidence_sSNV_v1.sort.HighConf.vcf
        cat high-confidence_sSNV_v1.sort.vcf | grep -v '#' | grep 'HighConf' >> high-confidence_sSNV_v1.sort.HighConf.vcf
        bedtools jaccard -a /data2/HCC1395/tumor_normal_WGS/SEQC2_VCF/bwa_VCF/High_confidence_datasets/SNVS_candidates/High_confidence_sSNV/high-confidence_sSNV_in_HC_regions_v1.2.1_sort.vcf.gz -b high-confidence_sSNV_v1.sort.HighConf.vcf > Jaccard_SNV_PGxvsSEQC2_HC 
        echo 'Filter Lavel1...'
        cd ../
        for i in `ls *.HCR.vcf`; do echo $i; cat $i | grep -v '##' | grep -v 'Level3\|Level4\|Level2' > ${i/vcf/Level1.vcf}; done
        python /data2/HCC1395/wux_WGS/somatic_SNVs/python_scripts/get_vcf_info.py -i Merge.TNseq.Indels.$num.filter.HCR.Level1.vcf,Merge.strelka.Indels.$num.filter.HCR.Level1.vcf,Merge.TNscope.Indels.$num.filter.HCR.Level1.vcf -t PGx -v Indel -c 3 -o High_confidence_sIndel_v1
        cd High_confidence_sIndel_v1
        docker run -i --rm -v /data2/HCC1395/wux_WGS/somatic_SNVs/VCF_noPASS/High_confidence_datasets_${num}/High_confidence_sIndel_v1:/data/ hap.py:v1 bcftools sort /data/high-confidence_sIndel_v1.vcf -Oz -o /data/high-confidence_sIndel_v1.sort.vcf
        bedtools jaccard -a /data2/HCC1395/tumor_normal_WGS/SEQC2_VCF/bwa_VCF/High_confidence_datasets/SNVS_candidates/High_confidence_sSNV/high-confidence_sSNV_in_HC_regions_v1.2.1_sort.vcf.gz -b high-confidence_sIndel_v1.sort.vcf > Jaccard_Indel_PGxvsSEQC2
        cat high-confidence_sIndel_v1.sort.vcf | grep '#' > high-confidence_sIndel_v1.sort.HighConf.vcf
        cat high-confidence_sIndel_v1.sort.vcf | grep -v '#' | grep 'HighConf' >> high-confidence_sIndel_v1.sort.HighConf.vcf
        bedtools jaccard -a /data2/HCC1395/tumor_normal_WGS/SEQC2_VCF/bwa_VCF/High_confidence_datasets/SNVS_candidates/High_confidence_sSNV/high-confidence_sSNV_in_HC_regions_v1.2.1_sort.vcf.gz -b high-confidence_sIndel_v1.sort.HighConf.vcf > Jaccard_Indel_PGxvsSEQC2_HC

        echo  'FInish...'
done
