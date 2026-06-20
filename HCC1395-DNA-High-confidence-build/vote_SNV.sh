cd /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/High_confidence_datasets_v2
for num in {11..32}
do
        cd /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/High_confidence_datasets_v2
        echo $num
        echo 'Vote...'
        python /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/python_scripts/vcf_vote_with_group.py -i modify_Deepsomatic/Deepsomatic.merged.snv.vcf.gz -t PGx -n $num -m /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/High_confidence_datasets_v2/WGS_hcc1395.sample.txt
        
        python /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/python_scripts/vcf_vote_with_group.py -i modify_Strelka/Strelka.merged.snv.vcf.gz -t PGx -n $num -m /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/High_confidence_datasets_v2/WGS_hcc1395.sample.strelka2.txt
        
        python /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/python_scripts/vcf_vote_with_group.py -i modify_Mutect2/Mutect2.merged.snv.vcf.gz -t PGx -n $num -m /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/High_confidence_datasets_v2/WGS_hcc1395.sample.txt

        mkdir High_confidence_datasets_${num} 

        echo 'Move file...'
        mv modify_Deepsomatic/*.$num.filter.vcf High_confidence_datasets_${num} 
        mv modify_Strelka/*.$num.filter.vcf High_confidence_datasets_${num} 
        mv modify_Mutect2/*.$num.filter.vcf High_confidence_datasets_${num} 

        cd High_confidence_datasets_${num} 
        echo 'Filter HCR...' 
        for i in `ls *filter.vcf`; do echo $i; bedtools intersect -a $i -b /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/High_confidence_datasets/high-confidence_dataset2_V2/PGx_High-Confidence_Regions_v1.6.sorted.bed -header > ${i/vcf/HCR.vcf}; done 
        echo 'Filter Level1/2...' 
        for i in `ls *.HCR.vcf`; do echo $i; cat $i | grep -v '##' | grep -v 'Level3\|Level4' > ${i/vcf/Level12.vcf}; done 

        echo 'Get HC SNV...' 
        python /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/python_scripts/get_vcf_info.py -i Deepsomatic.merged.snv.$num.filter.HCR.Level12.vcf,Mutect2.merged.snv.$num.filter.HCR.Level12.vcf,Strelka.merged.snv.$num.filter.HCR.Level12.vcf -t PGx -v SNV -c 3 -o High_confidence_sSNV

        echo 'Jaccard index...' 
        bcftools sort High_confidence_sSNV/high-confidence_sSNV_v1.vcf > High_confidence_sSNV/high-confidence_sSNV_v1.sorted.vcf
        
        bedtools jaccard -a /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/High_confidence_datasets/high-confidence_dataset2_V2/high-confidence_sSNV_v2.HCR.clean.sort.vcf.gz -b High_confidence_sSNV/high-confidence_sSNV_v1.sorted.vcf > Jaccard_SNV_PGxvsPGx_V2

        echo  'FInish...'
done
