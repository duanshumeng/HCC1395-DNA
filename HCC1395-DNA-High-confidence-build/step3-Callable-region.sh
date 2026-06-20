cd /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/SEQC2_high_confidence/High_confidence_by_PGx/Callable_region
python /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/python_scripts/Filter_callable.py -d /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/SEQC2_high_confidence/High_confidence_by_PGx/Callable_region -p WGS_WUX_2021,BGI_T10_WGS_2021,BGI_T20_WGS_2021 -o Callable_isec

#Take the intersection of high-confidence intervals among different platforms and score them
bedtools multiinter -header -i *.bed > 2021.callable.bed

#Extract the regions that are considered callable regions by all three platforms
awk '$4 == 3' 2021.callable.bed > 2021.HCR.3.bed

#Merge the BED intervals
bedtools merge -d 1 -i 2021.HCR.3.bed > 2021.HCR.3.merge.bed

#Only retain chr1-22，X，Y
cat 2021.HCR.3.merge.bed | grep -v '_' > 2021.HCR.3.merge.filter.bed


#Compare with high-confidenced region of SEQC2
bedtools jaccard -a 2021.HCR.3.merge.filter.bed -b /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/SEQC2_high_confidence/High-Confidence_Regions_v1.2.bed
#Jaccard Index
0.900563
