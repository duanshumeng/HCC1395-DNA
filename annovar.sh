####### ANNOVAR databases: HPC 55 /data/annovar

## 1. GnomAD update: v3.12 only available for hg38
# https://annovar.openbioinformatics.org/en/latest/user-guide/download/
perl ./annotate_variation.pl -downdb -webfrom annovar -buildver hg38 gnomad312_genome humandb/
perl annotate_variation.pl -downdb -buildver hg38 -webfrom annovar gnomad312_genome humandb/

## 2. clinvar_20221231
perl ./annotate_variation.pl -downdb -webfrom annovar -buildver hg38 clinvar_20221231 humandb/

## 3. Cosmic update: v98
# https://cancer.sanger.ac.uk/cosmic
### Get token
$ echo "19210700118@fudan.edu.cn:LYQ1998lyq@" | base64
MTkyMTA3MDAxMThAZnVkYW4uZWR1LmNuOkxZUTE5OThseXFACg==

### CosmicCodingMuts
$ curl -H "Authorization: Basic MTkyMTA3MDAxMThAZnVkYW4uZWR1LmNuOkxZUTE5OThseXFACg==" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v95/VCF/CosmicCodingMuts.normal.vcf.gz

$ nohup curl "https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v95/VCF/CosmicCodingMuts.normal.vcf.gz?AWSAccessKeyId=KRV7P7QR9DL41J9EWGA2&Expires=1638783027&Signature=kjL0QxayZrCz%2BuCIwgi%2Fnt%2F%2BzHs%3D" --output CosmicCodingMuts.normal.vcf.gz &

$ curl -H "Authorization: Basic MTkyMTA3MDAxMThAZnVkYW4uZWR1LmNuOkxZUTE5OThseXFACg==" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v95/CosmicMutantExport.tsv.gz

$ nohup curl "https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v95/CosmicMutantExport.tsv.gz?AWSAccessKeyId=KRV7P7QR9DL41J9EWGA2&Expires=1638776709&Signature=QFIsFY6x6ZIn268l%2F%2FEomEBj8Go%3D" --output CosmicMutantExport.tsv.gz &

### CosmicNonCodingMuts
$ curl -H "Authorization: Basic MTkyMTA3MDAxMThAZnVkYW4uZWR1LmNuOkxZUTE5OThseXFACg==" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v95/VCF/CosmicNonCodingVariants.normal.vcf.gz

$ nohup curl "https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v95/VCF/CosmicNonCodingVariants.normal.vcf.gz?AWSAccessKeyId=KRV7P7QR9DL41J9EWGA2&Expires=1638779745&Signature=c7j96CqwTq5Ou7HIF7OZaNY3fCE%3D" --output CosmicNonCodingVariants.normal.vcf.gz &

$ curl -H "Authorization: Basic MTkyMTA3MDAxMThAZnVkYW4uZWR1LmNuOkxZUTE5OThseXFACg==" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v95/CosmicNCV.tsv.gz

$ nohup curl "https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v95/CosmicNCV.tsv.gz?AWSAccessKeyId=KRV7P7QR9DL41J9EWGA2&Expires=1638779891&Signature=1BozOY1lEZNseWtYTb8ETOuJWbM%3D" --output CosmicNCV.tsv.gz &

### unzip
gunzip CosmicCodingMuts.normal.vcf.gz
gunzip CosmicMutantExport.tsv.gz
gunzip CosmicNonCodingVariants.normal.vcf.gz
gunzip CosmicNCV.tsv.gz

perl ../prepare_annovar_user.pl -dbtype cosmic CosmicMutantExport.tsv -vcf CosmicCodingMuts.normal.vcf > hg38_cosmic95_coding.txt
perl ../prepare_annovar_user.pl -dbtype cosmic CosmicNCV.tsv -vcf CosmicNonCodingVariants.normal.vcf > hg38_cosmic95_noncoding.txt





####### ANNOTATION
## 记得更新clinvar, cosmic, gnomad的版本
nt=$(nproc)
sample="HCC1395_PGx_SNV"; vcf="./vcf/high-confidence_sSNV_v1.sort.HighConf.vcf"
/data/annovar/table_annovar.pl ${vcf} /data/annovar/humandb -buildver hg38 \
  -out ${sample} -remove \
  -protocol refGene,ensGene,clinvar_20220320,cosmic95_coding,cosmic95_noncoding,gnomad30_genome,dbnsfp42c \
  -operation g,g,f,f,f,f,f -nastring . -vcfinput -thread $nt

sample="HCC1395_PGx_INDEL"; vcf="./vcf/high-confidence_sIndel_v1.sort.HighConf.vcf"
/data/annovar/table_annovar.pl ${vcf} /data/annovar/humandb -buildver hg38 \
  -out ${sample} -remove \
  -protocol refGene,ensGene,clinvar_20220320,cosmic95_coding,cosmic95_noncoding,gnomad30_genome,dbnsfp42c \
  -operation g,g,f,f,f,f,f -nastring . -vcfinput -thread $nt

sample="HCC1395_SEQC2_SNV"; vcf="./vcf/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf"
/data/annovar/table_annovar.pl ${vcf} /data/annovar/humandb -buildver hg38 \
  -out ${sample} -remove \
  -protocol refGene,ensGene,clinvar_20220320,cosmic95_coding,cosmic95_noncoding,gnomad30_genome,dbnsfp42c \
  -operation g,g,f,f,f,f,f -nastring . -vcfinput -thread $nt

sample="HCC1395_SEQC2_INDEL"; vcf="./vcf/high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf"
/data/annovar/table_annovar.pl ${vcf} /data/annovar/humandb -buildver hg38 \
  -out ${sample} -remove \
  -protocol refGene,ensGene,clinvar_20220320,cosmic95_coding,cosmic95_noncoding,gnomad30_genome,dbnsfp42c \
  -operation g,g,f,f,f,f,f -nastring . -vcfinput -thread $nt
