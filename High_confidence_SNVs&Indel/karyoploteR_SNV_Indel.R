#高置信SNV/Indel在基因组中的分布情况
library(karyoploteR)
library(rtracklayer)
setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/SEQC2_analysis/hcc1395_analysis_pipeline/High_confidence_SNVs&Indel/')

snv_indel_karyoplote <- function(bed_file,indel_file,snv_file,out_name){
  gr <- import(bed_file)
  mcols(gr) <- data.frame(y=0.1,y0=0.1,y1=0.1)
  indel <- read.table(indel_file,sep = '\t',skip =42)
  indel_sub <- indel[,c('V1','V2','V2','V3','V4','V5')]
  colnames(indel_sub) <- c('chr','start','end','strand','ref','alt')
  indel_sub$mut.type <- 'Indel'
  snv <- read.table(snv_file,sep = '\t')
  snv_sub <- snv[,c('V1','V2','V2','V3','V4','V5')]
  colnames(snv_sub) <- c('chr','start','end','strand','ref','alt')
  snv_sub$mut.type <- 'SNV'
  
  indel.gr <- toGRanges(indel_sub)
  snv.gr <- toGRanges(snv_sub)
  #kpPlotRainfall(kp, data = sm.gr)
  pdf(paste0(out_name,'_karyoplot.pdf'))
  kp <- plotKaryotype(genome = "hg38",plot.type = 1,chromosomes = paste('chr',seq(1:22),sep=''))
  kpSegments(kp, data=gr, col="red")
  kpPlotCoverage(kp, data = indel.gr,col="green",r0=0.2,r1=0.7)
  kpPlotCoverage(kp, data = snv.gr,col="orange",r0=0.6,r1=0.9)
  dev.off()
}
###PGx
bed_file.pgx <- 'PGx/PGx_High-Confidence_Regions_v1.5.bed'
indel_file.pgx <- "PGx/high-confidence_sIndel_v1.sort.HCR.0919.vcf"
snv_file.pgx <- "PGx/high-confidence_sSNV_v1.sort.HCR.0919.vcf"

snv_indel_karyoplote(bed_file.pgx,indel_file.pgx,snv_file.pgx,'PGx')

###SEQC2
bed_file.seqc2 <- "SEQC2/High-Confidence_Regions_v1.2.bed"
indel_file.seqc2 <- "SEQC2/high-confidence_sINDEL_in_HC_regions_v1.2.1._sort.vcf.gz"
snv_file.seqc2 <- "SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1_sort.vcf.gz"

gunzip(indel_file.seqc2, remove = F) 
gunzip(snv_file.seqc2, remove = F) 

snv_indel_karyoplote(bed_file.seqc2,
                     gsub('.gz','',indel_file.seqc2),
                     gsub('.gz','',snv_file.seqc2),'SEQC2')

