#SNVs/Indels High-confidence sets------------------
setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/SEQC2_analysis/hcc1395_analysis_pipeline/High_confidence_SNVs&Indel')
#最终结果展示
library(ggplot2)
library(dplyr)
snv <- data.frame(
    Category = c("HighConf", "MedConf"),
    Value = c(56381,2037)
)

indel <- data.frame(
  Category = c("HighConf", "MedConf"),
  Value = c(4601,415)
)

bar_plot <- function(data,v.type){
  p = ggplot(data, aes(x = Category, y = Value,fill=Category)) +
    geom_bar(stat = "identity",width=0.5) +
    geom_text(aes(label = Value), vjust = -0.5, size = 5)+
    labs( x = "Confidence level", y = v.type)+  
    theme_bw()+scale_fill_brewer(palette = "Set2")+scale_color_brewer(palette = "Set2")+
    theme(axis.text.x = element_text(angle=0,size = 25))+ 
    theme(axis.text.y = element_text(size = 25,color='black'))+
    theme(axis.title.y = element_text(size = 25,color='black'))+    
    theme(axis.title = element_text(size = 25),
          axis.text = element_text(size = 25),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 25))+guides(fill=FALSE)
  return(p)
}

snv_p = bar_plot(snv,'SNVs')
indel_p = bar_plot(indel,'Indels')

ggsave(snv_p,filename = 'HC_ref_summary.SNVs.png')
ggsave(indel_p,filename = 'HC_ref_summary.Indels.png')

#置信度
conf_df <- readxl::read_excel("SNV-Indel数量统计.xlsx" ,sheet = "置信度")
  
# 创建自定义的区间范围

df <- conf_df %>% mutate(Confidence_Range = cut(Confidence, breaks = interval_ranges, labels = FALSE))

conf_intervals <- cut(df$Confidence, breaks = seq(0.6,1,0.05))

result_snv <- df %>%
  group_by(Confidence_interval = conf_intervals) %>%
  summarize(Count = sum(SNV))

result_indel <- df %>%
  group_by(Confidence_interval = conf_intervals) %>%
  summarize(Count = sum(Indel))

conf_plot <- function(result,v_type){
  p = ggplot(result, aes(x = Confidence_interval, y = Count, fill = Confidence_interval)) +
    geom_bar(stat = "identity") +
    labs(x = "Confidence", y = v_type) +
    scale_fill_brewer(palette = "Set1")+    geom_text(aes(label = Count), vjust = -0.5, size = 5)+  
    theme_bw()+scale_fill_brewer(palette = "Set1")+scale_color_brewer(palette = "Set1")+
    theme(axis.text.x = element_text(angle=45,size = 15,hjust = 1))+ 
    theme(axis.text.y = element_text(size = 15,color='black'))+
    theme(axis.title.y = element_text(size = 15,color='black'))+    
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15))+guides(fill=FALSE)
  return(p)
}

p_snv_conf <- conf_plot(result_snv,'SNV')
p_indel_conf <- conf_plot(result_indel,'Indel')

ggsave(p_snv_conf,filename = 'HC_ref_confidence.SNVs.png')
ggsave(p_indel_conf,filename = 'HC_ref_confidence.Indels.png')


#-----SEQC2.vs.PGx------------------------------
##---差集SNV Type----------------------------
dir("PGx.vs.SEQC2/SNV")
dir("PGx.vs.SEQC2/Indel")
snv.pgx <- read.table("PGx.vs.SEQC2/SNV/Only.In.PGx.SNV.vcf",sep='\t')
snv.pgx$type = paste0(snv.pgx$V4,'>',snv.pgx$V5)
snv.pgx$type = gsub(',(.*)','',snv.pgx$type)
snv.pgx.type <- as.data.frame(table(snv.pgx$type) / length(snv.pgx$type) * 100 )
colnames(snv.pgx.type) <- c('SNV.Type','PGx')

snv.pgx.type <- snv.pgx.type[order(snv.pgx.type$PGx), ]

snv.seqc2 <- read.table("PGx.vs.SEQC2/SNV/Only.In.SEQC2.SNV.vcf",sep='\t')
snv.seqc2$type = paste0(snv.seqc2$V4,'>',snv.seqc2$V5)
snv.seqc2$type = gsub(',(.*)','',snv.seqc2$type)
snv.seqc2.type <- as.data.frame(table(snv.seqc2$type) / length(snv.seqc2$type) * 100 )
colnames(snv.seqc2.type) <- c('SNV.Type','SEQC2')


snv.type.df <- merge(snv.pgx.type,snv.seqc2.type,by='SNV.Type',all = F) %>% melt()
colnames(snv.type.df) <- c('SNV.Type','Batch','value')

library(reshape2)
p_snv_type = ggplot(snv.type.df, aes(x = SNV.Type,y=value,fill=Batch)) +
  labs(x = "SNV.Type", y = "Frequency") + xlim(snv.pgx.type$SNV.Type)+geom_bar(position = "dodge", stat = "identity", width = 0.5 )+
  theme_bw()+scale_fill_brewer(palette = "Set1")+scale_color_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle=45,size = 15,hjust = 1))+ 
  theme(axis.text.y = element_text(size = 15,color='black'))+
  theme(axis.title.y = element_text(size = 15,color='black'))+    
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))

ggsave(p_snv_type,filename = 'PGx.vs.SEQC2/SNVs.seqc2vspgx.png')


##----交集VAF相关性------------------


