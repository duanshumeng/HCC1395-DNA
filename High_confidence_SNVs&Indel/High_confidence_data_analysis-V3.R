setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/SEQC2_analysis/hcc1395_analysis_pipeline/High_confidence_SNVs&Indel')
#source('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/Quartet_Multiomics_Ratio_Code/utils/theme_nature.r')
my_theme <- theme_bw()+
  theme(plot.title = element_text(hjust=0.5,color = "black",face="plain"),
        axis.text.x =element_text(color = "black",size = 14,face = "bold",angle = 0),
        axis.text.y = element_text(color = "black",size = 14),
        axis.title.x = element_text(color = "black",face = "bold",size = 16),
        axis.title.y = element_text(color = "black",face = "bold",size = 16),
        #axis.ticks.x = element_text(color = "black",face = "bold",size = 14),
        legend.title = element_text(face = "bold",size = 16),
        legend.text = element_text(face = "bold",size = 16,colour = "black"),
        legend.position = "bottom",legend.background = element_blank(),
        panel.grid.minor  = element_line(colour = NA),
        panel.grid.major.x   = element_line(colour = NA),
        panel.background = element_rect(fill="transparent",colour = NA)
  )+theme(strip.background.x = element_rect(fill = "white", colour = "white")) +
  theme(strip.text.x = element_text(colour = "black",face = "bold",size = 16)) + 
  theme(strip.background.y = element_rect(fill = "white", colour = "white")) +
  theme(strip.text.y = element_text(colour = "black",face = "bold",size = 16)) + 
  theme(strip.placement = "inside") +
  theme(strip.switch.pad.grid = unit(1, "inch"))

var_colors = c('#8E1025','#598E44')
names(var_colors) <- c('SNVs','Indels')
library(patchwork)
library(stringr)
library(ggplot2)


hcr.SNV <- read.table("PGx_V2/high-confidence_sSNV_v1.sort.final.sort.vcf",sep = '\t')
hcr.SNV$tag <- paste0(hcr.SNV$V1,'&',hcr.SNV$V2,'&',hcr.SNV$V4,'&',hcr.SNV$V5)

hcr.Indel <- read.table("PGx_V2/high-confidence_sIndel_v1.sort.final.sort.vcf",sep = '\t')
hcr.Indel$tag <- paste0(hcr.Indel$V1,'&',hcr.Indel$V2,'&',hcr.Indel$V4,'&',hcr.Indel$V5)

hcr.SNV_Indel <- rbind(hcr.SNV,hcr.Indel)
#1. 统计数量-----------------------------
#当该变异位点的得分总和 Score_Sum 严格大于 24（不包含 24）时，被评级为 HighConf。
#当 Score_Sum 大于 22 且小于或等于 24 时，被评级为 MedConf。

hcr.Indel.df <- table(hcr.Indel$V7 %>% gsub('PASS;','',.)) %>% as.data.frame()
p.count.indel <- ggplot(data=hcr.Indel.df,aes(x=Var1,y=Freq,fill=Var1)) +
  geom_bar(stat = "identity",show.legend = FALSE)+
  my_theme+
  xlab(paste(''))+ylab(paste('sIndel','Count'))+
  xlim(c('HighConf','Rescued'))+ylim(c(0,max(hcr.Indel.df$Freq)+100))+
  geom_text(aes(label = Freq), vjust = -0.2)+
  scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")

hcr.SNV.df <- table(hcr.SNV$V7 %>% gsub('PASS;','',.)) %>% as.data.frame()
p.count.snv <- ggplot(data=hcr.SNV.df,aes(x=Var1,y=Freq,fill=Var1)) +
  geom_bar(stat = "identity",show.legend = FALSE)+
  my_theme+
  xlab(paste(''))+ylab(paste('sSNV','Count'))+
  xlim(c('HighConf','Rescued'))+ylim(c(0,max(hcr.SNV.df$Freq)+1000))+
  geom_text(aes(label = Freq), vjust = -0.2)+
  scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")

p.count.indel|p.count.snv
#2. 统计突变频率---------------------------

VAF_plot <- function(m.df,m.type,p.type){
  if (p.type=='SEQC2'){
    m.df$VAF <- stringr::str_split(m.df$V8,'\\;',simplify = T)[,18] %>% gsub("TVAF=","",.) %>% as.numeric()
    
  }else if (p.type=='PGx'){
    m.df$VAF <- stringr::str_split(m.df$V8,'\\;',simplify = T)[,5] %>% gsub("VAF=","",.) %>% as.numeric()
  }
  colnames(m.df)[7] <- "Type"
  m.df$Type <- gsub('PASS;','',m.df$Type)
  p_vaf_density = ggplot() +
    geom_histogram(data=m.df,aes(VAF,fill=Type),color="white",bins = 50)+
    my_theme+
    xlab(paste('VAF'))+ylab(paste(m.type,'Count'))+
    scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, by = 0.1))+
    scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")
  
  return(p_vaf_density)
}

p_vaf_density.hcr.Indel <- VAF_plot(hcr.Indel,'sIndel','PGx');p_vaf_density.hcr.Indel

p_vaf_density.hcr.SNV <- VAF_plot(hcr.SNV,'sSNV','PGx');p_vaf_density.hcr.SNV

p_summary <- ((p.count.snv / p.count.indel) | (p_vaf_density.hcr.SNV / p_vaf_density.hcr.Indel)) + 
  plot_layout(guides = 'collect',widths = c(1, 3))&theme(legend.position = 'bottom')

ggsave('PGx_V2/HC_variant_PGx.addLow.pdf',p_summary,width=8.15, height=5.7)


#Validated by WES---------------------------------
library(stringr)
library(ggplot2)
library(ggpubr)
#Element
wes = read.csv("Validated_by_WES/HCC1395_2023.isec2of3.validate.VAF.csv")
wes$tag <- paste(wes$X.CHROM,wes$POS,wes$REF,wes$ALT,sep='&') 
hcr.SNV_Indel$tag <- paste(hcr.SNV_Indel$V1,hcr.SNV_Indel$V2,hcr.SNV_Indel$V4,hcr.SNV_Indel$V5,sep='&')
hcr.SNV_Indel$TVAF <- str_split(hcr.SNV_Indel$V8,';',simplify = T)[,5] %>% gsub('VAF=','',.) %>% as.numeric()
hcr.SNV_Indel$Class <- gsub('PASS;','',hcr.SNV_Indel$V7)
validated.VAF <- merge(wes[,c('tag','VAF')],hcr.SNV_Indel[,c('tag','TVAF','Class')],by=c('tag'))
validated.VAF$VAF <- as.numeric(validated.VAF$VAF)

vaf.p <- ggplot(validated.VAF, aes(x = VAF, y = TVAF)) + theme_classic()+
  geom_point(aes(fill = Class),
             range_scale = .85,
             shape=21, size=3.5, color="white")+ scale_fill_brewer(palette = "Dark2")+
  geom_smooth(method = "lm", se = FALSE, color = "black") +stat_cor(method = "pearson",size = 5)+
  my_theme+coord_fixed()+
  labs(x = "VAF (Element WES set)", y = "VAF (WGS set)")

ggsave("Validated_by_WES/Validated.ELE.VAF.pdf",vaf.p,width = 5.4,height = 4.5)

#Illumina
wes_ilm = read.csv("Validated_by_WES/WES_HCC1395_ILM.isec2of3.validate.VAF.csv")
wes_ilm$tag <- paste(wes_ilm$X.CHROM,wes_ilm$POS,wes_ilm$REF,wes_ilm$ALT,sep='&') 
ilm.validated.VAF <- merge(wes_ilm[,c('tag','VAF')],
                           hcr.SNV_Indel[,c('tag','TVAF','Class')],by=c('tag'))

ilm.validated.VAF$VAF <- as.numeric(ilm.validated.VAF$VAF)

ilm.vaf.p <- ggplot(ilm.validated.VAF, aes(x = VAF, y = TVAF)) + theme_classic()+
  geom_point(aes(fill = Class),
             range_scale = .85,
             shape=21, size=3.5, color="white")+ scale_fill_brewer(palette = "Dark2")+
  geom_smooth(method = "lm", se = FALSE, color = "black") +stat_cor(method = "pearson",size = 5)+
  my_theme+coord_fixed()+
  labs(x = "VAF (Illumina WES set)", y = "VAF (WGS set)")

val_p <- vaf.p+ilm.vaf.p+plot_layout(guides='collect') & theme(legend.position='bottom',
                                                               plot.tag = element_text(color = "black",size = 18,face = "bold"))

ggsave("Validated_by_WES/Validated.ELE_ILM.VAF.pdf",val_p,,width=8.15, height=5.7*0.75)


#3. 统计注释类型-------------------
hcr_anno <- read.csv('PGx_V2/high-confidence_sSNV_Indel_v1.sort.final.sort.hg38_multianno.txt',sep='\t')

#4. 和SEQC2比较-------------------


