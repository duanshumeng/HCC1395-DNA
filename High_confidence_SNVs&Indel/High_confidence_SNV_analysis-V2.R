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
#SEQC2-------------

f1_seqc2 = read.csv('SEQC2/F1score_stats_total_SEQC2.csv')

f1_seqc2_indel = subset(f1_seqc2,type =='indels' & (!grepl('.somaticSniper.',f1_seqc2$source)))
f1_seqc2_snv = subset(f1_seqc2,type =='SNVs' & (!grepl('.norm.',f1_seqc2$source)))


f1_plot_seqc2 <- function(F1_df,v.type){
  F1_df$Caller = str_split(F1_df$source,'\\.',simplify = T)[,2]
  
  F1_df = subset(F1_df,!grepl('merge|sort',F1_df$Caller))
  
  p1 = ggplot(F1_df,aes(x=Caller,y=F1.score,fill=Caller))+theme_classic()+
    geom_boxplot(notch=FALSE)+
    stat_summary(fun = mean, geom = "point", shape = 23, size = 5, fill = "red")+
    geom_text(aes(label = sprintf("Mean: %.2f", ..y..)), 
              stat = "summary", fun = mean, vjust = -1,size = 3)+
    my_theme+
    scale_fill_manual(values = var_colors)+
    #theme(axis.title = element_text(size = 20),
    #      axis.text = element_text(size = 20),
    #      legend.text = element_text(size = 20),
    #      legend.title = element_text(size = 20),
    #      strip.text = element_text(size = 20))+
    labs(x='',y='F1.score')+labs(title = v.type)+
    scale_fill_brewer(palette = 'Accent')+guides(fill=FALSE)+ylim(0,1)
  
  return(p1)
  
}

seqc2_indel.p = f1_plot_seqc2(f1_seqc2_indel,'sIndels');seqc2_indel.p

seqc2_snv.p = f1_plot_seqc2(f1_seqc2_snv,'sSNVs');seqc2_snv.p

ggsave('SEQC2/F1Score_sIndel_SEQC2.pdf',seqc2_indel.p)

ggsave('SEQC2/F1Score_sSNV_SEQC2.pdf',seqc2_snv.p)

#PGx--------------------------
f1_pgx = read.csv('PGx/F1score_stats_total_0922.csv')
f1_pgx$caller = str_split(f1_pgx$source,'\\.',simplify = T)[,2]
f1_pgx$source = paste0(f1_pgx$caller,'.',f1_pgx$source)
f1_pgx_indel = subset(f1_pgx,type =='indels' & (!grepl('.norm.|merge',f1_pgx$source)))
f1_pgx_snv = subset(f1_pgx,type =='SNVs' & (!grepl('.norm.|merge',f1_pgx$source)))


f1_plot_pgx <- function(F1_df,v.type){
  if (v.type=="sIndels"){
    hline=0.7
  } else if (v.type == 'sSNVs'){
    hline=0.81
  }
  F1_df$Caller = str_split(F1_df$source,'\\.',simplify = T)[,2]
  
  F1_df = subset(F1_df,!grepl('merge|sort|Indel|SNV',F1_df$Caller))
  
  p1 = ggplot(F1_df,aes(x=Caller,y=F1.score,fill=Caller))+theme_classic()+
    geom_boxplot(notch=FALSE)+
    stat_summary(fun = mean, geom = "point", shape = 23, size = 5, fill = "red")+
    geom_text(aes(label = sprintf("Mean: %.2f", ..y..)), 
              stat = "summary", fun = mean, vjust = -1,size = 5)+
    geom_hline(aes(yintercept=hline), colour="Red", linetype="dashed")+
    theme_nature_xy()+
    scale_fill_manual(values = var_colors)+
    #theme(axis.title = element_text(size = 20),
    #      axis.text = element_text(size = 20),
    #      legend.text = element_text(size = 20),
    #      legend.title = element_text(size = 20),
    #      strip.text = element_text(size = 20))+
    labs(x='',y='F1.score')+labs(title = v.type)+
    scale_fill_brewer(palette = 'Accent')+guides(fill=FALSE)+ylim(0,1)+xlim(c("strelka","TNseq","TNscope"))
  
  return(p1)
  
}

pgx_indel.p = f1_plot_pgx(f1_pgx_indel,'sIndels')

pgx_snv.p = f1_plot_pgx(f1_pgx_snv,'sSNVs')

ggsave('PGx/F1Score_sIndel_PGx.pdf',pgx_indel.p)

ggsave('PGx/F1Score_sSNV_PGx.pdf',pgx_snv.p)

F1_df_make <- function(vcf_count){
  vcf_count$source <- gsub('WGS_WUX','WUX_WGS',vcf_count$source)
  vcf_count['Platform'] <- str_split(vcf_count$source,'\\_WGS',simplify = T)[,1] %>% gsub('strelka.|TNseq.|TNscope.','',.)
  vcf_count['Caller'] <- vcf_count$caller
  lib_p <- c('PCR-free','PCR-free','PCR','PCR-free','PCR-free','PCR-free','PCR-free')
  comp <- c('BGI_T20','BGI_T10','WUX','BGI_T7','WGE','ARD','NVG')
  tech <- c('DNBSEQ-T20','DNBSEQ-T10','NovaSeq','DNBSEQ-T7','DNBSEQ-T7','DNBSEQ-T7','NovaSeq')
  frag <-c('Enzyme','Enzyme','Ultrasounds','Enzyme','Ultrasounds','Enzyme','Ultrasounds')
  
  info_df <- data.frame(Lib_construction=lib_p,
                        Platform=comp,
                        Technology=tech,
                        Fragment=frag)
  
  vcf_count <- merge(vcf_count,info_df,by=c('Platform'),all = TRUE)
  
  vcf_count$Platform<- as.character(vcf_count$Platform)
  vcf_count$Type<- as.character(vcf_count$type)
  vcf_count$Lib_construction<- as.character(vcf_count$Lib_construction)
  vcf_count$Technology<- as.character(vcf_count$Technology)
  vcf_count$Fragment<- as.character(vcf_count$Fragment)
  vcf_count$Caller<- as.character(vcf_count$Caller)
  vcf_count$source <- paste0(vcf_count$Caller,'.',vcf_count$source)
  return(vcf_count)
}

#F1 score barplot
F1_barplot <- function(vcf_count,region){
  med_f1 = median(vcf_count$F1.score)
  print(med_f1)
  h1 <- ggplot(data=vcf_count, aes(x=source, y=F1.score, fill=Caller)) +
    geom_bar(stat = "identity", position = "dodge",color="black")+
    geom_hline(aes(yintercept=med_f1), colour="Red", linetype="dashed")+
    scale_fill_brewer(palette ='Pastel2')+theme_nature_border()+
    theme(axis.text = element_text(size = 13),
          axis.title = element_text(size = 20))+
    theme(axis.text.x = element_blank())+xlab('')+labs(title = region)
  
  
  h7 <- ggplot(vcf_count)+
    geom_bar(mapping = aes(x = source, y = 1, fill = Lib_construction), 
             stat = "identity", 
             width = 1)+
    scale_fill_brewer(palette ='Dark2')+
    theme_void()+
    theme(panel.spacing.x = unit(1, "mm"))
  
  
  h5 <- ggplot(vcf_count)+
    geom_bar(mapping = aes(x = source, y = 1, fill = Technology), 
             stat = "identity", 
             width = 1)+
    scale_fill_brewer(palette ='Paired')+
    theme_void()+
    theme(panel.spacing.x = unit(1, "mm"))
  
  h6 <- ggplot(vcf_count)+
    geom_bar(mapping = aes(x = source, y = 1, fill = Fragment), 
             stat = "identity", 
             width = 1)+
    scale_fill_brewer(palette ='Pastel1')+
    theme_void()+
    theme(panel.spacing.x = unit(1, "mm"))
  
  h2 <- ggplot(vcf_count)+
    geom_bar(mapping = aes(x = source, y = 1, fill = Caller), 
             stat = "identity", 
             width = 1)+
    scale_fill_brewer(palette ='Pastel2')+
    theme_void()+
    theme(panel.spacing.x = unit(1, "mm"))
  
  h4 <- ggplot(vcf_count)+
    geom_bar(mapping = aes(x = source, y = 1, fill = Platform), 
             stat = "identity", 
             width = 1)+
    scale_fill_brewer(palette ='Set1')+
    theme_void()+
    theme(panel.spacing.x = unit(1, "mm"))
  
  h3 <- ggplot(vcf_count)+
    geom_bar(mapping = aes(x = source, y = 1, fill = Batch), 
             stat = "identity", 
             width = 1)+
    scale_fill_brewer(palette ='Set2')+
    theme_void()+
    theme(panel.spacing.x = unit(1, "mm"))
  
  
  legend <- plot_grid(get_legend(h1),get_legend(h6), 
                      get_legend(h4),get_legend(h5),get_legend(h7),label_size = 5,
                      nrow = 1)
  h1 <- h1 + theme(legend.position = "none")
  #h2 <- h2 + theme(legend.position = "none")
  #h3 <- h3 + theme(legend.position = "none")
  h4 <- h4 + theme(legend.position = "none")
  h5 <- h5 + theme(legend.position = "none")
  h6 <- h6 + theme(legend.position = "none")
  h7 <- h7 + theme(legend.position = "none")
  plot <- plot_grid(h1,h4,h5,h6,h7,align = "v", ncol = 1, axis = "tb", rel_heights = c(10, 0.6,0.6,0.6,0.6))
  all_plot <- plot_grid(plot, legend,align = "v", ncol = 1, axis = "tb",rel_heights = c(2.5,1))
  return(all_plot)
  
}

f1_pgx_indel.df <- F1_df_make(f1_pgx_indel)
p.f1.indel <- F1_barplot(f1_pgx_indel.df,'Indels')
ggsave(p.f1.indel,filename = "F1score.median.Indel.pdf",width = 9.4,height = 6.7)

f1_pgx_snv.df <- F1_df_make(f1_pgx_snv)
p.f1.snv <- F1_barplot(f1_pgx_snv.df,'SNVs')
ggsave(p.f1.snv,filename = "F1score.median.SNV.pdf",width = 9.4,height = 6.7)

#Indel------------------------

#Check missing variant using SEQC datasets--------------
m.Indel <- read.table("sIndel_missing.vcf",sep = '\t')

m.SNV <- read.table("sSNV_missing.vcf",sep="\t")

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


p_vaf_density.indel <- VAF_plot(m.Indel,'sIndel','SEQC2')

p_vaf_density.snv <- VAF_plot(m.SNV,'sSNV','SEQC2')

ggsave('SEQC2/VAF_sIndel_SEQC2.missing.pdf',p_vaf_density.indel,width =4.4,height = 2.2,dpi = 300)

ggsave('SEQC2/VAF_sSNV_SEQC2.missing.pdf',p_vaf_density.snv,width =4.4,height = 2.2,dpi = 300)

#Check missing variant using SEQC datasets produced by PGx--------------
m.indel.pgx <- read.table("highconfidence_datasets_batch2021/high-confidence_sIndel_v1.sort.HighConf.HCR.missing.vcf",sep = '\t')
m.snv.pgx <- read.table("highconfidence_datasets_batch2021/high-confidence_sSNV_v1.sort.HighConf.HCR.missing.vcf",sep = '\t')
p_vaf_density.indel.pgx <- VAF_plot(m.indel.pgx,'sIndel','SEQC2')
p_vaf_density.snv.pgx <- VAF_plot(m.snv.pgx,'sSNV','SEQC2')

ggsave('highconfidence_datasets_batch2021/VAF_sIndel_SEQC2.byPGx.missing.pdf',p_vaf_density.indel.pgx,width =4.4,height = 2.2,dpi = 300)

ggsave('highconfidence_datasets_batch2021/VAF_sSNV_SEQC2.byPGx.missing.pdf',p_vaf_density.snv.pgx,width = 4.4,height = 2.2,dpi = 300)

#绘制VAF分布图------------------------------------------------------
#hcr.Indel <- read.table("PGx/high-confidence_sIndel_v1.sort.HCR.0919.vcf",sep = '\t',skip = 42)
#hcr.Indel.low <- read.table("PGx/BGI_Nova.VAFlt0.2.sIndel.HCR.vcf",sep = '\t',skip = 38)
#hcr.Indel <- rbind(hcr.Indel,hcr.Indel.low)

#hcr.SNV <- read.table("PGx/high-confidence_sSNV_v1.sort.HCR.0919.vcf",sep="\t")
#hcr.SNV.low <- read.table("PGx/BGI_Nova.VAFlt0.2.sSNV.HCR.vcf",sep="\t",skip = 125)
#hcr.SNV <- rbind(hcr.SNV,hcr.SNV.low)
hcr.Indel <- read.table("PGx/high-confidence_sIndel_v2.HCR.add_cred.vcf",sep = '\t')
hcr.SNV <- read.table("PGx/high-confidence_sSNV_v2.HCR.add_cred.vcf",sep="\t")

p_vaf_density.hcr.Indel <- VAF_plot(hcr.Indel,'sIndel','PGx');p_vaf_density.hcr.Indel

p_vaf_density.hcr.SNV <- VAF_plot(hcr.SNV,'sSNV','PGx');p_vaf_density.hcr.SNV


#统计SNV/Indel的数量------
hcr.Indel.df <- table(hcr.Indel$V7 %>% gsub('PASS;','',.)) %>% as.data.frame()
p.count.indel <- ggplot(data=hcr.Indel.df,aes(x=Var1,y=Freq,fill=Var1)) +
  geom_bar(stat = "identity",show.legend = FALSE)+
  my_theme+
  xlab(paste(''))+ylab(paste('sIndel','Count'))+
  xlim(c('HighConf','MedConf','LowConf','Unclassified'))+ylim(c(0,max(hcr.Indel.df$Freq)+100))+
  geom_text(aes(label = Freq), vjust = -0.2)+
  scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")

hcr.SNV.df <- table(hcr.SNV$V7 %>% gsub('PASS;','',.)) %>% as.data.frame()
p.count.snv <- ggplot(data=hcr.SNV.df,aes(x=Var1,y=Freq,fill=Var1)) +
  geom_bar(stat = "identity",show.legend = FALSE)+
  my_theme+
  xlab(paste(''))+ylab(paste('sSNV','Count'))+
  xlim(c('HighConf','MedConf','LowConf','Unclassified'))+ylim(c(0,max(hcr.SNV.df$Freq)+1000))+
  geom_text(aes(label = Freq), vjust = -0.2)+
  scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")

p.count.all <- ((p.count.indel/p.count.snv)|(p_vaf_density.hcr.Indel/p_vaf_density.hcr.SNV))+plot_layout(widths = c(1,2))+plot_layout(guides='collect') & theme(legend.position='bottom',
                                                                                        plot.tag = element_text(color = "black",size = 18,face = "bold"))

p.count.all
ggsave('PGx/HC_variant_PGx.addLow.pdf',p.count.all,width=8.15, height=5.7)

#ggsave('PGx/Count_sIndel_PGx.addLow.pdf',p.count.indel,width=8.1*0.6, height=2.7)

#ggsave('PGx/Count_sSNV_PGx.addLow.pdf',p.count.snv,width=8.1*0.6, height=2.7)

#Compare TMB----------------------------------------------
library(stringr)
tmb_df <- read.csv("TMB.HCC1395.csv")
tmb_df$source <- gsub('WGS_WUX','WUX_WGS',tmb_df$source)
tmb_df$Platform <- str_split(tmb_df$source,'_WGS_',simplify = T)[,1]
tmb_df$Batch <- str_split(str_split(tmb_df$source,'_WGS_',simplify = T)[,2],'\\_',simplify = T)[,1]
  
#同一批次的技术重复
#2021批次
p.2021 <- ggplot(subset(tmb_df,Batch=='2021'),aes(x=Platform,y= TMB,fill=Platform))+theme_classic()+
  geom_boxplot(notch=FALSE)+
  stat_summary(fun = mean, geom = "point", shape = 23, size = 5, fill = "red")+
  geom_text(aes(label = sprintf("Mean: %.2f", ..y..)), 
            stat = "summary", fun = mean, vjust = -1,size = 5)+
  theme_nature_border()+
  #theme(axis.title = element_text(size = 20),
  #      axis.text = element_text(size = 20),
  #      legend.text = element_text(size = 20),
  #      legend.title = element_text(size = 20),
  #      strip.text = element_text(size = 20))
  labs(x='',y='TMB')+ggtitle('Batch 2021')+
  scale_fill_brewer(palette = 'Accent')+guides(fill=FALSE)+ylim(15,18)

#2023批次
p.2023 <- ggplot(subset(tmb_df,Batch=='2023'),aes(x=Platform,y= TMB,fill=Platform))+theme_classic()+
  geom_boxplot(notch=FALSE)+
  stat_summary(fun = mean, geom = "point", shape = 23, size = 5, fill = "red")+
  geom_text(aes(label = sprintf("Mean: %.2f", ..y..)), 
            stat = "summary", fun = mean, vjust = -1,size = 5)+
  #theme(axis.title = element_text(size = 20),
  #      axis.text = element_text(size = 20),
  #      legend.text = element_text(size = 20),
  #      legend.title = element_text(size = 20),
  #      strip.text = element_text(size = 20))+
  theme_nature_border()+
  labs(x='',y='TMB')+ggtitle('Batch 2023')+
  scale_fill_brewer(palette = 'Accent')+guides(fill=FALSE)+ylim(15,18)+xlim(c('BGI_T7','BGI_T20','WUX','WGE','ARD','BGI_T10','NVG'))


#2021vs2023批次
p.2021vs2023 <- ggplot(subset(tmb_df,grepl('WUX|BGI_T20|BGI_T10',tmb_df$Platform)),aes(x=Batch,y= TMB,fill=Batch))+theme_classic()+
  geom_boxplot(notch=FALSE)+
  stat_summary(fun = mean, geom = "point", shape = 23, size = 5, fill = "red")+
  geom_text(aes(label = sprintf("Mean: %.2f", ..y..)), 
            stat = "summary", fun = mean, vjust = -1,size = 5)+
  theme_nature_border()+
  #theme(axis.title = element_text(size = 20),
  #      axis.text = element_text(size = 20),
  #      legend.text = element_text(size = 20),
  #      legend.title = element_text(size = 20),
  #      strip.text = element_text(size = 20))
  labs(x='',y='TMB')+ggtitle('Batch 2021vs2023')+
  scale_fill_brewer(palette = 'Accent')+guides(fill=FALSE)+ylim(15,18)


ggsave('PGx/TMB_boxplot_2021.pdf',p.2021)

ggsave('PGx/TMB_boxplot_2023.pdf',p.2023)

ggsave('PGx/TMB_boxplot_2021vs2023.pdf',p.2021vs2023)

#SEQC2得到的参考数据集与我们构建的参考数据集的F1 score比较-----------
library(ggpubr)
F1.seqc2 <- read.csv("SEQC2/F1score_stats_total_SEQC2.csv")
F1.seqc2$sample <- str_split(F1.seqc2$source,'.vcf.',simplify = T)[,1]
F1.seqc2.byPGx <- read.csv("SEQC2/F1score_stats_total_SEQC2.byPGx.csv")
F1.seqc2.byPGx$sample <- gsub('.reh','',F1.seqc2.byPGx$source) %>% gsub(
                          '.Indel.F1','',.) %>% gsub(
                            '.SNV.F1','',.)
F1.seqc2 <- subset(F1.seqc2,grepl('bwa',F1.seqc2$sample))
F1.seqc2$sample <- gsub('.bwa','',F1.seqc2$sample)

F1.seqc2.byPGx <- subset(F1.seqc2.byPGx,type!="records")

F1.combine <- merge(F1.seqc2.byPGx[,c('type','sample','F1.score')],
      F1.seqc2[,c('type','sample','F1.score')],
      by=c('type','sample'),suffixes =c('.PGx','.SEQC2'))
F1.combine$Caller <- str_split(F1.combine$sample,'\\.',simplify = T)[,2]
F1.combine$type <- gsub('indels','Indels',F1.combine$type)

F1.plot <- function(v_type){
  p_f1 <- ggplot(subset(F1.combine,type==v_type), aes(x = F1.score.PGx, y = F1.score.SEQC2)) + theme_classic()+
    geom_point(aes(color=Caller),size=10)+ 
    geom_smooth(method = "lm", se = FALSE, color = "black") +stat_cor(method = "pearson",p.accuracy = 0.0001,size = 5)+
    coord_fixed()+
    theme_nature_border()+
    #theme(
    #  panel.border = element_rect(color = "black", fill = NA, size = 1),
    #  panel.background = element_blank(),
    #  plot.margin = unit(c(1, 1, 1, 1), "cm"),
    #  text = element_text(size = 12),  # 设置默认文本大小
    #  axis.title = element_text(size = 14),  # 设置坐标轴标题字体大小
    #  axis.text = element_text(size = 12),  # 设置坐标轴刻度标签字体大小
    #  legend.title = element_text(size = 14),  # 设置图例标题字体大小
    #  legend.text = element_text(size = 12)  # 设置图例标签字体大小
    #)
    labs(x = "F1.score.PGx", y = "F1.score.SEQC2",title = v_type)+
    scale_color_brewer(palette = 'Accent')+
    xlim(c(0,1))+ylim(c(0,1))
}
p.f1.indel <- F1.plot('Indels')
p.f1.snv <- F1.plot('SNVs')

ggsave(p.f1.indel,filename = "SEQC2/F1_SEQC2vsPGx.Indel.pdf")
ggsave(p.f1.snv,filename = "SEQC2/F1_SEQC2vsPGx.SNV.pdf")

#Validated by WES---------------------------------
library(stringr)
library(ggplot2)
library(ggpubr)
#Element
wes = read.csv("Validated_by_WES/HCC1395_2023.isec2of3.validate.VAF.csv")
wes$tag <- paste(wes$X.CHROM,wes$POS,wes$REF,wes$ALT,sep='&') 
hc.indel.snv <- read.table("Validated_by_WES/high-confidence_sIndel_sSNV_v2.HCR.clean.sort.WES.vcf",sep = '\t')
hc.indel.snv$tag <- paste(hc.indel.snv$V1,hc.indel.snv$V2,hc.indel.snv$V4,hc.indel.snv$V5,sep='&')
hc.indel.snv$TVAF <- str_split(hc.indel.snv$V8,';',simplify = T)[,5] %>% gsub('VAF=','',.) %>% as.numeric()
hc.indel.snv$Class <- gsub('PASS;','',hc.indel.snv$V7)
validated.VAF <- merge(wes[,c('tag','VAF')],hc.indel.snv[,c('tag','TVAF','Class')],by=c('tag'))
validated.VAF$VAF <- as.numeric(validated.VAF$VAF)

vaf.p <- ggplot(validated.VAF, aes(x = VAF, y = TVAF)) + theme_classic()+
  geom_point(aes(fill = Class),
             range_scale = .85,
             shape=21, size=3.5, color="white")+ scale_fill_brewer(palette = "Dark2")+
  geom_smooth(method = "lm", se = FALSE, color = "black") +stat_cor(method = "pearson",size = 5)+
  my_theme+
  labs(x = "VAF (Element WES set)", y = "VAF (WGS set)")

ggsave("Validated_by_WES/Validated.ELE.VAF.pdf",vaf.p,width = 5.4,height = 4.5)

#Illumina
wes_ilm = read.csv("Validated_by_WES/WES_HCC1395_ILM.isec2of3.validate.VAF.csv")
wes_ilm$tag <- paste(wes_ilm$X.CHROM,wes_ilm$POS,wes_ilm$REF,wes_ilm$ALT,sep='&') 
ilm.validated.VAF <- merge(wes_ilm[,c('tag','VAF')],
                           hc.indel.snv[,c('tag','TVAF','Class')],by=c('tag'))

ilm.validated.VAF$VAF <- as.numeric(ilm.validated.VAF$VAF)

ilm.vaf.p <- ggplot(ilm.validated.VAF, aes(x = VAF, y = TVAF)) + theme_classic()+
  geom_point(aes(fill = Class),
             range_scale = .85,
             shape=21, size=3.5, color="white")+ scale_fill_brewer(palette = "Dark2")+
  geom_smooth(method = "lm", se = FALSE, color = "black") +stat_cor(method = "pearson",size = 5)+
  my_theme+
  labs(x = "VAF (Illumina WES set)", y = "VAF (WGS set)")

val_p <- vaf.p+ilm.vaf.p+plot_layout(guides='collect') & theme(legend.position='bottom',
                                                                                   plot.tag = element_text(color = "black",size = 18,face = "bold"))

ggsave("Validated_by_WES/Validated.ELE_ILM.VAF.pdf",val_p,,width=8.15, height=5.7*0.75)

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
  if (max(data$Value)<10000){
    y_max = max(data$Value) + 1000
  } else {
    y_max = max(data$Value) + 10000
  }
  
  p = ggplot(data, aes(x = Category, y = Value,fill=Category)) +
    geom_bar(stat = "identity",width=0.5) +
    geom_text(aes(label = Value), vjust = -0.5, size = 5)+
    labs( x = "Confidence level", y = v.type)+  
    theme_bw()+scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")+
    ylim(0,y_max)+
    theme_nature_border()+guides(fill=FALSE)
  return(p)
}

snv_p = bar_plot(snv,'SNVs')
indel_p = bar_plot(indel,'Indels')

ggsave(snv_p,filename = 'HC_ref_summary.SNVs.pdf',width=2.7, height=2.7)
ggsave(indel_p,filename = 'HC_ref_summary.Indels.pdf',width=2.7, height=2.7)

#置信度
#conf_df <- readxl::read_excel("SNV-Indel数量统计.xlsx" ,sheet = "置信度")
conf_df <- read.csv('Credibility_HC.csv')
# 创建自定义的区间范围

#df <- conf_df %>% mutate(Confidence_Range = cut(Confidence, breaks = interval_ranges, labels = FALSE))

conf_intervals <- cut(subset(conf_df,Type=='SNVs')$Credibility %>% sort(), breaks = seq(60,100,5))

result_snv <- subset(conf_df,Type=='SNVs') %>%
  group_by(Confidence_interval = conf_intervals) %>%
  summarize(SNVs = sum(Count))

conf_intervals.i <- cut(subset(conf_df,Type=='Indels')$Credibility %>% sort(), breaks = seq(60,100,5))

result_indel <- subset(conf_df,Type=='Indels')  %>%
  group_by(Confidence_interval = conf_intervals.i) %>%
  summarize(Indels = sum(Count))


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

ggsave(p_snv_conf,filename = 'HC_ref_confidence.SNVs.pdf')
ggsave(p_indel_conf,filename = 'HC_ref_confidence.Indels.pdf')


#-----SEQC2.vs.PGx------------------------------
##---差集SNV Type----------------------------
dir("PGx.vs.SEQC2")
dir("PGx.vs.SEQC2")
snv.pgx <- read.table("PGx.vs.SEQC2/PGx_only.vcf",sep='\t')
snv.pgx$type = paste0(snv.pgx$V4,'>',snv.pgx$V5)
snv.pgx$VAF= str_split(snv.pgx$V8,';',simplify = T)[,5] %>% gsub('VAF=','',.) %>% as.numeric()
snv.pgx$Class = 'PGx'
snv.pgx$type = gsub(',(.*)','',snv.pgx$type)
snv.pgx.type <- as.data.frame(table(snv.pgx$type) / length(snv.pgx$type) * 100 )
colnames(snv.pgx.type) <- c('SNV.Type','VAF')

snv.pgx.type <- snv.pgx.type[order(snv.pgx.type$VAF), ]


snv.seqc2 <- read.table("PGx.vs.SEQC2/SEQC2_only.vcf",sep='\t')
snv.seqc2$type = paste0(snv.seqc2$V4,'>',snv.seqc2$V5)
snv.seqc2_VAF <- lapply(snv.seqc2$V8,function(x){
  print(x)
  tvaf = grep("^TVAF=", unlist(strsplit(x, ";")), value = TRUE)
  print(tvaf)
  if (is.null(tvaf)){
    return(0)
  } else {
    return(gsub('TVAF=','',tvaf))
  }
}) %>% unlist()
snv.seqc2$VAF <- c(snv.seqc2_VAF,0) %>% as.numeric()
snv.seqc2$Class <- 'SEQC2'

snv.seqc2$type = gsub(',(.*)','',snv.seqc2$type)
snv.seqc2.type <- as.data.frame(table(snv.seqc2$type) / length(snv.seqc2$type) * 100 )
colnames(snv.seqc2.type) <- c('SNV.Type','SEQC2')


snv.type.df <- merge(snv.pgx.type,snv.seqc2.type,by='SNV.Type',all = F) %>% melt()
colnames(snv.type.df) <- c('SNV.Type','Batch','value')

library(reshape2)
p_snv_type = ggplot(snv.type.df, aes(x = SNV.Type,y=value,fill=Batch)) +
  labs(x = "SNV.Type", y = "Frequency (%)") + xlim(snv.pgx.type$SNV.Type)+
  geom_bar(position = "dodge", stat = "identity", width = 0.5,color = "black" )+
  theme_bw()+scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle=45,size = 15,hjust = 1))+ 
  theme(axis.text.y = element_text(size = 15,color='black'))+
  theme(axis.title.y = element_text(size = 15,color='black'))+    
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))

ggsave(p_snv_type,filename = 'PGx.vs.SEQC2/SNVs.type.diff.pdf')



##---交集SNV Type----------------------------
snv.pgx.isec <- read.table("PGx.vs.SEQC2/PGx_commonWith_SEQC2.vcf",sep='\t')

snv.pgx.isec$type = paste0(snv.pgx.isec$V4,'>',snv.pgx.isec$V5)
snv.pgx.isec$VAF= str_split(snv.pgx.isec$V8,';',simplify = T)[,5] %>% gsub('VAF=','',.) %>% as.numeric()
snv.pgx.isec$Class = 'Common'
snv.pgx.isec$type = gsub(',(.*)','',snv.pgx.isec$type)
snv.pgx.isec.type <- as.data.frame(table(snv.pgx.isec$type) / length(snv.pgx.isec$type) * 100 )
colnames(snv.pgx.isec.type) <- c('SNV.Type','PGx')

snv.pgx.isec.type <- snv.pgx.isec.type[order(snv.pgx.isec.type$PGx), ]

p_snv_isec_type = ggplot(snv.pgx.isec.type, aes(x = SNV.Type,y=PGx)) +
  labs(x = "SNV.Type", y = "Frequency (%)") + xlim(snv.pgx.isec.type$SNV.Type)+
  geom_bar(position = "dodge", stat = "identity", width = 0.5,fill = '#65B48E', color = "black")+
  theme_bw()+scale_fill_brewer(palette = "Set1")+scale_color_brewer(palette = "Set1")+
  my_theme

ggsave(p_snv_isec_type,filename = 'PGx.vs.SEQC2/SNVs.type.Isec.pdf')
df.pgx.seqc2 <- rbind(snv.pgx[,c('V1','VAF','Class')],snv.seqc2[,c('V1','VAF','Class')]) %>% rbind(.,snv.pgx.isec[,c('V1','VAF','Class')])
head(df.pgx.seqc2)
colnames(df.pgx.seqc2) <- c('Chr','VAF','Class')

df.pgx.seqc2.vaf.p <- ggplot() +
  geom_histogram(data=df.pgx.seqc2,aes(VAF,fill=Class),color="white",bins = 50)+
  theme_bw()+
  theme_nature_border()+
  xlab(paste('VAF'))+ylab(paste('Counts'))+
  scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, by = 0.1))+scale_fill_brewer(palette ="Set1")


# 百分比堆积条形图

df.pgx.seqc2.percent <- df.pgx.seqc2 %>%
  dplyr::group_by(Class,Chr) %>%
  dplyr::summarise(Frequency = dplyr::n()) %>%
  dplyr::group_by(Class) %>%
  dplyr::mutate(Frequency = Frequency / sum(Frequency))

df.pgx.seqc2.chr.p <- ggplot(df.pgx.seqc2.percent, aes(x = Chr,y=Frequency, fill = Class)) + geom_bar(position = "dodge", stat = "identity")+
  scale_fill_brewer(palette ="Set1")+theme_bw()+
  theme_nature_border()+xlim(paste0('chr',c(1:22)))

ggsave("PGx.vs.SEQC2/pgx.seqc2.VAF.pdf",df.pgx.seqc2.vaf.p,width = 5.4*2,height = 4.5)
ggsave("PGx.vs.SEQC2/pgx.seqc2.Chr.pdf",df.pgx.seqc2.chr.p,width = 5.4*2,height = 4.5)

#SEQC2 2023 与PGx 2023的差别----------
pgx.2023 <- read.table("PGx.vs.SEQC2/PGx_2023.only.Indel.SNV.vcf",sep='\t')
pgx.2023$VAF= str_split(pgx.2023$V8,';',simplify = T)[,5] %>% gsub('VAF=','',.) %>% as.numeric()
pgx.2023$Class = 'PGx.2023'
seqc2.2023 <- read.table("PGx.vs.SEQC2/SEQC2_2023.only.Indel.SNV.vcf",sep='\t')
seqc2.2023$VAF= str_split(seqc2.2023$V8,';',simplify = T)[,5] %>% gsub('VAF=','',.) %>% as.numeric()
seqc2.2023$Class = 'SEQC2.2023'
pgx.common.seqc2.2023 <- read.table("PGx.vs.SEQC2/SEQC_common_PGx_2023.Indel.SNV.vcf",sep='\t')
pgx.common.seqc2.2023$VAF= str_split(pgx.common.seqc2.2023$V8,';',simplify = T)[,5] %>% gsub('VAF=','',.) %>% as.numeric()
pgx.common.seqc2.2023$Class = 'Common.2023'
df.2023 <- rbind(pgx.2023[,c('V1','VAF','Class')],seqc2.2023[,c('V1','VAF','Class')]) %>% rbind(.,pgx.common.seqc2.2023[,c('V1','VAF','Class')])
head(df.2023)
colnames(df.2023) <- c('Chr','VAF','Class')

df.2023.vaf.p <- ggplot() +
  geom_histogram(data=df.2023,aes(VAF,fill=Class),color="white",bins = 50)+
  my_theme+
  xlab(paste('VAF'))+ylab(paste('Counts'))+
  scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, by = 0.1))+scale_fill_brewer(palette ="Set1")

# 百分比堆积条形图

df.2023.percent <- df.2023 %>%
  dplyr::group_by(Class,Chr) %>%
  dplyr::summarise(Frequency = dplyr::n()) %>%
  dplyr::group_by(Class) %>%
  dplyr::mutate(Frequency = Frequency / sum(Frequency))

df.2023.chr.p <- ggplot(df.2023.percent, aes(x = Chr,y=Frequency, fill = Class)) + geom_bar(position = "dodge", stat = "identity")+
  scale_fill_brewer(palette ="Set1")+
  my_theme+xlim(paste0('chr',c(1:22)))

ggsave("PGx.vs.SEQC2/2023.VAF.pdf",df.2023.vaf.p,width = 5.4*2,height = 4.5)
ggsave("PGx.vs.SEQC2/2023.Chr.pdf",df.2023.chr.p,width = 5.4*2,height = 4.5)

##----交集VAF相关性------------------
##SNV
get_Isec_VAF <- function(seqc.vcf,pgx.vcf,var.type){
  snv.isec.pgx <- read.table(pgx.vcf)
  snv.isec.pgx$tag <- paste(snv.isec.pgx$V1,snv.isec.pgx$V2,snv.isec.pgx$V4,snv.isec.pgx$V5,sep = '&')
  snv.isec.pgx$VAF.PGx <- lapply(snv.isec.pgx$V8,function(x){
    tvaf = grep("^VAF=", unlist(strsplit(x, ";")), value = TRUE)
    return(gsub('VAF=','',tvaf))
  }) %>% unlist()
  snv.isec.pgx$Class = gsub('PASS;','',snv.isec.pgx$V7)
  
  snv.isec.seqc2 <- read.table(seqc.vcf)
  
  snv.isec.seqc2$VAF.SEQC2 <- lapply(snv.isec.seqc2$V8,function(x){
    tvaf = grep("^TVAF=", unlist(strsplit(x, ";")), value = TRUE)
    return(gsub('TVAF=','',tvaf))
  }) %>% unlist()
  
  snv.isec.seqc2$tag <- paste(snv.isec.seqc2$V1,snv.isec.seqc2$V2,snv.isec.seqc2$V4,snv.isec.seqc2$V5,sep = '&')
  
  snv.isec <- merge(snv.isec.pgx[,c('tag','VAF.PGx','Class')],
                    snv.isec.seqc2[,c('tag','VAF.SEQC2')],by='tag',all=TRUE)
  
  snv.isec$VAF.PGx <- as.numeric(snv.isec$VAF.PGx)
  snv.isec$VAF.SEQC2 <- as.numeric(snv.isec$VAF.SEQC2)
  snv.isec$Type=var.type
  return(snv.isec)
}

snv_indel.isec.df <- get_Isec_VAF("PGx.vs.SEQC2/SEQC2_commonwith_PGx.vcf","PGx.vs.SEQC2/PGx_commonWith_SEQC2.vcf",'All')

#snv.isec.df <- get_Isec_VAF("PGx.vs.SEQC2/SNV/Isec.in.SEQC2.SNV.vcf","PGx.vs.SEQC2/SNV/Isec.in.PGx.SNV.vcf",'SNV')
#indel.isec.df <- get_Isec_VAF("PGx.vs.SEQC2/Indel/Isec.in.SEQC2.Indel.vcf","PGx.vs.SEQC2/Indel/Isec.in.PGx.Indel.vcf",'Indel')


snv.isec.p <- ggplot(snv_indel.isec.df, aes(x = VAF.SEQC2, y = VAF.PGx)) + theme_classic()+
  geom_point(aes(fill=Class),shape = 21, size = 1.5, color = "white")+ scale_fill_brewer(palette = "Dark2")+
  geom_smooth(method = "lm", se = FALSE, color = "black") +stat_cor(method = "pearson",size = 5) +
  my_theme+coord_fixed()+xlim(c(0,1))+ylim(c(0,1))+
  labs(x = "VAF (SEQC2)", y = "VAF (PGx)")

indel.isec.p <- ggplot(indel.isec.df, aes(x = VAF.SEQC2, y = VAF.PGx)) + theme_classic()+
  geom_point(aes(fill=Class),shape = 21, size = 1.5, color = "white")+ scale_fill_brewer(palette = "Dark2")+
  geom_smooth(method = "lm", se = FALSE, color = "black") +stat_cor(method = "pearson",size = 5) +
  my_theme+coord_fixed()+xlim(c(0,1))+ylim(c(0,1))+
  labs(x = "VAF (SEQC2)", y = "VAF (PGx)")


isec.p <- ggplot(snv_indel.isec.df, aes(x = VAF.SEQC2, y = VAF.PGx)) + 
  geom_point(aes(fill=Class),shape = 21, size = 2, color = "white")+ scale_fill_brewer(palette = "Dark2")+
  geom_smooth(method = "lm", se = FALSE, color = "black") +stat_cor(method = "pearson",size = 5) +
  my_theme+coord_fixed()+xlim(c(0,1))+ylim(c(0,1))+
  labs(x = "VAF (SEQC2)", y = "VAF (PGx)")

ggsave(isec.p,filename = 'PGx.vs.SEQC2/All.isec.VAF.pdf')

#PGx和SEQC2特有的变异位点的染色体分布-----
common_df <- read.table("PGx.vs.SEQC2/PGx_commonwith_SEQC2.vcf")
pgx.only.df <- read.table("PGx.vs.SEQC2/PGx_only.vcf")
seqc2.only.df <- read.table("PGx.vs.SEQC2/SEQC2_only.vcf")

chr.common <- common_df$V1 %>% table() %>% as.data.frame()
chr.common$Class = 'Common'
chr.common$Percent = (chr.common$Freq / sum(chr.common$Freq))*100 
chr.pgx <- pgx.only.df$V1 %>% table() %>% as.data.frame()
chr.pgx$Class = 'PGx'
chr.pgx$Percent = (chr.pgx$Freq / sum(chr.pgx$Freq))*100 
chr.seqc2 <- seqc2.only.df$V1 %>% table() %>% as.data.frame()
chr.seqc2$Class = 'SEQC2'
chr.seqc2$Percent = (chr.seqc2$Freq / sum(chr.seqc2$Freq))*100 

chr.df <- rbind(chr.common,chr.pgx) %>% rbind(.,chr.seqc2)
colnames(chr.df) <- c('Chr','Freq','Class','Percent')



ggplot(data = chr.df,aes(x=Chr,y=Percent,fill=Class))+
  geom_bar(stat = "identity",
           # "fill"
           position = "stack")

ggplot(chr.df, aes(x = Chr, y = Freq, fill = Class)) +
  geom_bar(stat = "identity", width = 1) +
  #coord_polar("x", start = 0) +
  theme_void() +
  labs(title = "Percentage Bar Chart", y = "Percentage")


#F1 and Reproducibility information------------------
f1_df <- read.csv('F1score_stats_total.csv')
f1_df$source <- gsub('WGS_WUX','WUX_WGS',f1_df$source)
f1_df <- subset(f1_df,!grepl('2021',f1_df$source))
f1_df$Batch <- str_split(f1_df$source,'\\_WGS',simplify = T)[,1]
f1.mean_sd <- f1_df %>% dplyr::group_by(type,Batch) %>% dplyr::summarise(
  mean_f1 = mean(F1.score,na.rm=TRUE),
  sd_f1 = sd(F1.score,na.rm=TRUE)
)
f1.mean_sd$type <- gsub('indels','INDELs',f1.mean_sd$type)
f1.mean_sd <- subset(f1.mean_sd,f1.mean_sd$type !='records')

rep_df <- read.csv('TNseq_VCF_PASS_JaccardIndex.csv')
rep_df$source <- gsub('WGS_WUX','WUX_WGS',rep_df$source)
rep_df$smp1 <- str_split(rep_df$source,'VS',simplify = T)[,1]
rep_df$batch1 <- str_split(rep_df$smp1,'_WGS',simplify = T)[,1]
rep_df$smp2 <- str_split(rep_df$source,'VS',simplify = T)[,2]
rep_df$batch2 <- str_split(rep_df$smp2,'_WGS',simplify = T)[,1]
rep_df <- subset(rep_df,batch1==batch2)
rep_df$Batch <- rep_df$batch1
rep.mean_sd <- rep_df %>% group_by(type,Batch) %>%
  dplyr::summarise(
    mean_ji = mean(JI,na.rm=TRUE),
    sd_ji = sd(JI,na.rm=TRUE)
  )
rep.mean_sd$type <- gsub('indels','INDELs',rep.mean_sd$type) %>% gsub('snps','SNVs',.)

f1_rep.mean_sd <- merge(f1.mean_sd,rep.mean_sd,by=c('type','Batch'),all.x=T)
f1_rep.mean_sd <- subset(f1_rep.mean_sd,!grepl('merge_test',f1_rep.mean_sd$Batch))
f1_rep.mean_sd[,-c(1,2)] <- f1_rep.mean_sd[,-c(1,2)] %>% round(.,3)
write.csv(f1_rep.mean_sd,'SNV_INDEL.performance.info.csv',row.names = F)




#标准参考数据集的注释信息-----
indel_anno <- read.csv('PGx/high-confidence_sIndel_v2.HCR.clean.sort.hg38_multianno.txt',sep = '\t')
#coding
indel_anno[grepl('lung|intestine',indel_anno$cosmic95_coding),]


#noncoding
indel_anno[grepl('lung|intestine',indel_anno$cosmic95_noncoding),] 

snv_anno <- read.csv('PGx/high-confidence_sSNV_v2.HCR.clean.sort.hg38_multianno.txt',sep = '\t')
#coding
snv_anno[grepl('lung|intestine',snv_anno$cosmic95_coding),]

#noncoding
snv_anno[grepl('lung|intestine',snv_anno$cosmic95_noncoding),] 

aim_var.coding <- rbind(indel_anno[grepl('lung|intestine',indel_anno$cosmic95_coding),],
      snv_anno[grepl('lung|intestine',snv_anno$cosmic95_coding),])

write.csv(aim_var.coding,'PGx/lung_intestine_var.coding.csv',row.names = F)

aim_var.uncoding <- rbind(indel_anno[grepl('lung|intestine',indel_anno$cosmic95_noncoding),],
                         snv_anno[grepl('lung|intestine',snv_anno$cosmic95_noncoding),])

write.csv(aim_var.uncoding,'PGx/lung_intestine_var.uncoding.csv',row.names = F)

