library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggsignif)
library(patchwork)
library(ggpubr)
library(cowplot)
library(ggsci)
options(scipen=200)
source('/Users/duanshumeng/生物信息/PGx_lab/毕业论文/Scripts/my_theme.R')
setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/SEQC2_analysis/hcc1395_analysis_pipeline/Pipeline_validation')



#统计每个caller识别到的PASS的变异位点数,可视化-----------------------
vcf_count = read.csv("VCF_counts.vcf_stats.csv")
vcf_count_hcr = read.csv("VCF_counts.HCR.vcf_stats.csv")

vcf_count = vcf_count[,c('type','count','source')]
colnames(vcf_count) <- c('type','Total_count','source')

vcf_count_hcr = vcf_count_hcr[,c('type','count','source')]
colnames(vcf_count_hcr) <- c('type','HCR_count','source')
vcf_merge <- merge(vcf_count,vcf_count_hcr,by=c('type','source'))
vcf_merge$source <- gsub('_stralka','.strelka',vcf_merge$source)
vcf_merge['caller'] <- str_split(vcf_merge$source,'\\.',simplify = TRUE)[,2]
vcf_merge$source <- paste0(str_split(vcf_merge$source,'\\.',simplify = TRUE)[,2],'.',str_split(vcf_merge$source,'\\.',simplify = TRUE)[,1])

for (i in (vcf_merge$type %>% unique())){
  print(i)
  vcf_sub = subset(vcf_merge,type==i)
  data = melt(vcf_sub)
  colnames(data) <- c('type','source','Caller','Count_type','Number of variants')
  class(data$`Number of variants`)
  
  #combn(unique(data$Caller),2,simplify=FALSE)
  
  p1 <- ggplot(data,aes(source,`Number of variants`,fill=Caller,alpha=Count_type))+
    geom_bar(stat = "identity",color="black",position = position_dodge(),width = 0.7,size=0.25)+theme_classic()+
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          axis.text.x = element_text(angle = 45,hjust = 1),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15))+labs(y=paste0('Number of ',i),x='')+scale_fill_brewer(palette = "Accent")
  
  ggsave(p1,filename = paste0(i,'_count_barplot.png'))
  ggsave(p1,filename = paste0(i,'_count_barplot.pdf'))
  
}



#SEQC2原始vcf中PASS的count数
vcf_count_S = read.csv("VCF_counts.vcf_stats_SEQC2.csv")
vcf_count_S <- subset(vcf_count_S,grepl('bwa',vcf_count_S$source) & grepl('WGS_EA_1|WGS_FD_1|WGS_NS_1',vcf_count_S$source))
vcf_count_S$source <- gsub('.bwa','',vcf_count_S$source) %>% gsub('_1','',.) %>% gsub('PASS.vcf.gz','.SEQC2',.)
vcf_count_S <- subset(vcf_count_S,count != 0)

#vcf_count_S = subset(vcf_count_S,count != 0)

vcf_count_hcr_S = read.csv("VCF_counts.HCR.vcf_stats_SEQC2.csv")
vcf_count_hcr_S <- subset(vcf_count_hcr_S,grepl('bwa',vcf_count_hcr_S$source) & grepl('WGS_EA_1|WGS_FD_1|WGS_NS_1',vcf_count_hcr_S$source))
vcf_count_hcr_S$source <- gsub('.bwa','',vcf_count_hcr_S$source) %>% gsub('_1','',.) %>% gsub('PASS.vcf.gz','.SEQC2',.)
vcf_count_hcr_S <- subset(vcf_count_hcr_S,count != 0)

vcf_count = read.csv("VCF_counts.vcf_stats.csv")
vcf_count$source <- gsub('_stralka','.strelka',vcf_count$source) %>% gsub('strelka','strelka',.) %>% gsub('.PASS.vcf.gz','.PGx',.)
vcf_count_hcr = read.csv("VCF_counts.HCR.vcf_stats.csv")
vcf_count_hcr$source <- gsub('_stralka','.strelka',vcf_count_hcr$source) %>% gsub('strelka','strelka',.) %>% gsub('.PASS.vcf.gz','.PGx',.)

vcf_count = vcf_count[,c('type','count','source')]
vcf_count_hcr = vcf_count_hcr[,c('type','count','source')]

vcf_count
compare_count <- function(vcf,vcf_s,region){
  vcf['batch'] = 'PGx'
  vcf_s['batch'] = 'SEQC2'
  vcf_compare = rbind(vcf[,c('type','count','batch','source')],
                      vcf_s[,c('type','count','batch','source')])
  vcf_compare['Caller'] <- str_split(vcf_compare$source,'\\.',simplify = T)[,2]
  vcf_compare$source <- paste0(str_split(vcf_compare$source,'\\.',simplify = TRUE)[,2],'.',
                              str_split(vcf_compare$source,'\\.',simplify = TRUE)[,1],',',
                              str_split(vcf_compare$source,'\\.',simplify = TRUE)[,3])
  bar_p <- ggplot(data=vcf_compare, aes(x=source, y=count, fill=Caller)) +theme_classic()+
    geom_bar(stat="identity",color="black")+
    scale_fill_brewer(palette ='Accent')+
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 60,hjust = 1),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 20),
          strip.text = element_text(size = 20))+xlab('')+
    theme(panel.spacing.x = unit(1, "mm"))+facet_grid(rows = vars(type))
  ggsave(bar_p,filename = paste0('Count_barplot_',region,'.png'),width = 14)
  ggsave(bar_p,filename = paste0('Count_barplot_',region,'.pdf'),width = 14)
  
  boxdot_plot <- function(data,tag){
    #color_set = ifelse(tag == 'Whole_Genome','Pastel1','Pastel2')
    #data = subset(data_row,type==snv_type)
    comp_list = list(c('SEQC2','PGx'))
    
    p1 = ggplot(data,aes(x=batch,y=count,fill=batch))+theme_classic()+
      geom_boxplot(notch=FALSE,alpha=0.5,)+
      geom_point(alpha=0.5,size = 1,colour = "black")+geom_signif(comparisons = comp_list,test.args =c(exact=FALSE),
                                                                  map_signif_level = function(p) sprintf("p = %0.2f", p),textsize = 7)+
      theme(axis.title = element_text(size = 25),
            axis.text = element_text(size = 20),
            legend.text = element_text(size = 20),
            legend.title = element_text(size = 25),
            strip.text = element_text(size = 20))+
      labs(x='',y='Count')+
      scale_fill_brewer(palette ='Pastel1')+facet_grid(cols = vars(type))
    return(p1)
  }
  
  p_sl = boxdot_plot(vcf_compare,region)
  
  ggsave(p_sl,filename = paste0('Count_compare_',region,'.png'),width = 14)
  ggsave(p_sl,filename = paste0('Count_compare_',region,'.pdf'),width = 14)
  
  
}

#vcf_count vs. vcf_count_S
compare_count(vcf_count,vcf_count_S,'Whole_Genome')
compare_count(vcf_count_hcr,vcf_count_hcr_S,'High_Confidence_Region')

#比较VAF的分布------------------------
vcf_list = dir(pattern = '.VAF.csv')

P_plot <- function(vcf_df,tag){
  p = ggplot(data=vcf_df,aes(x=VAF))+
    geom_histogram(aes(fill=Type),
                   bins = 50)+theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = c(0.8,0.8))+
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15))+
    scale_x_continuous(expand = expansion(mult = c(0,0)),
                       breaks = seq(0,1,0.2))+
    scale_y_continuous(expand = expansion(mult = c(0,0)))+scale_fill_brewer(palette = "Accent")+ggtitle(tag)+ylim(0,40000)
  return(p)
  
}
vcf_df_1 = read.csv(vcf_list[1])
p_raw = ggdraw(P_plot(vcf_df_1,gsub('.PASS.VAF.csv','',vcf_list[1])))

for (i in vcf_list[-1]){
    print(i)
    vcf_df <- read.csv(i)
    p = ggdraw(P_plot(vcf_df,gsub('.PASS.VAF.csv','',i)))
    p_raw = p_raw + p 

}
ggsave(p_raw,filename = paste0('VAF_count_histogram.png'),width = 10)
ggsave(p_raw,filename = paste0('VAF_count_histogram.pdf'),width = 10)

#找出共有的变异位点：
#1. 样本内：至少有两个caller鉴定到的位点，作为该样本的可鉴定变异位点
#2. 样本间：至少在两个样本中被检测到的位点
smps = str_split(vcf_list,'\\.',simplify = T)[,1] %>% unique()
total_lst = list()

for (i in smps){
  tag_lst = c()
  for (a in vcf_list[grepl(i,vcf_list)]){
    vcf_df <- read.csv(a)
    tag_lst = append(tag_lst,vcf_df$tag)
    tag_df = table(tag_lst) %>% as.data.frame()
    tag_df = tag_df[tag_df$Freq > 2,]
  }
  tag_lst[[i]] = tag_df
}

#SEQC2高置信数据集的VAF
library(data.table)
seqc2_vaf <- fread("SEQC2_HRC_VAF.txt",sep='\t',stringsAsFactors = FALSE) %>% as.data.frame()
seqc2_vaf$VAF <- gsub('TVAF=','',seqc2_vaf$TVAF) %>% as.numeric()
seqc2_vaf$VAF95 <- str_split(seqc2_vaf$TVAF95,'\\=',simplify = TRUE)[,2]
p = ggplot(data=seqc2_vaf,aes(x=VAF))+
  geom_histogram(fill="#88ada6", color="#fffbf0",alpha=.5,bins = 50)+theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8))+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+
  scale_x_continuous(expand = expansion(mult = c(0,0)),
                     breaks = seq(0,1,0.2))+
  scale_y_continuous(expand = expansion(mult = c(0,0)))+scale_fill_brewer(palette = "Accent")

VAF_list = c()
nm_list = c()
for (i in vcf_list){
  print(i)
  vcf_df <- read.csv(i)
  nm = str_split(i,'\\.',simplify = TRUE)[,2]
  VAF_list= append(VAF_list,vcf_df[vcf_df$Type=='HRC','VAF'])
  nm_list = append(nm_list,rep(nm, length(vcf_df[vcf_df$Type=='HRC','VAF'])))
}

VAF_list = append(VAF_list,seqc2_vaf$VAF)
nm_list = append(nm_list,rep('SEQC2',length(seqc2_vaf$VAF)))

VAF_df <- data.frame(VAF=VAF_list,Type=nm_list)

p_seqc2 = P_plot(VAF_df,'SEQC2_HRC')
ggsave(p_seqc2,filename = paste0('VAF_compre_histogram.png'),width = 10)
ggsave(p_seqc2,filename = paste0('VAF_compre_histogram.pdf'),width = 10)

#https://www.jianshu.com/p/bf5ca7f3c5c5
violin_plot_F1 <- function(data){
  comp_list = combn(data$Type %>% unique(),2,simplify = FALSE)
  #comp_list = list(c("TNscope","SEQC2"))
  p1 = ggplot(data,aes(x=Type,y=VAF,fill=Type))+theme_classic()+
    geom_violin(cex=1.2)+
    geom_boxplot(width=0.1,cex=1.2)+
    geom_signif(comparisons = comp_list,test.args =c(exact=FALSE),y_position = c(1.00, 1.04, 1.08),
                                                                map_signif_level = function(p) sprintf("p = %0.2f", p),textsize = 5)+
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15))+labs(x='',y='Variant allele frequency')+
    scale_fill_brewer(palette = "Accent")+guides(fill=FALSE)+ylim(0,1.15)
  return(p1)
}


ggsave(violin_plot_F1(VAF_df),filename = paste0('VAF_compre_violin.png'),width = 10)
ggsave(violin_plot_F1(VAF_df),filename = paste0('VAF_compre_violin.pdf'),width = 10)

#与SEQC2比较的F1 score------------------------
f1_s = read.csv("stats_total_compares_SEQC2.csv")
f1_s$source <- gsub('.bwa','',f1_s$source) %>% gsub('_1','',.) %>% paste0(.,'.SEQC2')
#f1 = read.csv("F1score_stats_total.csv")
f1 = read.csv("SEQC2_F1score_stats_total.validation.csv")
f1$source <- gsub('_stralka','.strelka',f1$source) %>% 
  gsub('strelka','strelka',.) %>% 
  gsub('.PASS','.PGx',.) %>% gsub('_deepsomatic','.deepsomatic',.) %>% gsub('_mutect2','.mutect2',.)

#f1$source <- str_split(f1$source,'\\.',simplify = TRUE)[,2]

f1_s['F1.score'] <- 2 * (f1_s['precision'] * f1_s['recall']) / (f1_s['precision'] + f1_s['recall'])
#f1['F1.score'] <- 2 * (f1['precision'] * f1['recall']) / (f1['precision'] + f1['recall'])

#合并SEQC2和PGx的数据绘制散点图
f1_s.sub <- f1_s[,c('type','source','F1.score')]
colnames(f1_s.sub) <- c('type','source','SEQC2.F1.score')
f1_s.sub$Sample <- str_split(f1_s.sub$source,'\\.',simplify = T)[,1]
f1_s.sub$Caller <- str_split(f1_s.sub$source,'\\.',simplify = T)[,2]
subset(f1_s.sub,f1_s.sub$type!='records',c("type","SEQC2.F1.score","Sample","Caller"))

#f1_s.sub$source <- sub('.strelka.SEQC2','_1',f1_s.sub$source) %>% sub('.muTect2.SEQC2','_2',.)  %>% sub('.somaticSniper.SEQC2','_3',.) 
f1.sub <- f1[,c('type','source','F1.score')] 
colnames(f1.sub) <- c('type','source','PGx.F1.score')
f1.sub$Sample <- str_split(f1.sub$source,'\\.',simplify = T)[,1]
f1.sub$Caller <- str_split(f1.sub$source,'\\.',simplify = T)[,2]
subset(f1.sub,f1.sub$type!='records',c("type","PGx.F1.score","Sample","Caller"))
#f1.sub$source <- sub('.strelka.PGx','_1',f1.sub$source) %>% sub('.TNseq.PGx','_2',.) %>% sub('.TNscope.PGx','_3',.)
f1.merge <- merge(f1_s.sub,f1.sub,by=c('type','source'))

F1.cor.snv <- ggplot(subset(f1.merge,type=='SNVs'), aes(x = SEQC2.F1.score, y = PGx.F1.score)) + theme_classic()+
  geom_point(size=5,alpha=0.75,color='#8E1025')+
  geom_smooth(method = "lm", se = FALSE, color = "black") +stat_cor(method = "pearson",size = 5)+
  my_theme+
  labs(x = "F1.score (SEQC2)", y = "F1.score (PGx)");F1.cor.snv


F1.cor.indel <- ggplot(subset(f1.merge,type=='indels'), aes(x = SEQC2.F1.score, y = PGx.F1.score)) + theme_classic()+
  geom_point(size=5,alpha=0.75,color='#8E1025')+
  geom_smooth(method = "lm", se = FALSE, color = "black") +stat_cor(method = "pearson",size = 5)+
  my_theme+
  labs(x = "F1.score (SEQC2)", y = "F1.score (PGx)");F1.cor.indel

F1.cor.p <- F1.cor.indel+F1.cor.snv+plot_annotation(tag_levels = 'A')+
  plot_layout(guides='collect') & theme(legend.position='',
                                        plot.tag = element_text(color = "black",size = 18,face = "bold"))

ggsave('F1_score_compare_point.pdf',F1.cor.p,width=8.15, height=5.7*0.7,dpi = 300)

#barplot
f1_combine = rbind(f1,f1_s)
f1_combine = subset(f1_combine,type != 'records')
f1_combine['Caller'] = str_split(f1_combine$source,'\\.',simplify = T)[,2]
f1_combine['Platform'] = str_split(f1_combine$source,'\\.',simplify = T)[,3]
f1_combine$type = gsub('indels','Indels',f1_combine$type)
f1_combine$source <- paste0(str_split(f1_combine$source,'\\.',simplify = TRUE)[,2],'.',
                            str_split(f1_combine$source,'\\.',simplify = TRUE)[,1],',',
                            str_split(f1_combine$source,'\\.',simplify = TRUE)[,3])

h1 <- ggplot(data=f1_combine, aes(x=source, y=F1.score, fill=Caller))+theme_classic() +
  geom_bar(stat="identity",color="black")+
  scale_fill_brewer(palette ='Accent')+
  geom_hline(aes(yintercept=0.85), colour="Red", linetype="dashed")+
  geom_hline(aes(yintercept=0.70), colour="Blue", linetype="dashed")+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 60,hjust = 1),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 20))+xlab('')+
  theme(panel.spacing.x = unit(1, "mm"))+facet_grid(rows = vars(type))

ggsave(paste0('F1_score_compare_bar.png'),plot = h1,width = 10)
ggsave(paste0('F1_score_compare_bar.pdf'),plot = h1,width = 10)
#scatter_plot

scatter_plot <- function(data,tag){
  data$recall <- as.numeric(data$recall)
  data$precision <- as.numeric(data$precision)
  data$source <- as.factor(data$source)
  p1 = ggplot(data,aes(x=precision,y=recall,colour=source,
                       shape=type,fill=source))+theme_classic()+
    geom_point(size=5,alpha = 0.5)+scale_color_igv()+scale_fill_igv()+
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 20))+
    labs(x = "Precision", y = "Recall")+ylim(0,1)+xlim(0,1)+ggtitle(paste0(tag,' pipeline'))
  ggsave(paste0(tag,'_compare_scatter.png'),plot = p1)
  ggsave(paste0(tag,'_compare_scatter.pdf'),plot = p1)
}
f1_s = subset(f1_s,type != 'records')
f1 = subset(f1,type != 'records')
scatter_plot(f1_s,'SEQC2')
scatter_plot(f1,'PGx')

#比较SEQC2与PGx的流程计算出的变异的可重复性----------------
rep.df <- read.csv('PGx.vs.SEQC2.reproducibility.csv')
library(ComplexHeatmap)
library(circlize)
#Reproducibility-------------------
rep.df$type <- gsub('indels','Indels',rep.df$type) %>% gsub('snps','SNVs',.)
add_batch <- function(df){
  df['batch'] <- ''
  df[(str_count(df$source,'bwa')==2),'Batch'] = 'SEQC2'
  df[(str_count(df$source,'bwa')==0),'Batch'] = 'PGx'
  df[(str_count(df$source,'bwa')==1),'Batch'] = 'PGx_vs_SEQC2'
  return(df)
}


rep.df <- add_batch(rep.df)
rep.df$source <- sub('_HCR','',rep.df$source) %>% sub('_1','',.) %>% sub('_strelka','',.)


#JI_df_T20 <- subset(JI_df,Platform=='BGI_T20')
#JI_df <- subset(JI_df,Platform %in% c('BGI_T20','WUX','BGI_T10'))
make_matrix <- function(df,snv_type,M){
  #M = matrix(1,6,6)
  df_sub <- subset(df, type==snv_type)
  s1 <- stringr::str_split(df_sub$source,'VS',simplify = T)[,1] %>% unique() 
  s2 <- stringr::str_split(df_sub$source,'VS',simplify = T)[,2] %>% unique() 
  smp_df <- c(s1,s2) %>% unique() %>% sort() 
  print(smp_df)
  colnames(M) <- smp_df
  rownames(M) <- smp_df
  for (i in 1:dim(df_sub)[1]){
    smp <- df_sub[i,'source']
    smp1 = stringr::str_split(smp,'VS',simplify = T)[,1] 
    smp2 = stringr::str_split(smp,'VS',simplify = T)[,2]
    ji = df_sub[i,'Reproducibility']
    M[smp1,smp2] = ji  
    M[smp2,smp1] = ji
  }
  return(M)
}

library(RColorBrewer)

rep_mat_indel <- make_matrix(rep.df,'Indels',matrix(1,15,15))
rep_mat_snvs <- make_matrix(rep.df,'SNVs',matrix(1,15,15))

anno_df <- data.frame(source=colnames(rep_mat_snvs),
                      Sample=str_split(colnames(rep_mat_snvs),'\\.',simplify = T)[,1],
                      Batch=ifelse(grepl('bwa',colnames(rep_mat_snvs)),'SEQC2','PGx'))

Sample_col = brewer.pal(n =length(unique(anno_df$Sample)) , name = "Set1")[1:length(unique(anno_df$Sample))]
names(Sample_col) = unique(anno_df$Sample)


Batch_col = brewer.pal(n =3 , name = "Dark2")[1:length(unique(anno_df$Batch))]
names(Batch_col) = unique(anno_df$Batch)
ha = HeatmapAnnotation(df=anno_df[,-1],
                       col = list(Batch = Batch_col,
                                  Sample = Sample_col))

ji_col_fun = colorRamp2(c(0.5, 0.75, 1), c("#377EB8", "white", "#E41A1C"))
p_snvs <- Heatmap(as.matrix(rep_mat_snvs),show_row_names = FALSE,show_column_names = FALSE, 
                  top_annotation = ha, col = ji_col_fun,
                  show_column_dend = FALSE, column_title = "SNVs",name = "Reproducibility")
pdf("SEQC2.vs.PGx.reproducibility_snvs.pdf",width = 8.15, height = 8.15)
p_snvs
dev.off()

p_indels <- Heatmap(as.matrix(rep_mat_indel),show_row_names = FALSE,show_column_names = FALSE, 
                  top_annotation = ha, col = ji_col_fun,
                  show_column_dend = FALSE, column_title = "Indels",name = "Reproducibility")
pdf("SEQC2.vs.PGx.reproducibility_indels.pdf",width = 8.15, height = 5.7)
p_indels
dev.off()

#比较F1_score的一致性--------------------
f1_s_sub = f1_s[,c('type','F1.score')]
f1_s_sub['batch'] = 'SEQC2'

f1_sub = f1[,c('type','F1.score')]
f1_sub['batch'] = 'PGx'

f1_combine <- rbind(f1_s_sub,f1_sub)
f1_combine$type = gsub('indels','Indels',f1_combine$type)
boxdot_plot_F1 <- function(data){
  #data = subset(data_row,type==snv_type)
  #color_set = ifelse(snv_type == 'SNVs','Pastel1','Pastel2')
  comp_list = list(c('SEQC2','PGx'))
  
  p1 = ggplot(data,aes(x=batch,y=F1.score,fill=batch))+theme_classic()+
    geom_boxplot(notch=FALSE,alpha=0.5,)+
    geom_point(alpha=0.5,size = 1,colour = "black")+geom_signif(comparisons = comp_list,test.args =c(exact=FALSE),
                                                                map_signif_level = function(p) sprintf("p = %0.2f", p),textsize = 7)+
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          strip.text = element_text(size = 20))+labs(x='',y='F1.score')+
    scale_fill_brewer(palette = 'Pastel1')+guides(fill=FALSE)+ylim(0,1)+facet_grid(cols = vars(type))
  return(p1)
}

p_sl = boxdot_plot_F1(subset(f1_combine,type!='records'))
#p_s = boxdot_plot_F1(subset(f1_combine,type!='records'),'SNVs')
#p_l = boxdot_plot_F1(subset(f1_combine,type!='records'),'indels')

ggsave(paste0('F1_score_compare_boxplot.png'),plot = p_sl,width = 10)
ggsave(paste0('F1_score_compare_boxplot.pdf'),plot = p_sl,width = 10)

##比较一个caller，两个caller组合，三个caller组合后，与SEQC2的一致性
call_single <- read.csv("F1score_stats_total.csv")
call_single$source <- gsub('_stralka','.strelka',call_single$source)
call_single['Caller'] <- str_split(call_single$source,'\\.',simplify = T)[,2]
call_isec <- read.csv("F1score_stats_total_caller_intersect.csv")
call_isec['Caller'] <- gsub('WGS_EA_|WGS_FD_|WGS_NS_|.isec|.PASS','',call_isec$source)

call_combine <- rbind(call_single,call_isec)
call_combine <- subset(call_combine,type!='records')
call_combine$type <- gsub('indels','Indels',call_combine$type)
call_combine['F1.score'] <- 2 * (call_combine['precision'] * call_combine['recall']) / (call_combine['precision'] + call_combine['recall'])
boxdot_plot <- function(data_row,snv_type){
  x_lims = c('TNscope','TNseq','strelka',"TNscope_strelka","TNseq_strelka","TNseq_TNscope","TNseq_TNscope_strelka")
  y_position = c(0.95,1.00,1.05,1.10,1.15,1.20)
  data = subset(data_row,type==snv_type)
  print(data$Caller)
  p1 = ggplot(data,aes(x=Caller,y=F1.score,fill=Caller))+theme_classic()+
    geom_boxplot(notch=FALSE,alpha=0.5,)+
    geom_point(alpha=0.5,size = 1,colour = "black")+
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          axis.text.x = element_text(angle = 60,hjust = 1),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15))+labs(x='',y='F1 score')+scale_fill_brewer(palette = "Pastel1")+ggtitle(paste0(snv_type))+
    ylim(c(0,1))+xlim(x_lims)+guides(fill=FALSE)
  return(p1)
  
}
sv_21 = boxdot_plot(call_combine,'SNVs')
in_21 = boxdot_plot(call_combine,'Indels')

ggsave(paste0('F1_score_call_combine_boxplot.png'),plot = sv_21+in_21,width = 14)
ggsave(paste0('F1_score_call_combine_boxplot.pdf'),plot = sv_21+in_21,width = 14)

