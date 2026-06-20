library(dplyr)
source('/Users/duanshumeng/生物信息/PGx_lab/毕业论文/Scripts/my_theme.R')
setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/SEQC2_analysis/hcc1395_analysis_pipeline/High_confidence_SNVs&Indel/Validated_by_RNA')

dna_df <- read.table("sSNV_sIndel_v2.in_rna_regions.validates.vcf")

dna_df$VAF <- str_split(dna_df$V8,'\\;',simplify = T)[,5] %>% gsub('VAF=','',.) %>% as.numeric()
dna_df$Tag <- paste0(dna_df$V1,':',dna_df$V2,':',dna_df$V4,':',dna_df$V5)

rna_df <- read.table("HCC1395.RNA.isec.validates.vcf")
library(tidyverse)

# 1. 假设你的 dataframe 名称为 df
# 为了演示，先处理 V10 到 V15 这几列
df_vaf <- rna_df %>%
  mutate(across(starts_with("V"), function(x) {
    # 步骤 A: 处理缺失值标识 ".:.:.:.:."，将其统一为 NA
    x[x == ".:.:.:.:."] <- NA
    
    # 步骤 B: 以 ":" 拆分字符串
    # 我们需要第二个元素 (AD) 和 第三个元素 (DP)
    parts <- str_split_fixed(x, ":", 5)
    ad_full <- parts[, 2]
    dp <- as.numeric(parts[, 3])
    
    # 步骤 C: 对 AD 以 "," 拆分，获取 ALT (第二个元素)
    # str_split_fixed 处理 "REF,ALT" 结构
    ad_parts <- str_split_fixed(ad_full, ",", 2)
    alt <- as.numeric(ad_parts[, 2])
    
    # 步骤 D: 计算 VAF (ALT / DP)
    vaf <- alt / dp
    
    return(vaf)
  }))

# 查看处理后的结果
head(df_vaf)

df_vaf$VAF <- rowMeans(df_vaf,na.rm = TRUE)

rna_df$VAF <- df_vaf$VAF

rna_df$Tag <- paste0(rna_df$V1,':',rna_df$V2,':',rna_df$V4,':',rna_df$V5)

rna_df[,c('Tag','VAF')] %>% head()

dna_df[,c('Tag','VAF')] %>% head()

library(tidyverse)
library(ggpubr)
library(ggExtra)
plot_df <- inner_join(dna_df[,c('Tag','VAF')], rna_df[,c('Tag','VAF')], by = "Tag", suffix = c("_DNA", "_RNA"))

p <- ggplot(plot_df, aes(x = VAF_DNA, y = VAF_RNA)) +
  # 绘制散点
  geom_point(alpha = 0.8, size = 1.5,shape = 21,fill="#1F77B4",
   color = "white",
   stroke = 0.5)+
  # 添加 1:1 理想拟合线 (虚线)
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  # 添加线性回归拟合线 (蓝色实线)
  geom_smooth(method = "lm", color = "#E15759", fill = "lightgray") +
  # 设置坐标轴范围 (VAF 0到1)
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  # 添加标签和主题
  labs(x = "VAF (DNA)", y = "VAF (RNA)")+
  # 在图中添加相关系数 (R) 
  stat_cor(method = "pearson", label.x = 0.1, label.y = 0.9)+
  my_theme

# 将散点图转换为带边际分布的组合图
final_plot <- ggExtra::ggMarginal(p, 
                                  type = "histogram",    # 也可以选 "density"
                                  fill = "#1F77B4", 
                                  alpha = 0.6,
                                  color = "white",
                                  margins = "both")

# 显示图片


ggsave('DNA_RNA_VAF.cor.pdf',final_plot,width = 5.7*0.8,height = 5.7*0.8)
