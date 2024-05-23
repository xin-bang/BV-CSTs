#参考：https://www.jianshu.com/p/31a92d6dc746
#参考：https://community.rstudio.com/t/vaginal-microbiome-cst/140076/2
#参考：https://www.nature.com/articles/s41522-022-00336-6#Sec7
#绘图可使用：https://biit.cs.ut.ee/clustvis/

#clustering_distance_rows:计算类间距的方式：默认：euclidean（欧式)
#clustering_method :cluster的方法：默认complete
#cutree_rows: 若对行进行了cluster,则可以指定cluster个数(基于hclust)
#kmeans_cols：the number of kmeans clusters to make, if we want to aggregate the rows before drawing heatmap 
#需要先使用Valencia对矩阵进行分型。
#! /home/luxingbang/miniconda3/envs/mamba/bin/Rscript
rm(list=ls())
set.seed(12)
suppressPackageStartupMessages({
  library(argparse)
  library(pheatmap)
  library(cowplot)
  library(ggplot2)
  library(grid)
  library(readxl)
  library(dplyr)
  library(ecodist)
})



# args = list(
#   f1 = "Valencia.CSTs.out.csv",
#   fc = 0.01
# )

# # #参数定义：
parser <- ArgumentParser(description="用BV-CSTs分型聚类，可以同时进行多种不同的聚类方法")
parser$add_argument("--f1", help="输入Valencia分析的矩阵。矩阵格式：横列为样本，纵列为菌群名，数值为read数目，具体格式见README")
parser$add_argument("--fc", help="过滤标准，样本中的菌群占比低于该标准的会被过滤掉，默认是0.01",default=0.01,type="numeric")
args <- parser$parse_args()     # 解析参数


df = read.csv(args$f1,header = TRUE)
#数据预处理：
df = df %>% select(-None)   #这里是否该去除None？
row.names(df) = df$sampleID
df = df %>% select(-sampleID)
df = df[, !grepl("_sim$", names(df))]
df_valencia = df %>% select(CST,subCST)
df_valencia$sampleID = rownames(df_valencia)
df = df %>% select(-c(subCST,CST,score))

#以valencia的分析结果确定其余几种分类方法的簇数
k = df_valencia$CST %>% unique() %>% length() 


##从原始矩阵中计算样本中每一种菌群的百分占比
df2 <- df %>%
  mutate_at(vars(-matches("read_count")), ~./read_count)
df2 = df2 %>% select(-read_count)
# print(rowSums(df2))


##剔除每一列均小于0.01的列
column_filter <- apply(df2, 2, function(col) all(col < args$fc))
data_cluster = df2[, !column_filter]



###使用bcdist计算距离
hclust_result <- hclust(bcdist(data_cluster),method = "ward.D2")
# k = args$k  # 设定簇数
cluster_labels <- cutree(hclust_result, k)
cluster_labels = as.data.frame(cluster_labels)
colnames(cluster_labels) = c("bcdist")
cluster_labels$sampleID <- rownames(cluster_labels)



###使用manhattan计算距离 :ward.D2
hclust_result2 = hclust(dist(data_cluster,method = "manhattan"),method = "ward.D2")
# k <- args$k  # 设定簇数
cluster_labels2 <- cutree(hclust_result2, k)
cluster_labels2 = as.data.frame(cluster_labels2)
colnames(cluster_labels2) = c("ward.D2")
cluster_labels2$sampleID <- rownames(cluster_labels2)



###使用manhattan计算距离 :complete
hclust_result3 = hclust(dist(data_cluster,method = "manhattan"),method = "complete")
# k <- args$k  # 设定簇数
cluster_labels3 <- cutree(hclust_result3, k)
cluster_labels3 = as.data.frame(cluster_labels3)
colnames(cluster_labels3) = c("complete")
cluster_labels3$sampleID <- rownames(cluster_labels3)


merged_df = cluster_labels %>% full_join(cluster_labels2,by="sampleID")
merged_df = merged_df %>% full_join(cluster_labels3,by = "sampleID")
merged_df = merged_df %>% full_join(df_valencia,by="sampleID")


row.names(merged_df) = merged_df$sampleID
merged_df = merged_df %>% select(-sampleID)
merged_df$CST <- factor(merged_df$CST)
merged_df$subCST <- factor(merged_df$subCST)
merged_df$bcdist <- factor(merged_df$bcdist)
merged_df$ward.D2 <- factor(merged_df$ward.D2)
merged_df$complete <- factor(merged_df$complete)


##查看数据类型
# str(data_cluster) #查看每一列  
# apply(data_cluster, 1, function(row) sapply(row, class))  #查看每一行

custom_color_range <- colorRampPalette(c("#f9ed69", "#f08a5d", "#b83b5e"))(100)

data_cluster = as.data.frame(t(data_cluster))
data_cluster[] <- lapply(data_cluster, as.numeric)

# annotation_colors <- list(
#   complete= c("1" = "#e9f1f6", "2" = "#ffa631", "3" = "#afdd22", "4" = "#ed5736", "5" = "#1685a9"),
#   ward.D2= c("1" = "#fff2df", "2" = "#ffc773", "3" = "#0eb83a", "4" = "#f9906f", "5" = "#3b2e7e"),
#   bcdist= c("1" = "#feff89", "2" = "#ff9f68", "3" = "#f85959", "4" = "#ff304f", "5" = "#7c203a"),
#   CST= c("CST1" = "#107a8b", "CST2" = "#2cb978", "CST3" = "#83e85a", "CST4" = "#3b5441", "CST5" = "#071a52"),
#   Valencia= c("CST1" = "#fff591", "CST2" = "#ff8a5c", "CST3" = "#f5587b", "CST4" = "#e41749", "CST5" = "red"),
# )

merged_df = merged_df %>% select(-subCST)
p5_annotation <- pheatmap(
  data_cluster,
  clustering_method = "ward.D2",
  clustering_distance_cols = "manhattan",
  cutree_cols = k,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = merged_df,
  color = custom_color_range
  # annotation_colors = annotation_colors  # 设置 annotation 颜色
)

write.csv(merged_df,file = "一致性.csv")
ggsave ("BV-CSTs.pdf", plot=p5_annotation,device = cairo_pdf,width = 380,height = 280, units = "mm")

###结论：
#bcdist : 使用bcdist计算距离（参考文章），聚类方法使用ward.D2
#ward.D2 : 使用manhattan计算距离（参考文章），聚类方法使用 ward.D2
#complete : 使用 manhattan 计算距离，聚类方法使用 complete
#输入数据使用 病原频率表




