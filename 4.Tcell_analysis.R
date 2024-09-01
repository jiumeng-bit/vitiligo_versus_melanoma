library(Seurat)
library(stringr)
library(dplyr)
library(future)
library(future.apply)
library(msigdbr)
library(clusterProfiler)
library(devtools)
library(harmony)
library(clustree)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(cowplot)
library(org.Hs.eg.db)
library(tidyverse)
library(dplyr)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(DESeq2)
library(AUCell)
library(scRNAtoolVis)

#define the color
library(ggsci)
cors <- pal_npg()(10) #定义颜色

#读取数据
setwd("E:/【科研学习】/【皮肤科】/白癜风课题/白癜风+黑色素瘤/思路3原始结果/QC_scRNA_data")
scRNA_T=readRDS('./scRNA_T&NK_new.rds')
table(scRNA_T$T_celltype,scRNA_T$group)


#### 提取T亚群 ####
scRNA_T=subset(scRNA,celltype=='T & NK')

# 提T细胞亚群,重新降维聚类
scRNA_T<- FindVariableFeatures(scRNA_T, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNA_T)
scRNA_T <- ScaleData(scRNA_T, features = scale.genes)
scRNA_T<- RunPCA(scRNA_T, features = VariableFeatures(scRNA_T))
DimPlot(scRNA_T, reduction = "pca", group.by = "orig.ident")
ElbowPlot(scRNA_T)

## 重新harmony
library(harmony)
set.seed(1000)
scRNA_T <- RunHarmony(scRNA_T, group.by.vars = "orig.ident")
DimPlot(scRNA_T, reduction = "harmony", group.by = "orig.ident")
ElbowPlot(scRNA_T,reduction = 'harmony')

scRNA_T <- FindNeighbors(scRNA_T, reduction = 'harmony',dims = 1:10)
scRNA_T <- FindClusters(scRNA_T, resolution = 1)
scRNA_T <- RunUMAP(scRNA_T, reduction = 'harmony',dims = 1:10)
scRNA_T <- RunTSNE(scRNA_T, reduction = 'harmony',dims = 1:10)
DimPlot(scRNA_T,reduction = "umap",label = T) 
dev.off()

# T细胞亚群注释参考1
# https://blog.csdn.net/weixin_52505487/article/details/126687526
Idents(scRNA_T)=scRNA_T$seurat_clusters
T_marker <- list('CD4+ Naive'=c("CCR7","LEF1", "TCF7",'SELL','KLF2'),
                 'CD4+ Teff'=c("KLRG1", "CX3CR1", "GZMH"),
                 'CD4+ Tcm'=c('PTGER2','ICAM2','ANXA2'),
                 'CD4+ Trm'=c('PTGER4','CXCR6','MYADM'),
                 'CD4+ Tfh'=c("BCL6", "CXCR5","ICA1"),
                 'CD4+ Th17'=c('IL23R',"CCR6",'RORC','IL17A'),
                 'CD4+ Treg'=c('FOXP3','IL2RA','IL1R2'),
                 'CD8+ Tem'=c('GZMK','CXCR3','CXCR4','CD44'),
                 'CD8+ Tcm'=c('IL7R','CD27','CD28','GZMA','CCL5','GPR183','S1PR1'),
                 'CD8+ Trm'=c('CD6','XCL1','XCL2','CAPG', 'RORA','CD69','ITGAE'),
                 'CD8+ Tex'=c('HAVCR2','PDCD1','LAG3'),
                 'NK'=c('GNLY','NKG7','PRF1'))
# T细胞亚群注释参考2
# https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-04221-8/MediaObjects/41586_2021_4221_MOESM1_ESM.pdf
Idents(scRNA_T)=scRNA_T$seurat_clusters
T_marker <- c("KLRB1","CD40LG", "CCR6",'RPL12','RORA', #CD4+ Teff
              'FOXP3','LAIR2','CTLA4','TIGIT','CUL9', #CD4+ Treg
              'CCL4','NKG7','GZMB','GZMH','GZMK',  #CD8+ Tcyt
              'ZNF683','TRGC2','CTSW','CCL5','FXYD2')  #CD8+ Tres
DotPlot(scRNA_T, 
        features = c('CD4'),
        group.by = "seurat_clusters") + coord_flip()
DotPlot(object = scRNA_T, 
        features = T_marker, 
        assay = "RNA") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))

# 粗分
T_celltype0=c('CD8+ T',
              'CD8+ T','CD8+ T','CD8+ T','CD8+ T','CD8+ T',
              'NK','CD4+ T','CD4+ T','CD4+ T','CD8+ T',
              'CD8+ T','NK','CD4+ T','CD8+ T')
# 第一种
T_celltype=c('CD8+ Tcm',
             'CD8+ Tem','CD8+ Trm','CD8+ Trm','CD8+ Trm','CD8+ Tcm',
             'NK','CD4+ Tcm','CD4+ Treg','CD4+ Treg','CD8+ Tem',
             'CD8+ Trm','NK','CD4+ Tfh','CD4+ Naive')
Idents(scRNA_T) <- scRNA_T@meta.data$seurat_clusters
names(T_celltype) <- levels(scRNA_T)
scRNA_T<- RenameIdents(scRNA_T, T_celltype)
scRNA_T@meta.data$T_celltype <- Idents(scRNA_T)
# 删除unknown
scRNA_T <- subset(scRNA_T, T_celltype != 'Unknown')
scRNA_T$T_celltype <- Idents(scRNA_T) #完全去除Unknown
table(scRNA_T$T_celltype)
#自定义celltype顺序
cell = c('CD4+ Naive','CD4+ Tcm',"CD4+ Tfh","CD4+ Treg",
         "CD8+ Tem","CD8+ Tcm",'CD8+ Trm','NK') 
#cell = c('CD4+ Teff','CD4+ Treg','CD8+ Tcyt','CD8+ Tres')
scRNA_T$T_celltype <- factor(scRNA_T$T_celltype,levels = cell)
Idents(scRNA_T) <- scRNA_T$T_celltype
table(scRNA_T$T_celltype)

saveRDS(scRNA_T, 'scRNA_T&NK_new.rds')

## tsne图
DimPlot(scRNA_T,reduction = "tsne",group.by = 'T_celltype',label = T,cols = cors) 

library(scRNAtoolVis)
pdf(file=paste0('T_celltype_tsne.pdf'),width = 7,height = 5)
clusterCornerAxes(object = scRNA_T, reduction = 'tsne',
                  clusterCol = "T_celltype", #分类依据
                  #noSplit = F, groupFacet = 'group', #分组时使用
                  cellLabel = T, cellLabelSize = 3.5)+  #标签
  #addCircle = TRUE, cicAlpha = 0.1, nbin = 200) + #画圈
  scale_color_npg() + scale_fill_npg()
dev.off()




##绘制细胞marker气泡图
pdf(file=paste0('T_celltype_dotplot.pdf'),width = 6,height = 8)
DotPlot(scRNA_T, features = T_marker, assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 60, hjust=1,vjust=1))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF')) #颜色
dev.off()
#美化气泡图  12*8
jjDotPlot(object = scRNA_T,gene = T_marker,
          id = 'T_celltype',  #分组设置
          ytree = F,
          dot.col = c('#3C5488FF','white','#E64B35FF'),
          rescale.min = 0,rescale.max = 2,midpoint = 0) #+ coord_flip()


####T细胞比例####
library(reshape2)
library(ggplot2)
#define the color
library(ggsci)
cors <- pal_npg()(10) #定义颜色

#准备绘图数据
prop_df <- table(scRNA_T$T_celltype,scRNA_T$group) %>% melt()
colnames(prop_df) <- c("Cluster","Sample","Number")
prop_df$Cluster <- factor(prop_df$Cluster)
prop_df$Proportion <- ave(prop_df$Number, prop_df$Sample, 
                          FUN = function(x) x/sum(x))
#比例图1  尺寸6*3
pdf(file=paste0('T_proportion_barplot.pdf'),width = 6,height = 3)
ggplot(data = prop_df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=cors) + #配色
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45))
dev.off()

#批量散点箱式图
#计算各组中各样本不同细胞群比例
Cellratio <- prop.table(table(Idents(scRNA_T), scRNA_T$orig.ident), margin = 2)
Cellratio <- data.frame(Cellratio)
colnames(Cellratio) <- c('T_celltype','Sample','Freq')
Cellratio$group <- ifelse(grepl('H',Cellratio$Sample),'Healthy',
                          ifelse(grepl('M',Cellratio$Sample),'Melanoma','Vitiligo'))

my_comparisons <- list(c("Healthy","Melanoma"),c("Healthy","Vitiligo"),c("Melanoma","Vitiligo"))

pdf(file=paste0('scRNA_T_proportion_boxplot.pdf'),width = 8,height = 4)
ggboxplot(Cellratio, x="group", y="Freq", width = 0.6, #按group分组
          color = "black",#轮廓颜色
          fill="group",#填充
          palette = cors,
          xlab = F, #不显示x轴的标签
          bxp.errorbar=T,#显示误差条
          bxp.errorbar.width=0.5, #误差条大小
          size=0.5, #箱型图边线的粗细
          outlier.shape=NA, #不显示outlier
          legend = "right") + 
  geom_jitter(width = 0.3, size = 1) + # 添加散点
  facet_wrap(~T_celltype, ncol=4, scales = "free_y") + # 按照细胞类型分面
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.border = element_rect(colour = "black", fill=NA)) #+
  #stat_compare_means(comparisons=my_comparisons,
                     #symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1),symbols=c("***", "**", "*", "ns")),label="p.signif") # 添加t检验
dev.off()

####CD8T组间DEG####
scRNA_CD4T <- scRNA_T[,!scRNA_T$T_celltype %in% c("CD8+ Tem","CD8+ Tcm","CD8+ Trm","NK")]
scRNA_CD8T <- subset(scRNA_T, T_celltype %in% c("CD8+ Tem","CD8+ Tcm","CD8+ Trm"))

## 1)pseudobulks差异分析####
bs = split(colnames(scRNA_CD8T),scRNA_CD8T$orig.ident)
ct = do.call(
  cbind,lapply(names(bs), function(x){ 
    kp =colnames(scRNA_CD8T) %in% bs[[x]]
    rowSums( as.matrix(scRNA_CD8T@assays$RNA@counts[, kp]  ))
  })
)
colnames(ct) <- names(bs)
phe = unique(scRNA_CD8T@meta.data[,c('orig.ident','group')])#样本&信息，自行修改
group_list = phe[match(names(bs),phe$orig.ident),'group']
exprSet = ct
exprSet=exprSet[apply(exprSet,1, function(x) sum(x>1) > 1),]
table(group_list)
#第一步构建DEseq对象
colData <- data.frame(row.names=colnames(exprSet),group_list=group_list)
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ group_list)
# 第二步，进行差异表达分析
dds2 <- DESeq(dds)
table(group_list)
tmp <- results(dds2,contrast=c("group_list","Healthy","Vitiligo")) #分组自行修改
DEG_DESeq2 <- as.data.frame(tmp[order(tmp$padj),])
DEG_DESeq2 = na.omit(DEG_DESeq2)
#添加上下调信息
DEG_DESeq2 <- DEG_DESeq2 %>%
  mutate(Type = if_else(pvalue > 0.05, "ns", #pvalue/padj
                        if_else(abs(log2FoldChange) < 0.25, "ns",
                                if_else(log2FoldChange >= 0.25, "up", "down")))) %>%
  arrange(desc(abs(log2FoldChange))) %>% rownames_to_column("Gene_Symbol")
table(DEG_DESeq2$Type)
write.csv(DEG_DESeq2,"CD8T_DEG_H_vs_V.csv")

## 2)FindMarkers组间DEG####
Idents(scRNA_CD8T)='group'
all_markers <- FindAllMarkers(object = scRNA_CD8T,
                              only.pos = T, #only.pos改为T则只输出高表达gene
                              min.pct = 0.25,logfc.threshold = 0.25)
top10_markers <- all_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top10_markers <- top10_markers[!duplicated(top10_markers$gene),]
# 各亚群平均表达量提取
genes <- unique(top10_markers$gene)
aver_dt <- AverageExpression(scRNA_CD8T,
                             features = genes,
                             group.by = 'group',
                             slot = 'data')
aver_dt <- as.data.frame(aver_dt$RNA)

#热图
library(ComplexHeatmap)
library(cols4all)
#热图配色自定义
mycol <- colorRampPalette(c("#4DBBD5FF", "white", "#E64B35FF"))(50)
#行列注释配色自定义
celltype_col <- c4a('10', 3) #此处改为细胞类型数目
gene_anno <- data.frame(gene_anno = top10_markers$cluster,#行注释：Top5marker对应celltype
                        row.names = top10_markers$gene)
cell_anno <- data.frame(cell_anno = colnames(aver_dt),#列注释：celltype
                        row.names = colnames(aver_dt))
names(celltype_col) <- cell_anno$cell_anno
anno_col <- list(cell_anno = celltype_col,gene_anno = celltype_col)

pdf(file="scRNA_CD8T_top10marker_heatmap.pdf", width=4, height=6)
pheatmap(as.matrix(aver_dt),
         scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col = cell_anno,
         annotation_colors = anno_col, #注释配色
         color = mycol, #热图配色
         #gaps_row = c(10,20), #描边颜色
         border_color = 'white') 
dev.off()

##更便捷方法
#法一：平均热图
library(scRNAtoolVis)
AverageHeatmap(object = scRNA_CD8T,
               markerGene = top10_markers$gene,
               clusterAnnoName = F, annoCol = TRUE, myanCol = cors[1:3],
               htCol = c("#3C5488FF", "white", "#DC0000FF"))
#法二：原始热图
pdf(file=paste0('scRNA_CD4Treg_DEG_heatmap.pdf'),width = 5,height = 4)
DoHeatmap(scRNA_CD8T, top10_markers$gene, 
          group.by = 'group', group.colors = cors) +
  scale_fill_gradientn(colors = c("#3C5488FF", "white", "#DC0000FF"))
dev.off()


####3)GSVA####
library(GSVA)
# 各亚群平均表达量提取
genes <- unique(rownames(scRNA_CD8T@assays$RNA@counts))
aver_dt <- AverageExpression(scRNA_CD8T,
                             features = genes,
                             group.by = 'group',
                             slot = 'data')
aver_dt <- as.data.frame(aver_dt$RNA)
aver_dt <- aver_dt[rowSums(aver_dt)>0,]  #过滤细胞表达量全为零的基因

#设置参考水平
group_list <- c('Healthy','Melanoma','Vitiligo')
group_list = factor(group_list,levels = c('Healthy','Melanoma','Vitiligo'))

##热图
library(pheatmap) 
annotation_col = data.frame(group=group_list)
rownames(annotation_col) <- colnames(exp_gsva)

pdf('CD8T_group_GSVA_interest.pdf', width = 5, height = 4.5)
pheatmap::pheatmap(exp_gsva,
                   show_colnames = T, 
                   cluster_rows = F, cluster_cols = F, 
                   color = colorRampPalette(c("#3C5488FF", "white", "#DC0000FF"))(50),
                   annotation_col = annotation_col)
dev.off()


####4)取趋势相反DEG交集####
#变化趋势相反的DEG
Idents(scRNA_CD8T)='group'
DEG_Viti <- FindMarkers(object = scRNA_CD8T,ident.1 = 'Vitiligo',ident.2 = 'Healthy',
                              only.pos = F, #only.pos改为T则只输出高表达gene
                              min.pct = 0.25,logfc.threshold = 0.25)
DEG_Mel <- FindMarkers(object = scRNA_CD8T,ident.1 = 'Melanoma',ident.2 = 'Healthy',
                               only.pos = F, #only.pos改为T则只输出高表达gene
                               min.pct = 0.25,logfc.threshold = 0.25)
DEG_Viti <- DEG_Viti %>%
  mutate(change = if_else(p_val_adj > 0.05, "ns", #pvalue/padj
                          if_else(abs(avg_log2FC) < 0.25, "ns",
                                  if_else(avg_log2FC >= 0.25, "up", "down"))))

DEG_Vitiligo_up = subset(DEG_Viti, change == 'up')
DEG_Vitiligo_down = subset(DEG_Viti, change == 'down')
DEG_Melanoma_up = subset(DEG_Mel, change == 'up')
DEG_Melanoma_down = subset(DEG_Mel, change == 'down')

inverse_DEG = union(intersect(rownames(DEG_Vitiligo_up), rownames(DEG_Melanoma_down)), 
                    intersect(rownames(DEG_Vitiligo_down), rownames(DEG_Melanoma_up)))
write.csv(inverse_DEG, 'CD8T_inverse_DEG.csv',row.names = F, col.names = F)

##韦恩图
library(ggvenn)
venn_dat <- list(
  "DEG_Vitiligo_up" = rownames(DEG_Vitiligo_up),
  "DEG_Melanoma_down" = rownames(DEG_Melanoma_down),
  "DEG_Vitiligo_down" = rownames(DEG_Vitiligo_down),
  "DEG_Melanoma_up" = rownames(DEG_Melanoma_up))

ggvenn(data = venn_dat,         # 数据列表
       columns = NULL,           # 对选中的列名绘图，最多选择4个，NULL为默认全选
       show_elements = F,        # 当为TRUE时，显示具体的交集情况，而不是交集个数
       label_sep = "\n",         # 当show_elements = T时生效，分隔符 \n 表示的是回车的意思
       show_percentage = F,      # 显示每一组的百分比
       digits = 1,               # 百分比的小数点位数
       fill_color = c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF"), # 填充颜色
       fill_alpha = 0.5,         # 填充透明度
       stroke_color = "white",   # 边缘颜色
       stroke_alpha = 0.5,       # 边缘透明度
       stroke_size = 0.5,        # 边缘粗细
       stroke_linetype = "solid", # 边缘线条 # 实线：solid  虚线：twodash longdash 点：dotdash dotted dashed  无：blank
       set_name_color = "black", # 组名颜色
       set_name_size = 4,        # 组名大小
       text_color = "black",     # 交集个数颜色
       text_size = 4)            # 交集个数文字大小

##差异基因可视化
#气泡图
pdf(file=paste0('CD8T_inverseDEG_dotplot.pdf'),width = 4.5,height = 5)
DotPlot(scRNA_CD8T, features = inverse_DEG, group.by = 'group',assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colors = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF')) #颜色
dev.off()
#关注ZNF683和ATP5I


####CD4Treg/CD8Trm组间DEG####
library(ggVolcano)
library(ggpubr)
library(ggthemes)
#提取某种细胞子集
scRNA_Trm <- subset(scRNA_T, T_celltype == 'CD8+ Trm')
#FindMarkers计算组间DEG
Trm_markers <- FindMarkers(scRNA_Trm,
                           ident.1 = 'Melanoma',ident.2 = 'Vitiligo', 
                           only.pos = F,min.pct = 0.1,logfc.threshold = 0) #0.25,0.25
Trm_markers <- Trm_markers %>%
  filter(p_val_adj > 0) %>%  #把p等于0的删除掉
  mutate(change = if_else(p_val_adj > 0.05, "ns", #pvalue/padj
                        if_else(abs(avg_log2FC) < 0.25, "ns",
                                if_else(avg_log2FC >= 0.25, "up", "down"))))
Trm_markers$logP <- -log10(Trm_markers$p_val_adj)
#添加火山图的基因标签
Trm_markers$Label = ""   #新加一列label
Trm_markers <- Trm_markers[order(Trm_markers$p_val_adj), ]  #对差异基因的p值进行从小到大的排序
Trm_markers$Gene <- rownames(Trm_markers)
up.genes <- head(Trm_markers$Gene[which(Trm_markers$change == "up")], 5) #高表达的基因中选择fdr值最小的5个
down.genes <- head(Trm_markers$Gene[which(Trm_markers$change == "down")], 5) #低表达的基因中选择fdr值最小的5个
#将up.genes和down.genes合并，加入到Label中
DEG.top5.genes <- c(as.character(up.genes), as.character(down.genes))
Trm_markers$Label[match(DEG.top5.genes, Trm_markers$Gene)] <- DEG.top5.genes
##火山图
pdf('CD8Trm_marker_volcano.pdf', width = 5, height = 4)
gradual_volcano(Trm_markers, x = "avg_log2FC", y = "p_val_adj",
                log2FC_cut = 0.5, FDR_cut = -log10(0.05),
                label = "Label", label_number = nrow(Trm_markers), output = FALSE)
dev.off()


####对角散点图####
#参考Onenote笔记


####亚群marker富集分析####
#只选择CD8 T分析
scRNA_CD8T = subset(scRNA_T,T_celltype %in% c('CD8+ Tem','CD8+ Tcm','CD8+ Trm'))
####0)提取亚群marker####
all_markers <- FindAllMarkers(object = scRNA_CD8T,
                              only.pos = T, #only.pos改为T则只输出高表达gene
                              min.pct = 0.25,logfc.threshold = 0.25) #画火山图时logfc.threshold = 0
top5_markers <- all_markers %>% group_by(cluster) %>% top_n(5, avg_log2FC) #输出每个组的前5个高的gene
top30_markers <- all_markers %>% group_by(cluster) %>% top_n(30, avg_log2FC)
# 各亚群平均表达量提取
genes <- unique(top5_markers$gene)
aver_dt <- AverageExpression(scRNA_CD8T,
                             features = genes,
                             group.by = 'T_celltype',
                             slot = 'data')
aver_dt <- as.data.frame(aver_dt$RNA)

##热图
library(ComplexHeatmap)
library(cols4all)
#热图配色自定义
mycol <- colorRampPalette(c("#4DBBD5FF", "white", "#E64B35FF"))(50)
#行列注释配色自定义
celltype_col <- c4a('10', 3) #此处改为细胞类型数目
gene_anno <- data.frame(gene_anno = top5_markers$cluster,#行注释：Top5marker对应celltype
                        row.names = top5_markers$gene)
cell_anno <- data.frame(cell_anno = colnames(aver_dt),#列注释：celltype
                        row.names = colnames(aver_dt))
names(celltype_col) <- cell_anno$cell_anno
anno_col <- list(cell_anno = celltype_col,gene_anno = celltype_col)

pdf(file="CD8T_markertop5_heatmap.pdf", height=4, width=4.5)
pheatmap(as.matrix(aver_dt),
         scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col = cell_anno,
         annotation_row = gene_anno,
         annotation_colors = anno_col, #注释配色
         color = mycol, #热图配色
         border_color = 'white') #描边颜色
dev.off()


##火山图
Trm_markers = all_markers[all_markers$cluster=='CD8+ Trm',]
Trm_markers <- Trm_markers %>%
  mutate(change = if_else(p_val_adj > 0.05, "ns", #pvalue/padj
                          if_else(abs(avg_log2FC) < 0.25, "ns",
                                  if_else(avg_log2FC >= 0.25, "up", "down"))))
Trm_markers$logP <- -log10(Trm_markers$p_val_adj)
#添加火山图的基因标签
Trm_markers$Label = ""   #新加一列label
Trm_markers <- Trm_markers[order(Trm_markers$p_val_adj), ]  #对差异基因的p值进行从小到大的排序
Trm_markers$Gene <- rownames(Trm_markers)
up.genes <- head(Trm_markers$Gene[which(Trm_markers$change == "up")], 5) #高表达的基因中选择fdr值最小的5个
down.genes <- head(Trm_markers$Gene[which(Trm_markers$change == "down")], 5) #低表达的基因中选择fdr值最小的5个
#将up.genes和down.genes合并，加入到Label中
DEG.top5.genes <- c(as.character(up.genes), as.character(down.genes))
Trm_markers$Label[match(DEG.top5.genes, Trm_markers$Gene)] <- DEG.top5.genes
#绘图
pdf('CD8Trm_marker_volcano.pdf', width = 5, height = 4)
gradual_volcano(Trm_markers, x = "avg_log2FC", y = "p_val_adj",
                log2FC_cut = 0.5,
                label = "Label", label_number = nrow(Trm_markers), output = FALSE)
dev.off()


####1)GO分析####
Tem_markers = all_markers[all_markers$cluster=='CD8+ Tem',]$gene
Tcm_markers = all_markers[all_markers$cluster=='CD8+ Tcm',]$gene
Trm_markers = all_markers[all_markers$cluster=='CD8+ Trm',]$gene
#产生ENTREZID
genelist <- bitr(Tem_markers, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')

ego1 <- enrichGO(gene = genelist$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "BP",            #此处选择BP
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, qvalueCutoff =0.2,
                readable = TRUE)
ego1_res <- ego1@result
ego1_res$Group <- 'CD8+ Tem'
ego1_res$LogP <- -log(ego1_res$p.adjust) #计算-logP
##ego2_res,ego3_res同理计算

#只提取与xxx有关的基因集   JAK|IFN|cytokine|chemokine|oxidative|pathway|T cell|immune
ego1_res = ego1_res  %>%
  dplyr::filter(grepl('pathway|signal|T cell|immune|oxidative|antigen|mitochondrion', Description)) %>%
  dplyr::arrange(dplyr::desc(LogP),dplyr::desc(Description)) %>%
  mutate(Description =forcats:: fct_inorder(Description))

## 气泡图
barplot(ego1, showCategory = 10)
# 个性化柱状图
pdf('CD8Tem_GO_BP_barplot.pdf', width = 6, height = 5)
ggplot(ego1_res[1:10,],aes(Count,Description))+ #只展示前十条通路
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=LogP),stat='identity')+
  scale_fill_gradient2(low="#486b98", mid="#f5f2b1", high="#b93735", midpoint=8)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=0.5))
dev.off()

# 合并富集结果
ego_res = rbind(ego1_res[1:10,],ego2_res[1:10,],ego3_res[1:10,])

# 分组柱状图
pdf('scRNA_CD8T_GOBP_barplot.pdf', width = 6, height = 6)
ggplot(ego_res, aes(Count, Description)) +
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=LogP),stat='identity') +
  scale_fill_gradient(low="#F39B7FFF",high="#DC0000FF")+
  facet_wrap(~Group, scales="free", ncol=1)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=1))
dev.off()

# 分组气泡图
pdf('scRNA_CD8T_GOBP_dotplot.pdf', width = 6.5, height = 4)
ggplot(ego_res, aes(Group, Description)) +
  geom_point(aes(color=LogP, size=Count))+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=1))+
  scale_color_gradient(low='#4DBBD5FF',high='#DC0000FF')+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1))
dev.off()


####2)KEGG分析####
kk1 <- enrichKEGG(gene = genelist$ENTREZID,
                 keyType = 'kegg',
                 organism = 'hsa',
                 pvalueCutoff = 0.1,
                 qvalueCutoff =0.1)
kk1_res <- kk1@result
kk1_res$Group <- 'CD8+ Tem'
kk1_res$LogP <- -log(kk1_res$p.adjust)
##kk2_res,kk3_res同理计算
##绘图同理



####3)GSVA分析####
library(GSVA)
# 各亚群平均表达量提取
genes <- unique(rownames(scRNA_CD8T@assays$RNA@counts))
aver_dt <- AverageExpression(scRNA_CD8T,
                             features = genes,
                             group.by = 'celltype.group',
                             slot = 'data')
aver_dt <- as.data.frame(aver_dt$RNA)
aver_dt <- aver_dt[rowSums(aver_dt)>0,]  #过滤细胞表达量全为零的基因

##法一：选择msigdb的KEGG基因集
#msigdbr包教程：https://www.jianshu.com/p/f2febb3123d8
geneSet = msigdbr(species = "Homo sapiens")  #category = "C2"
geneSet %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
#挑选感兴趣的基因集  主流通路/NO/ROS/衰老/自噬/线粒体
geneSet_name = c("KEGG_JAK_STAT_SIGNALING_PATHWAY",
                 'KEGG_MAPK_SIGNALING_PATHWAY',
                 "KEGG_HEDGEHOG_SIGNALING_PATHWAY",
                 'KEGG_NOTCH_SIGNALING_PATHWAY',
                 'KEGG_WNT_SIGNALING_PATHWAY',
                 'KEGG_TGF_BETA_SIGNALING_PATHWAY',
                 'KEGG_OXIDATIVE_PHOSPHORYLATION',
                 'KEGG_P53_SIGNALING_PATHWAY',
                 'KEGG_MTOR_SIGNALING_PATHWAY',
                 'KEGG_CHEMOKINE_SIGNALING_PATHWAY',
                 #'KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION',
                 'REACTOME_STING_MEDIATED_INDUCTION_OF_HOST_IMMUNE_RESPONSES',
                 'REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM',
                 'REACTOME_CELLULAR_SENESCENCE',
                 'BIOCARTA_NOS1_PATHWAY',
                 'BIOCARTA_IFNG_PATHWAY', #REACTOME_INTERFERON_GAMMA_SIGNALING
                 'GOBP_RESPONSE_TO_OXIDATIVE_STRESS',
                 #'GOBP_T_CELL_MEDIATED_CYTOTOXICITY',
                 'GOBP_AUTOPHAGY_OF_MITOCHONDRION',
                 'GOBP_RESPONSE_TO_TYPE_I_INTERFERON',
                 'WP_FERROPTOSIS',
                 'BIOCARTA_CYTOKINE_PATHWAY',
                 'GOBP_IMMUNE_RESPONSE')
#趋化/迁移相关
geneSet_name = c('GOBP_REGULATION_OF_T_CELL_CHEMOTAXIS',
                 'GOBP_T_CELL_EXTRAVASATION',  #T cell migration
                 'GOBP_REGULATION_OF_MACROPHAGE_CHEMOTAXIS',
                 'GOBP_MONOCYTE_EXTRAVASATION')

my_geneSet = subset(geneSet, gs_name %in% geneSet_name)
my_geneSet = my_geneSet %>% split(x = .$gene_symbol, f = .$gs_name)#基因集是list

#只提取与xxx有关的基因集（不运行）
geneSet = geneSet  %>%
  dplyr::filter(stringr::str_detect(pattern = "IL", gene_symbol))
geneSet = geneSet %>% split(x = .$gene_symbol, f = .$gs_name)#基因集是list
geneSet = split(geneSet$gene, geneSet$term)

##法二：读取免疫/代谢相关通路基因集
#教程1：https://blog.csdn.net/arrhythmia08/article/details/130772130
#教程2：https://www.jianshu.com/p/461c0dc1297a
#KEGG官网：https://www.kegg.jp/kegg/pathway.html
##提取KEGG所有的通路
library(KEGGREST)
hsa_path <- keggLink("pathway","hsa")
length(unique(names(hsa_path))) #KEGG共收集8443个相关通路的基因
length(unique(hsa_path)) #涉及通路353个
#代谢相关通路在KEGG中以00开头，免疫相关通路以046开头
pathway=unique(hsa_path)[grepl('hsa046',unique(hsa_path))]
length(pathway)
#提取这些通路的基因
hsa_info <- lapply(pathway, keggGet)
#提取基因名字
gene_symbol = unlist(lapply( hsa_info , function(x) {
  g = x[[1]]$GENE
  str_split(g[seq(2,length(g),by=2)],';',simplify = T)[,1]
}))
#提取基因通路
gene_pathway = unlist(lapply( hsa_info , function(x) {
  g = x[[1]]$PATHWAY_MAP
}))
gene_num = NULL
for (i in 1:length(gene_pathway)){
  gene_num = c(gene_num,length(hsa_info[[i]][[1]]$GENE)/2)
}
gene_pathway = rep(gene_pathway, times = gene_num)
#合成一个数据框
genelist <- data.frame(gene_pathway,gene_symbol)
genelist <- genelist[!duplicated(genelist$gene_symbol),] #去重复
##读取基因集
colnames(genelist) <- c('term','gene')
geneSet = split(genelist$gene, genelist$term)
str(head(geneSet))


##GSVA分析 
exp_gsva <- gsva(as.matrix(aver_dt),  #表达矩阵需转换为matrix格式
                 gset.idx.list = my_geneSet) #method=c("gsva","ssgsea","zscore","plage")

#选择相应表达水平中位绝对偏差中位数排名前20位的通路（看情况运行）
mad_scores <- apply(exp_gsva, 1, mad)
top_genes <- order(mad_scores, decreasing = TRUE)[1:20]
exp_gsva_top <- exp_gsva[top_genes, ]

#调整exp_gsva行顺序
rownames(exp_gsva)
exp_gsva <- exp_gsva[c(6,2,1,3:5,7:20),]
rownames(exp_gsva) <- c('Type I IFN pathway',
                        'IFN-γ pathway',
                        'Cytokine network',
                        'NOS1 pathway',
                        'Autophagy of mitochondrion',
                        'Response to oxidative stress',
                        'Chemokine signaling',
                        "HEDGEHOG signaling",
                        "JAK STAT signaling",
                        'MAPK signaling',
                        'MTOR signaling',
                        'NOTCH signaling',
                        'Oxidative Phosphorylation',
                        'P53  signaling',
                        'TGF BETA signaling',
                        'WNT signaling',
                        'Cellular senescence ',
                        'Immune system cytokine signaling',
                        'STING mediated immune response',
                        'Ferroptosis')


#设置参考水平
group_list <- c('CD8+ Tem','CD8+ Tcm','CD8+ Trm')
group_list = factor(group_list,levels = c('CD8+ Tem','CD8+ Tcm','CD8+ Trm'))

##热图
library(pheatmap) 
annotation_col = data.frame(group=group_list)
rownames(annotation_col) <- colnames(exp_gsva)

pdf('CD8T_celltype_group_GSVA_interest.pdf', width = 7, height = 4.5)
pheatmap::pheatmap(exp_gsva,
                   show_colnames = F, 
                   cluster_rows = F, cluster_cols = F, 
                   #gaps_row = c(3,5,14),  #添加行分隔
                   gaps_col = c(3,6,9),  #添加列分隔
                   color = colorRampPalette(c("#3C5488FF", "white", "#DC0000FF"))(50),
                   annotation_col = annotation_col)
dev.off()



####AUCell关键通路评分####
###AUCell与GSVA方法结果可能有不一致！
##通过R获取msigdbr数据集
#geneSet和GSVA一样  I/II型干扰素/细胞毒性/迁移

cells_AUC <- AUCell_run(scRNA_CD8T@assays$RNA@data, my_geneSet)
#提取PATHWAY
geneSet <- 'GOBP_RESPONSE_TO_TYPE_I_INTERFERON'
AUCell_auc <- as.numeric(getAUC(cells_AUC)[geneSet, ])
scRNA_CD8T$AUCell <- AUCell_auc #添加至metadata中
head(scRNA_CD8T@meta.data)

##绘图可视化
#Seurat自带小提琴图
pdf(file='CD8T_CYTOTOXICITY_vlnplot.pdf',width = 5,height = 4)
VlnPlot(scRNA_CD8T,features = 'AUCell', #features也可改为AUCell
        pt.size = 0,group.by = "T_celltype",col = cors) #按细胞类型分组
dev.off()

#箱式图
my_comparisons = list(c("Healthy","Melanoma"),c("Healthy","Vitiligo"),c("Melanoma","Vitiligo"))

pdf(file='CD8T_typeI_IFN_boxplot.pdf',width = 7,height = 3)
ggboxplot(scRNA_CD8T@meta.data, x="group", y="AUCell", width = 0.6, #按group分组
          color = "black",#轮廓颜色
          fill="group",#填充
          palette = cors,
          xlab = F, ylab = "typeI_IFN",  #标题
          bxp.errorbar=T,#显示误差条
          bxp.errorbar.width=0.5, #误差条大小
          size=0.5, #箱型图边线的粗细
          outlier.shape=NA, #不显示outlier
          legend = "right") + 
  #ylim(0,0.1) +
  facet_wrap(~T_celltype,ncol=4,scales="free_y") + # 按照细胞类型分面
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  stat_compare_means(comparisons=my_comparisons,
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1),symbols=c("***", "**", "*", "ns")),label="p.signif") # 添加t检验
dev.off()

## 小提琴图
ggviolin(scRNA@meta.data, x="group", y="GENE", width = 0.6, #按group分组
         color = "black",#轮廓颜色
         fill="group",#填充
         palette = cors,
         xlab = F, ylab = "Normalized Expression", #不显示x轴的标签
         bxp.errorbar=T,#显示误差条
         bxp.errorbar.width=0.5, #误差条大小
         size=0.5, #箱型图边线的粗细
         outlier.shape=NA, #不显示outlier
         legend = "right",
         add = "boxplot", 
         add.params = list(fill = "white", #填充颜色
                           width = 0.1,#宽度
                           linetype = 1)) + 
  facet_wrap(~celltype,ncol=4,scales="free_y") + # 按照细胞类型分面
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  stat_compare_means(comparisons=my_comparisons,
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1),symbols=c("***", "**", "*", "ns")),label="p.signif") # 添加t检验


## tsne图
# 法二  https://zhuanlan.zhihu.com/p/482523999
library(ggrepel)
#提取tsne坐标数据
tsne <- data.frame(scRNA_CD8T@meta.data, scRNA_CD8T@reductions$tsne@cell.embeddings)
library(ggplot2)
scRNA_CD8T$tsne_1 <- tsne$tSNE_1
scRNA_CD8T$tsne_2 <- tsne$tSNE_2
mydata <- FetchData(scRNA_CD8T,vars = c("tsne_1","tsne_2","AUCell"))
#绘图
pdf(file='CD8T_CYTOTOXICITY_tsne.pdf',width = 4,height = 3)
ggplot(mydata,aes(x = tsne_1,y =tsne_2,colour = AUCell))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('#333366',"#6666FF",'#CC3333','#FFCC33')) + 
  theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                       panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"))
dev.off()


###感兴趣通路基因dotplot
#趋化因子，炎症因子
genes <- unique(geneSet$GOBP_T_CELL_EXTRAVASATION)

pdf(file=paste0('CD8T_migration_dotplot.pdf'),width = 5,height = 4)
DotPlot(scRNA_CD8T, features = genes, group.by = 'T_celltype', assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF')) #颜色
dev.off()


####Treg细分亚群####
scRNA_Treg=subset(scRNA_T,T_celltype=='CD4+ Treg')
# 提Treg细胞亚群,重新降维聚类
scRNA_Treg<- FindVariableFeatures(scRNA_Treg, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNA_Treg)
scRNA_Treg <- ScaleData(scRNA_Treg, features = scale.genes)
scRNA_Treg<- RunPCA(scRNA_Treg, features = VariableFeatures(scRNA_Treg))
DimPlot(scRNA_Treg, reduction = "pca", group.by = "orig.ident")
ElbowPlot(scRNA_Treg)
## 重新harmony
library(harmony)
set.seed(1000)
scRNA_Treg <- RunHarmony(scRNA_Treg, group.by.vars = "orig.ident")
DimPlot(scRNA_Treg, reduction = "harmony", group.by = "orig.ident")
ElbowPlot(scRNA_Treg,reduction = 'harmony')

scRNA_Treg <- FindNeighbors(scRNA_Treg, reduction = 'harmony',dims = 1:10)
scRNA_Treg <- FindClusters(scRNA_Treg, resolution = 1)
scRNA_Treg <- RunUMAP(scRNA_Treg, reduction = 'harmony',dims = 1:10)
scRNA_Treg <- RunTSNE(scRNA_Treg, reduction = 'harmony',dims = 1:10)
DimPlot(scRNA_Treg,reduction = "tsne",label = T) 
dev.off()

# Treg细胞亚群注释参考
# https://blog.csdn.net/weixin_52505487/article/details/126687526
Idents(scRNA_Treg)=scRNA_Treg$seurat_clusters
T_marker <- list('Th1-like Treg'=c("IFNG","CXCR3", "CCR5"),
                 'Th2-like Treg'=c("IL13", "IL4", "IRF4",'GATA3'),
                 'Th17-like Treg'=c('IL17A','CCR6','RORC'),
                 'Tfh-like Treg'=c('BCL6','CXCR5'))
#https://mp.weixin.qq.com/s/M_ZhR_u1HBihAQYJl7F7Ww
T_marker <- list('Treg-FOXP3'=c('IL2RA', 'IKZF2'),
                 'Treg-LAG3'=c('LAG3', 'HAVCR2', 'PDCD1', 'TIGIT'))
#https://www.cyagen.com/cn/zh-cn/community/frontier/information-20171121-1.html
T_marker <- c('IL35','CCR7',  #IL-35-Treg
              'IL10','CCR5','CCR4')   #IL-10-Treg

#区分nTreg还是iTreg  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3708155/  
DotPlot(scRNA_Treg, 
        features = 'FOXP3', 
        group.by = "seurat_clusters") + coord_flip()
DotPlot(object = scRNA_Treg, 
        features = T_marker, 
        assay = "RNA") + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

# 第一种
Treg_celltype=c('iTreg',
                'iTreg','nTreg','nTreg','nTreg','nTreg',
                'nTreg','nTreg','nTreg')
Idents(scRNA_Treg) <- scRNA_Treg@meta.data$seurat_clusters
names(Treg_celltype) <- levels(scRNA_Treg)
scRNA_Treg<- RenameIdents(scRNA_Treg, Treg_celltype)
scRNA_Treg@meta.data$Treg_celltype <- Idents(scRNA_Treg)
table(scRNA_Treg$Treg_celltype)

DotPlot(scRNA_Treg, 
        features = c('IL10','IL35','TGFB1'), 
        group.by = "group") + coord_flip()



####NK细分亚群####
scRNA_NK=subset(scRNA_T,T_celltype=='NK')
# 提Treg细胞亚群,重新降维聚类
scRNA_NK<- FindVariableFeatures(scRNA_NK, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNA_NK)
scRNA_NK <- ScaleData(scRNA_NK, features = scale.genes)
scRNA_NK<- RunPCA(scRNA_NK, features = VariableFeatures(scRNA_NK))
DimPlot(scRNA_NK, reduction = "pca", group.by = "orig.ident")
ElbowPlot(scRNA_NK)
## 重新harmony
library(harmony)
set.seed(1000)
scRNA_NK <- RunHarmony(scRNA_NK, group.by.vars = "orig.ident")
DimPlot(scRNA_NK, reduction = "harmony", group.by = "orig.ident")
ElbowPlot(scRNA_NK,reduction = 'harmony')

scRNA_NK <- FindNeighbors(scRNA_NK, reduction = 'harmony',dims = 1:10)
scRNA_NK <- FindClusters(scRNA_NK, resolution = 0.2)
scRNA_NK <- RunUMAP(scRNA_NK, reduction = 'harmony',dims = 1:10)
scRNA_NK <- RunTSNE(scRNA_NK, reduction = 'harmony',dims = 1:10)
DimPlot(scRNA_NK,reduction = "tsne",split.by = 'group',label = T) 
dev.off()

# 细胞亚群注释
#准备基因集
geneSet = msigdbr(species = "Homo sapiens")  #category = "C2"
geneSet %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
geneSet_name = 'WP_OXIDATIVE_STRESS_RESPONSE'
my_geneSet = subset(geneSet, gs_name %in% geneSet_name)
my_geneSet = my_geneSet %>% split(x = .$gene_symbol, f = .$gs_name)
#AUCell评分
cells_AUC <- AUCell_run(scRNA_NK@assays$RNA@data, my_geneSet)
geneSet <- "WP_OXIDATIVE_STRESS_RESPONSE"   
AUCell_auc <- as.numeric(getAUC(cells_AUC)[geneSet, ])
scRNA_NK$AUCell <- AUCell_auc 
#绘图可视化
VlnPlot(scRNA_NK,features = 'AUCell', #features也可改为AUCell
        pt.size = 0,group.by = "seurat_clusters",col = cors) 
ggboxplot(scRNA_NK@meta.data, x="group", y="AUCell", width = 0.6, #按group分组
          color = "black",#轮廓颜色
          fill="group",#填充
          palette = cors,
          xlab = F, #不显示x轴的标签
          bxp.errorbar=T,#显示误差条
          bxp.errorbar.width=0.5, #误差条大小
          size=0.5, #箱型图边线的粗细
          outlier.shape=NA, #不显示outlier
          legend = "right") + 
  #ylim(0,0.2) +
  #facet_wrap(~T_celltype,ncol=4,scales="free_y") + # 按照细胞类型分面
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  stat_compare_means(comparisons=my_comparisons,
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1),symbols=c("***", "**", "*", "ns")),label="p.signif") # 添加t检验



####T细胞活化相关基因比较####
scRNA_T$celltype.group <- paste(scRNA_T$T_celltype, scRNA_T$group, sep = "_")

geneset_activate = c('IFNG','TNF','IL2', #Th1
                      'IL4','IL5','IL6','IL10','IL13', #Th2
                      'IL17A','IL21','IL22','IL26')  #Th17

geneset_cytokine = c('TNF', 'IL1B', 'IL6', 'IL12', 'IL18', #细胞因子
                     'IL10', 'TGFB1',
                     'IL2RA', 'IL4R', 'IL6R', 'IL10RA', 'IL12RB1', 'TGFBRII',
                     'CCL2','CXCL8', 'CCR5', 'CXCR1')

geneset_check = c('CTLA4', 'PD1', 'PDL1', 'LAG3', 'TIM3', 'TIGIT', #检查点免疫通路
                  'PIK3CD', 'PTEN', 'SHP2')

geneset_exhauste = c('PDCD1', 'LAG3', 'HAVCR2', 'CD244','CD160')

geneset_ILR = c('IL2RA','IL2RB','IL2RG','NGFR','GYPE',
                #'IL17RA','IL17RB','IL17RC','IL17RD','IL17RE',
                'CCR3','CXCR3','SDC4','ACKR3')

geneset_CXCLR = c('CCR3','CXCR3','SDC4')

scRNA_CD8Tcm = subset(scRNA_T, T_celltype == 'CD8+ Tcm')
scRNA_CD8Tem = subset(scRNA_T, T_celltype == 'CD8+ Tem')
scRNA_CD8Trm = subset(scRNA_T, T_celltype == 'CD8+ Trm')

p1=DotPlot(scRNA_CD8T, features = geneset_exhauste, group.by = 'group', assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust =1,vjust=1))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF'))
p2=DotPlot(scRNA_CD8Tem, features = geneset_check, group.by = 'group', assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust =1,vjust=1))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF'))
p3=DotPlot(scRNA_CD8Trm, features = geneset_check, group.by = 'group', assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust =1,vjust=1))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF'))
p1+p2+p3
#美化气泡图  9*7
jjDotPlot(object = scRNA_CD8T,gene = geneset_ILR,
          id = 'group',  #分组设置
          ytree = F,
          dot.col = c('#3C5488FF','white','#E64B35FF'),
          rescale.min = 0,rescale.max = 2,midpoint = 0)  + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust =1,vjust=1))+ coord_flip()

DotPlot(scRNA_CD8T, features = geneset_ILR, group.by = 'group', assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust =1,vjust=1))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF'))


#### 桑基图 ####
## 获取T细胞的相互作用对
library(ggplot2)
library(ggalluvial)
library(celltalker)

molecules <- ramilowski_pairs
molecules <- molecules %>% filter(ligand %in% c('CXCL9','CXCL10','IL2'))

# 创建配体和受体组合的标识
#molecules$pair <- with(molecules, paste(ligand, receptor, sep = "_"))

# 使用ggplot2绘制桑基图
ggplot(data = molecules,
       aes(axis1 = ligand, axis2 = receptor)) + #weight可省略
  geom_alluvium(aes(fill = pair)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_void() + scale_color_igv() + scale_fill_igv()
ggsave("IL-ILR.pdf", width = 5, height = 4)
