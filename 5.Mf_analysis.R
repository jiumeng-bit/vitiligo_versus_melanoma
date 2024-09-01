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

#define the color
library(ggsci)
cors <- pal_npg()(10) #定义颜色

#读取数据
setwd("E:/【科研学习】/【皮肤科】/白癜风课题/白癜风+黑色素瘤/思路3原始结果/QC_scRNA_data")
scRNA_Myeloid=readRDS('./scRNA_Myeloid.rds')
table(scRNA_Myeloid$Mono_celltype,scRNA_Myeloid$group)

#### 提取单核吞噬细胞 ####
#https://zhuanlan.zhihu.com/p/375318689 朗格汉斯细胞只是DC的一个亚群！
scRNA_Mf=subset(scRNA,celltype=='Mono phagocyte') 
# 提细胞亚群,重新降维聚类
scRNA_Mf <- FindVariableFeatures(scRNA_Mf, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNA_Mf)
scRNA_Mf <- ScaleData(scRNA_Mf, features = scale.genes)
scRNA_Mf<- RunPCA(scRNA_Mf, features = VariableFeatures(scRNA_Mf))
DimPlot(scRNA_Mf, reduction = "pca", group.by = "orig.ident")
ElbowPlot(scRNA_Mf)

## 重新harmony
library(harmony)
set.seed(1000)
scRNA_Mf <- RunHarmony(scRNA_Mf, group.by.vars = "orig.ident")
DimPlot(scRNA_Mf, reduction = "harmony", group.by = "orig.ident")
ElbowPlot(scRNA_Mf,reduction = 'harmony')

scRNA_Mf <- FindNeighbors(scRNA_Mf, reduction = 'harmony',dims = 1:10)
sce_res <- scRNA_Mf
for (i in c(0.01, 0.05, 0.1, 0.15, 0.2, 0.3,0.4, 0.5,0.8,1)){
  sce_res <- FindClusters(sce_res,resolution = i)
}
scRNA_Mf <- FindClusters(scRNA_Mf, resolution = 0.8)
scRNA_Mf <- RunUMAP(scRNA_Mf, reduction = 'harmony',dims = 1:10)
scRNA_Mf <- RunTSNE(scRNA_Mf, reduction = 'harmony',dims = 1:10)
DimPlot(scRNA_Mf,reduction = "tsne",group.by = "seurat_clusters",label = T) 

## 单核吞噬细胞亚群注释参考
# marker：https://www.jianshu.com/p/16280b2da028,https://zhuanlan.zhihu.com/p/520281423
# Mf亚群marker：https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9882988/
# 背景知识：https://zhuanlan.zhihu.com/p/653640358
Idents(scRNA_Mf)=scRNA_Mf$seurat_clusters
# 最全注释
Mono_marker <- c('CXCL5','CCL23','CCL3','CCL20','CCL4', # Mφ_CCL3
                 'SPP1','MMP9','FPR3', # Mφ_SPP1
                 "C1QA","C1QB",'TREM2','APOE','APOC1','GPNMB',  # Mφ_APOE
                 'IL1B','IL1RN','FCN1','DUSP6', # Mono_FCN1
                 'FCGR3A','CX3CR1','TCF7L2','LRRC25','FCGR3B') # Mono_CX3CR1
# 简易版注释
Mono_marker <- c('CXCL8','APOE','FCN1')
DotPlot(scRNA_Mf, 
        features = Mono_marker,
        group.by = "seurat_clusters") + coord_flip()

##如果文献中提供的marker效果不好，则直接按各cluster的top marker编号
all_markers <- FindAllMarkers(object = scRNA_Mf,
                              only.pos = T, #only.pos改为T则只输出高表达gene
                              min.pct = 0.25,logfc.threshold = 0.5)
# 提取各亚群Top5Marker
top10_markers <- all_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC) #输出每个组的前5个高的gene

Mono_celltype=c('Mφ_CCL3',
                'Mono_FCN1','Mφ_CCL3','Mφ_MMP9','Mφ_APOE','Mφ_APOE',
                'Mφ_CCL3','Mφ_MMP9','Mφ_MMP9','Mφ_CCL3','Unknown',
                'Unknown','Unknown','Unknown','Mono_FCN1','Unknown')

Idents(scRNA_Mf) <- scRNA_Mf@meta.data$seurat_clusters
names(Mono_celltype) <- levels(scRNA_Mf)
scRNA_Mf <- RenameIdents(scRNA_Mf, Mono_celltype)
scRNA_Mf$Mono_celltype <- Idents(scRNA_Mf)
# 删除unknown
scRNA_Mf <- subset(scRNA_Mf, Mono_celltype != 'Unknown')
scRNA_Mf$Mono_celltype <- Idents(scRNA_Mf) #完全去除Unknown
table(scRNA_Mf$Mono_celltype)
#自定义celltype顺序
cell = c('Mφ_CCL3','Mφ_APOE','Mφ_MMP9','Mono_FCN1') 
scRNA_Mf$Mono_celltype <- factor(scRNA_Mf$Mono_celltype,levels = cell)
Idents(scRNA_Mf) <- scRNA_Mf$Mono_celltype
saveRDS(scRNA_Mf, 'scRNA_Mf_new.rds')


#### 提取Langerhans细胞 ####
scRNA_DC=subset(scRNA,celltype=='Langerhans')
# 提细胞亚群,重新降维聚类(各做)
scRNA_DC <- FindVariableFeatures(scRNA_DC, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNA_DC)
scRNA_DC <- ScaleData(scRNA_DC, features = scale.genes)
scRNA_DC<- RunPCA(scRNA_DC, features = VariableFeatures(scRNA_DC))
DimPlot(scRNA_DC, reduction = "pca", group.by = "orig.ident")
ElbowPlot(scRNA_DC)

## 重新harmony
library(harmony)
set.seed(1000)
scRNA_DC <- RunHarmony(scRNA_DC, group.by.vars = "orig.ident")
DimPlot(scRNA_DC, reduction = "harmony", group.by = "orig.ident")
ElbowPlot(scRNA_DC,reduction = 'harmony')

scRNA_DC <- FindNeighbors(scRNA_DC, reduction = 'harmony',dims = 1:10)
scRNA_DC <- FindClusters(scRNA_DC, resolution = 0.8)
scRNA_DC <- RunUMAP(scRNA_DC, reduction = 'harmony',dims = 1:10)
scRNA_DC <- RunTSNE(scRNA_DC, reduction = 'harmony',dims = 1:10)
DimPlot(scRNA_DC,reduction = "tsne",group.by = "seurat_clusters",label = T) 
dev.off()

Idents(scRNA_DC)=scRNA_DC$seurat_clusters
# Langerhans注释 
 #https://pubmed.ncbi.nlm.nih.gov/34508661/
DC_marker <- c('CD207','CD1A', #LC1
               'CD1C', 'CD1B','CLEC10A', #LC2
               'CD83','MMP9','PIM3', #activated LCs (aLC)
               'CCR7')#migratory LCs (migLC)
# 传统DC注释
DC_marker <- c("GZMB","TSPAN13","LILRA4","ITM2C", # pDC
               "CLEC9A","CPNE3", "HLA-DRB1", # cDC1
               'FCER1A',"CD1C",'CD1A','CD1E','HLA-DQA1', # cDC2
               "LAMP3","FSCN1","CD83","CSF2RA","CCR7") # cDC3

DotPlot(scRNA_DC, 
        features = DC_marker,
        group.by = "seurat_clusters") + coord_flip()

Mono_celltype=c('aLC',
                'LC1','LC1','aLC','LC2','LC1',
                'LC1','LC1','LC1','aLC','Unknown',
                'LC1')

Idents(scRNA_DC) <- scRNA_DC$seurat_clusters
names(Mono_celltype) <- levels(scRNA_DC)
scRNA_DC <- RenameIdents(scRNA_DC, Mono_celltype)
scRNA_DC$Mono_celltype <- Idents(scRNA_DC)
# 删除unknown
scRNA_DC <- subset(scRNA_DC, Mono_celltype != 'Unknown')
scRNA_DC$Mono_celltype <- Idents(scRNA_DC) #完全去除Unknown
table(scRNA_DC$Mono_celltype)
#自定义celltype顺序
cell = c('LC1',"LC2","aLC") 
scRNA_DC$Mono_celltype <- factor(scRNA_DC$Mono_celltype,levels = cell)
Idents(scRNA_DC) <- scRNA_DC$Mono_celltype
DimPlot(scRNA_DC,reduction = "tsne",label = T)
saveRDS(scRNA_DC, 'scRNA_DC_new.rds')


#### 合并Myeloid亚群 ####
scRNA_Myeloid = merge(scRNA_Mf,scRNA_DC)
table(scRNA_Myeloid$Mono_celltype)
#自定义celltype顺序
cell = c('Mφ_CXCL8','Mφ_APOE','Mono_FCN1','LC1',"LC2","aLC") 
scRNA_Myeloid$Mono_celltype <- factor(scRNA_Myeloid$Mono_celltype,levels = cell)
Idents(scRNA_Myeloid) <- scRNA_Myeloid$Mono_celltype

# 提细胞亚群,重新降维聚类
scRNA_Myeloid<- FindVariableFeatures(scRNA_Myeloid, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNA_Myeloid)
scRNA_Myeloid <- ScaleData(scRNA_Myeloid, features = scale.genes)
scRNA_Myeloid<- RunPCA(scRNA_Myeloid, features = VariableFeatures(scRNA_Myeloid))
DimPlot(scRNA_Myeloid, reduction = "pca", group.by = "orig.ident")
ElbowPlot(scRNA_Myeloid)

## 重新harmony
library(harmony)
set.seed(1000)
scRNA_Myeloid <- RunHarmony(scRNA_Myeloid, group.by.vars = "orig.ident")
DimPlot(scRNA_Myeloid, reduction = "harmony", group.by = "orig.ident")
ElbowPlot(scRNA_Myeloid,reduction = 'harmony')

scRNA_Myeloid <- FindNeighbors(scRNA_Myeloid, reduction = 'harmony',dims = 1:10)
scRNA_Myeloid <- FindClusters(scRNA_Myeloid, resolution = 0.8)
scRNA_Myeloid <- RunUMAP(scRNA_Myeloid, reduction = 'harmony',dims = 1:10)
scRNA_Myeloid <- RunTSNE(scRNA_Myeloid, reduction = 'harmony',dims = 1:10)
Idents(scRNA_Myeloid) <- 'Mono_celltype'
DimPlot(scRNA_Myeloid,reduction = "tsne",label = T) 
dev.off()

saveRDS(scRNA_Myeloid, 'scRNA_Myeloid_new.rds')


#### 注释后可视化 ####
## tsne图
DimPlot(scRNA_Myeloid,reduction = "tsne",group.by = 'Mono_celltype',label = T,cols = cors) 

library(scRNAtoolVis)
pdf(file=paste0('Myeloid_celltype_tsne.pdf'),width = 7,height = 5)
clusterCornerAxes(object = scRNA_Myeloid, reduction = 'tsne',
                  clusterCol = "Mono_celltype", #分类依据
                  #noSplit = F, groupFacet = 'group', #分组时使用
                  cellLabel = T, cellLabelSize = 3.5)+  #标签
  #addCircle = TRUE, cicAlpha = 0.1, nbin = 200) + #画圈
  scale_color_npg() + scale_fill_npg()
dev.off()

##绘制细胞marker umap图
pdf(file=paste0('Myeloid_celltype_markertsne.pdf'),width = 12,height = 6)
FeatureCornerAxes(object = scRNA_Myeloid, reduction = 'tsne',
                  groupFacet = NULL,
                  relLength = 0.5,
                  relDist = 0.2,
                  aspect.ratio = 1,
                  features = c("CXCL8","APOE","FCN1",
                               "CD207","CD1B","CD83"))
dev.off()

##绘制细胞marker气泡图
Myeloid_marker <- c('CXCL8','CCL3','CCL20','CCL4', # Mφ_CCL3
                    "C1QA","C1QB",'APOE','APOC1',  # Mφ_APOE
                    'IL1B','IL1RN','FCN1','DUSP6', # Mono_FCN1
                    'CD207','CD1A', #LC1
                    'CD1C', 'CD1B','CLEC10A', #LC2
                    'CD83','MMP9','PIM3') # aLC
Myeloid_marker <- top5_markers$gene                

pdf(file=paste0('Myeloid_celltype_dotplot.pdf'),width = 6,height = 8)
DotPlot(scRNA_Myeloid, features = Myeloid_marker, assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 60, hjust=1,vjust=1))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF')) #颜色
dev.off()
#美化气泡图  10*8
jjDotPlot(object = scRNA_Myeloid,gene = Myeloid_marker,
          id = 'Mono_celltype',  #分组设置
          ytree = F,
          dot.col = c('#3C5488FF','white','#E64B35FF'),
          rescale.min = 0,rescale.max = 2,midpoint = 0) #+ coord_flip()


####髓细胞比例####
library(reshape2)
library(ggplot2)
#define the color
library(ggsci)
cors <- pal_npg()(10) #定义颜色

#准备绘图数据
prop_df <- table(scRNA_Myeloid$Mono_celltype,scRNA_Myeloid$group) %>% melt()
colnames(prop_df) <- c("Cluster","Sample","Number")
prop_df$Cluster <- factor(prop_df$Cluster)
prop_df$Proportion <- ave(prop_df$Number, prop_df$Sample, 
                          FUN = function(x) x/sum(x))
#比例图1  尺寸6*3
pdf(file=paste0('Myeloid_proportion_barplot.pdf'),width = 5.5,height = 3.5)
ggplot(data = prop_df, aes(x = Sample, y = Number, fill = Cluster)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  scale_fill_manual(values = cors) +
  theme_bw() +
  labs(x="",y="Number") +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1),
        plot.margin = margin(10, 10, 10, 10))
dev.off()

#批量散点箱式图
#计算各组中各样本不同细胞群比例
Cellratio <- prop.table(table(Idents(scRNA_Myeloid), scRNA_Myeloid$orig.ident), margin = 2)
Cellratio <- data.frame(Cellratio)
colnames(Cellratio) <- c('Mono_celltype','Sample','Freq')
Cellratio$group <- ifelse(grepl('H',Cellratio$Sample),'Healthy',
                          ifelse(grepl('M',Cellratio$Sample),'Melanoma','Vitiligo'))

my_comparisons <- list(c("Healthy","Melanoma"),c("Healthy","Vitiligo"),c("Melanoma","Vitiligo"))

pdf(file=paste0('Myeloid_proportion_boxplot.pdf'),width = 7,height = 4)
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
  facet_wrap(~Mono_celltype, ncol=3, scales = "free_y") + # 按照细胞类型分面
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.border = element_rect(colour = "black", fill=NA)) #+
#stat_compare_means(comparisons=my_comparisons,
#symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1),symbols=c("***", "**", "*", "ns")),label="p.signif") # 添加t检验
dev.off()

####髓细胞组间DEG####
## 1)pseudobulks差异分析####
bs = split(colnames(scRNA_Myeloid),scRNA_Myeloid$orig.ident)
ct = do.call(
  cbind,lapply(names(bs), function(x){ 
    kp =colnames(scRNA_Myeloid) %in% bs[[x]]
    rowSums( as.matrix(scRNA_Myeloid@assays$RNA@counts[, kp]  ))
  })
)
colnames(ct) <- names(bs)
phe = unique(scRNA_Myeloid@meta.data[,c('orig.ident','group')])#样本&信息，自行修改
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
scRNA_DC=subset(scRNA_DC,group!='Melanoma')
Idents(scRNA_DC)='group'
all_markers <- FindAllMarkers(object = scRNA_DC,
                              only.pos = T, #only.pos改为T则只输出高表达gene
                              min.pct = 0.25,logfc.threshold = 0.25)
top10_markers <- all_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top10_markers <- top10_markers[!duplicated(top10_markers$gene),]
# 各亚群平均表达量提取
genes <- unique(top10_markers$gene)
aver_dt <- AverageExpression(scRNA_DC,
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
celltype_col <- c4a('10', 2) #此处改为细胞类型数目
gene_anno <- data.frame(gene_anno = top10_markers$cluster,#行注释：Top5marker对应celltype
                        row.names = top10_markers$gene)
cell_anno <- data.frame(cell_anno = colnames(aver_dt),#列注释：celltype
                        row.names = colnames(aver_dt))
names(celltype_col) <- cell_anno$cell_anno
anno_col <- list(cell_anno = celltype_col,gene_anno = celltype_col)

pdf(file="DC_top10marker_heatmap.pdf", width=4, height=6)
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


####3)富集分析####
###GSVA分析 (也可使用AUCell，参考后面)
library(GSVA)
# 各亚群平均表达量提取
genes <- unique(rownames(scRNA_Myeloid@assays$RNA@counts))
aver_dt <- AverageExpression(scRNA_Myeloid,
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


####Mf_APOE组间DEG####
library(ggVolcano)
library(ggpubr)
library(ggthemes)
scRNA_Myeloid = subset(scRNA_Myeloid,Mono_celltype %in% c('CD8+ Trm'))
Mf_APOE_markers = all_markers[all_markers$cluster=='CD8+ Trm',]
Mf_APOE_markers <- Mf_APOE_markers %>%
  filter(p_val_adj > 0) %>%  #把p等于0的删除掉
  mutate(change = if_else(p_val_adj > 0.05, "ns", #pvalue/padj
                          if_else(abs(avg_log2FC) < 0.25, "ns",
                                  if_else(avg_log2FC >= 0.25, "up", "down"))))

Mf_APOE_markers$logP <- -log10(Mf_APOE_markers$p_val_adj)
#添加火山图的基因标签
Mf_APOE_markers$Label = ""   #新加一列label
Mf_APOE_markers <- Mf_APOE_markers[order(Mf_APOE_markers$p_val_adj), ]  #对差异基因的p值进行从小到大的排序
Mf_APOE_markers$Gene <- rownames(Mf_APOE_markers)
up.genes <- head(Mf_APOE_markers$Gene[which(Mf_APOE_markers$change == "up")], 5) #高表达的基因中选择fdr值最小的5个
down.genes <- head(Mf_APOE_markers$Gene[which(Mf_APOE_markers$change == "down")], 5) #低表达的基因中选择fdr值最小的5个
#将up.genes和down.genes合并，加入到Label中
DEG.top5.genes <- c(as.character(up.genes), as.character(down.genes))
Mf_APOE_markers$Label[match(DEG.top5.genes, Mf_APOE_markers$Gene)] <- DEG.top5.genes

pdf('stressMyeloid_marker_volcano.pdf', width = 5, height = 4)
gradual_volcano(Mf_APOE_markers, x = "avg_log2FC", y = "p_val_adj",
                log2FC_cut = 1,
                label = "Label", label_number = nrow(Mf_APOE_markers), output = FALSE)
dev.off()


####亚群marker富集分析####
####0)提取亚群marker####
all_markers <- FindAllMarkers(object = scRNA_Myeloid,
                              only.pos = T, #only.pos改为T则只输出高表达gene
                              min.pct = 0.25,logfc.threshold = 0.25) #画火山图时logfc.threshold = 0
top5_markers <- all_markers %>% group_by(cluster) %>% top_n(5, avg_log2FC) #输出每个组的前5个高的gene
top5_markers <- subset(top5_markers, !duplicated(gene))
# 各亚群平均表达量提取
genes <- unique(top5_markers$gene)
aver_dt <- AverageExpression(scRNA_Myeloid,
                             features = genes,
                             group.by = 'Mono_celltype',
                             slot = 'data')
aver_dt <- as.data.frame(aver_dt$RNA)

##热图
library(ComplexHeatmap)
library(cols4all)
#热图配色自定义
mycol <- colorRampPalette(c("#4DBBD5FF", "white", "#E64B35FF"))(50)
#行列注释配色自定义
celltype_col <- c4a('10', 6) #此处改为细胞类型数目
gene_anno <- data.frame(gene_anno = top5_markers$cluster,#行注释：Top5marker对应celltype
                        row.names = top5_markers$gene)
cell_anno <- data.frame(cell_anno = colnames(aver_dt),#列注释：celltype
                        row.names = colnames(aver_dt))
names(celltype_col) <- cell_anno$cell_anno
anno_col <- list(cell_anno = celltype_col,gene_anno = celltype_col)

pdf(file="Myeloid_markertop5_heatmap.pdf", height=5, width=4.5)
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

##简易热图
#法一：平均热图  #4.5*5
library(scRNAtoolVis)
AverageHeatmap(object = scRNA_Myeloid, 
               markerGene = top5_markers$gene,
               clusterAnnoName = F, annoCol = TRUE, myanCol = cors[1:6],
               htCol = c("#3C5488FF", "white", "#DC0000FF"))
#法二：原始热图
pdf(file=paste0('scRNA_Myeloid_DEG_heatmap.pdf'),width = 5,height = 4)
DoHeatmap(scRNA_Myeloid, top5_markers$gene, 
          group.by = 'Mono_celltype', group.colors = cors) +
  scale_fill_gradientn(colors = c("#3C5488FF", "white", "#DC0000FF"))
dev.off()


##火山图
library(ggVolcano)
Mf_APOE_markers = all_markers[all_markers$cluster=='aLC',]
Mf_APOE_markers <- Mf_APOE_markers %>%
  mutate(change = if_else(p_val_adj > 0.05, "ns", #pvalue/padj
                          if_else(abs(avg_log2FC) < 0.25, "ns",
                                  if_else(avg_log2FC >= 0.25, "up", "down"))))
Mf_APOE_markers$logP <- -log10(Mf_APOE_markers$p_val_adj)
#添加火山图的基因标签
Mf_APOE_markers$Label = ""   #新加一列label
Mf_APOE_markers <- Mf_APOE_markers[order(Mf_APOE_markers$p_val_adj), ]  #对差异基因的p值进行从小到大的排序
up.genes <- head(Mf_APOE_markers$gene[which(Mf_APOE_markers$change == "up")], 5) #高表达的基因中选择fdr值最小的5个
down.genes <- head(Mf_APOE_markers$gene[which(Mf_APOE_markers$change == "down")], 5) #低表达的基因中选择fdr值最小的5个
#将up.genes和down.genes合并，加入到Label中
DEG.top5.genes <- c(as.character(up.genes), as.character(down.genes))
Mf_APOE_markers$Label[match(DEG.top5.genes, Mf_APOE_markers$gene)] <- DEG.top5.genes
#绘图
pdf('Mf_APOE_marker_volcano.pdf', width = 6, height = 4)
gradual_volcano(Mf_APOE_markers, x = "avg_log2FC", y = "p_val_adj",
                log2FC_cut = 0.5,
                label = "Label", label_number = nrow(Mf_APOE_markers), output = FALSE)
dev.off()


####1)GO分析####
Mf_CXCL8_markers = all_markers[all_markers$cluster=='Mφ_CXCL8',]$gene
Mf_APOE_markers = all_markers[all_markers$cluster=='Mφ_APOE',]$gene
Mono_FCN1_markers = all_markers[all_markers$cluster=='Mono_FCN1',]$gene
#产生ENTREZID
genelist <- bitr(Mf_CXCL8_markers, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')

ego1 <- enrichGO(gene = genelist$ENTREZID,
                 OrgDb = org.Hs.eg.db, 
                 ont = "BP",            #此处选择BP
                 pAdjustMethod = "BH",
                 minGSSize = 1,
                 pvalueCutoff =0.05, qvalueCutoff =0.2,
                 readable = TRUE)
ego1_res <- ego1@result
ego1_res$Group <- 'Mφ_CXCL8'
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
pdf('stressMyeloid_GO_BP_barplot.pdf', width = 5, height = 2.5)
ggplot(ego1_res[1:10,],aes(Count,Description))+ #只展示前十条通路
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=LogP),stat='identity')+
  scale_fill_gradient(low="#F39B7FFF",high="#DC0000FF")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=0.5))
dev.off()

# 合并富集结果
ego_res = rbind(ego1_res[1:10,],ego2_res[1:10,],ego3_res[1:10,])

# 分组柱状图
pdf('Myeloid_GOBP_barplot.pdf', width = 6, height = 6)
ggplot(ego_res, aes(Count, Description)) +
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=LogP),stat='identity') +
  scale_fill_gradient2(low="#486b98", mid="#f5f2b1", high="#b93735", midpoint=25)+
  facet_wrap(~Group, scales="free", ncol=1)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=1))
dev.off()

# 分组气泡图
pdf('Myeloid_GOBP_dotplot.pdf', width = 6.5, height = 4)
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
genes <- unique(rownames(scRNA_Myeloid@assays$RNA@counts))
aver_dt <- AverageExpression(scRNA_Myeloid,
                             features = genes,
                             group.by = 'Mono_celltype',
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
                 'BIOCARTA_IFNG_PATHWAY',
                 'GOBP_RESPONSE_TO_OXIDATIVE_STRESS',
                 #'GOBP_T_CELL_MEDIATED_CYTOTOXICITY',
                 'GOBP_AUTOPHAGY_OF_MITOCHONDRION',
                 'GOBP_RESPONSE_TO_TYPE_I_INTERFERON',
                 'WP_FERROPTOSIS',
                 'BIOCARTA_CYTOKINE_PATHWAY')
#趋化因子
geneSet_name = c('GOBP_REGULATION_OF_T_CELL_CHEMOTAXIS',
                 'GOBP_REGULATION_OF_MACROPHAGE_CHEMOTAXIS',
                 'GOBP_REGULATION_OF_CELL_CHEMOTAXIS_TO_FIBROBLAST_GROWTH_FACTOR',
                 'GOBP_MONOCYTE_EXTRAVASATION',
                 'GOBP_T_CELL_EXTRAVASATION')

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
group_list <- levels(scRNA_Myeloid$Mono_celltype)
group_list = factor(group_list,levels = levels(scRNA_Myeloid$Mono_celltype))
#选择相应表达水平中位绝对偏差中位数排名前20位的通路（看情况运行）
mad_scores <- apply(exp_gsva, 1, mad)
top_genes <- order(mad_scores, decreasing = TRUE)[1:20]
exp_gsva_top <- exp_gsva[top_genes, ]

##热图
library(pheatmap) 
annotation_col = data.frame(group=group_list)
rownames(annotation_col) <- colnames(exp_gsva)
ann_colors = list(Cell_group = c("Mφ_CXCL8"="#003366","Mφ_APOE"="#993333","Mono_FCN1"="#00A087FF",
                                 "LC1"="#3C5488FF","LC2"="#F39B7FFF","aLC"="#B09C85FF")) 

pdf('Myeloid_celltype_GSVA_interest.pdf', width = 5.5, height = 4.5)
pheatmap::pheatmap(exp_gsva,
                   show_colnames = T, 
                   cluster_rows = F, cluster_cols = F, 
                   color = colorRampPalette(c("#4DBBD5FF", "white", "#E64B35FF"))(50),
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors)
dev.off()


####AUCell关键通路评分####
###AUCell与GSVA方法结果可能有不一致！
##通过R获取msigdbr数据集
#geneSet和GSVA一样 抗原递呈
geneSet_name = c('WP_OXIDATIVE_STRESS_RESPONSE',
                 'GOBP_RESPONSE_TO_TYPE_I_INTERFERON',
                 'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_OR_POLYSACCHARIDE_ANTIGEN_VIA_MHC_CLASS_II',
                 'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I',
                 'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_ENDOGENOUS_PEPTIDE_ANTIGEN',
                 'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_ENDOGENOUS_ANTIGEN')

cells_AUC <- AUCell_run(scRNA_Myeloid@assays$RNA@data, my_geneSet)
#提取PATHWAY
geneSet <- geneSet_name[2]   
AUCell_auc <- as.numeric(getAUC(cells_AUC)[geneSet, ])
scRNA_Myeloid$AUCell <- AUCell_auc #添加至metadata中
head(scRNA_Myeloid@meta.data)

##绘图可视化
#Seurat自带小提琴图
pdf(file='Myeloid_antigen_process_vlnplot.pdf',width = 5,height = 4)
VlnPlot(scRNA_Myeloid,features = 'AUCell', #features也可改为AUCell
        pt.size = 0,group.by = "Mono_celltype",col = cors) #按细胞类型分组
dev.off()

#箱式图
my_comparisons = list(c("Healthy","Melanoma"),c("Healthy","Vitiligo"),c("Melanoma","Vitiligo"))

pdf(file='Myeloid_antigen_process_boxplot.pdf',width = 6,height = 3.5)
ggboxplot(scRNA_Myeloid@meta.data, x="group", y="AUCell", width = 0.6, #按group分组
          color = "black",#轮廓颜色
          fill="group",#填充
          palette = cors,
          xlab = F, #不显示x轴的标签
          bxp.errorbar=T,#显示误差条
          bxp.errorbar.width=0.5, #误差条大小
          size=0.5, #箱型图边线的粗细
          outlier.shape=NA, #不显示outlier
          legend = "right") + 
  facet_wrap(~Mono_celltype,ncol=3,scales="free_y") + # 按照细胞类型分面
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.border = element_rect(colour = "black", fill=NA)) #+
  #stat_compare_means(comparisons=my_comparisons,
                     #symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1),symbols=c("***", "**", "*", "ns")),label="p.signif") # 添加t检验
dev.off()

## tsne图
# 法二  https://zhuanlan.zhihu.com/p/482523999
library(ggrepel)
#提取tsne坐标数据
tsne <- data.frame(scRNA_Myeloid@meta.data, scRNA_Myeloid@reductions$tsne@cell.embeddings)
library(ggplot2)
scRNA_Mf$tSNE_1 <- tsne$tSNE_1
scRNA_Mf$tSNE_2 <- tsne$tSNE_2
mydata <- FetchData(scRNA_Myeloid,vars = c("tSNE_1","tSNE_2","AUCell"))
#绘图
pdf(file='CD8T_CYTOTOXICITY_tsne.pdf',width = 4,height = 3)
ggplot(mydata,aes(x = tSNE_1,y = tSNE_2,colour = AUCell))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('#333366',"#6666FF",'#CC3333','#FFCC33')) + 
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"))
dev.off()


###感兴趣通路基因dotplot
#趋化因子，炎症因子
genes <- unique(my_geneSet$GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_ENDOGENOUS_ANTIGEN)

pdf(file=paste0('aLC_antigen_process_dotplot.pdf'),width = 4,height = 5)
DotPlot(scRNA_aLC, features = genes, group.by = 'group', assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF')) #颜色
dev.off()


####巨噬细胞M1/M2样评分####
scRNA_Mf <- subset(scRNA_Myeloid, celltype=='Mono phagocyte')

# 根据基因集评分分组
M1_features <- list(c('IL1A','IL1B','IL6','NOS2','TLR2','TLR4','CD80','CD86')) #需要真是的基因symbol
M2_features <- list(c('CD115','CD206','PPARG','ARG1','CD163','CD301','IL4R','PDL2','CD200R','MGL1','MGL2'))

### AddModuleScore 打分
scRNA_Mf <- AddModuleScore(scRNA_Mf,
                             features = M2_features,
                             ctrl = 100,
                             name = "M2_features")
head(scRNA_Mf@meta.data)
#这里就得到了基因集评分结果
colnames(scRNA_Mf@meta.data)[14] <- 'M2_Score'

## tsne图


#### 氧化应激评分分群再分析 ####
##原因：之前的传统分群找不出白癜风和健康组的差异，因此换个角度分群再分析！
#准备基因集
geneSet = msigdbr(species = "Homo sapiens")  #category = "C2"
geneSet %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
geneSet_name = 'WP_OXIDATIVE_STRESS_RESPONSE'
my_geneSet = subset(geneSet, gs_name %in% geneSet_name)
my_geneSet = my_geneSet %>% split(x = .$gene_symbol, f = .$gs_name)
#AUCell评分
cells_AUC <- AUCell_run(scRNA_Myeloid@assays$RNA@data, my_geneSet)
geneSet <- "WP_OXIDATIVE_STRESS_RESPONSE"   
AUCell_auc <- as.numeric(getAUC(cells_AUC)[geneSet, ])
scRNA_Myeloid$AUCell <- AUCell_auc 
##绘图可视化
#Seurat自带小提琴图、箱式图、tsne图同理绘制

##根据氧化应激评分分群(选择分高的Fib1,Fib2,Fib4)
Myeloid_celltype=c('Stressed Mφ','Stressed Mφ','Mono_FCN1','LC1','LC2','aLC')
Idents(scRNA_Myeloid) <- scRNA_Myeloid$Mono_celltype
names(Myeloid_celltype) <- levels(scRNA_Myeloid)
scRNA_Myeloid <- RenameIdents(scRNA_Myeloid, Myeloid_celltype)
scRNA_Myeloid$Myeloid_celltype1 <- Idents(scRNA_Myeloid)
#自定义celltype顺序
cell = c('Stressed Mφ','Mono_FCN1','LC1','LC2','aLC') 
scRNA_Myeloid$Myeloid_celltype1 <- factor(scRNA_Myeloid$Myeloid_celltype1,levels = cell)
Idents(scRNA_Myeloid) <- scRNA_Myeloid$Myeloid_celltype1
saveRDS(scRNA_Myeloid, 'scRNA_Myeloid.rds')

##再次可视化

##再次绘制细胞比例

##再次亚群marker&富集分析
Mf_stress_markers = FindMarkers(scRNA_Myeloid,ident.1 = 'Stressed Mφ', 
                              verbose = FALSE, test.use = 'wilcox',
                              min.pct = 0.25, logfc.threshold = 0)

##再次拟时序分析、细胞通讯、转录因子分析