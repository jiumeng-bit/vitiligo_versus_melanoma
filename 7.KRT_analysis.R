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
scRNA_KRT=readRDS('./scRNA_KRT.rds')
table(scRNA_KRT$KRT_celltype,scRNA_KRT$group)

#### 提取角质细胞 ####
#https://zhuanlan.zhihu.com/p/375318689 朗格汉斯细胞只是DC的一个亚群！
scRNA_KRT=subset(scRNA,celltype=='Keratinocyte') 
# 提细胞亚群,重新降维聚类
scRNA_KRT <- FindVariableFeatures(scRNA_KRT, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNA_KRT)
scRNA_KRT <- ScaleData(scRNA_KRT, features = scale.genes)
scRNA_KRT<- RunPCA(scRNA_KRT, features = VariableFeatures(scRNA_KRT))
DimPlot(scRNA_KRT, reduction = "pca", group.by = "orig.ident")
ElbowPlot(scRNA_KRT)

## 重新harmony
library(harmony)
set.seed(1000)
scRNA_KRT <- RunHarmony(scRNA_KRT, group.by.vars = "orig.ident")
DimPlot(scRNA_KRT, reduction = "harmony", group.by = "orig.ident")
ElbowPlot(scRNA_KRT,reduction = 'harmony')

scRNA_KRT <- FindNeighbors(scRNA_KRT, reduction = 'harmony',dims = 1:10)
scRNA_KRT <- FindClusters(scRNA_KRT, resolution = 0.1)
scRNA_KRT <- RunUMAP(scRNA_KRT, reduction = 'harmony',dims = 1:10)
scRNA_KRT <- RunTSNE(scRNA_KRT, reduction = 'harmony',dims = 1:10)
DimPlot(scRNA_KRT,reduction = 'tsne',split.by = 'group')
DimPlot(scRNA_KRT,reduction = "tsne",group.by = "seurat_clusters",label = T) 

## 角质细胞亚群注释参考
# marker：https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9310536/#sd
Idents(scRNA_KRT)=scRNA_KRT$seurat_clusters
KRT_marker <- c('KRT16','KRT6A','KRT6B','S100A8','S100A9') #Stressed KRT

DotPlot(scRNA_KRT, 
        features = KRT_marker,
        group.by = "seurat_clusters") + coord_flip()

scRNA_KRT <- AddModuleScore(scRNA_KRT,   # AddModuleScore 打分
                             features = list(KRT_marker),
                             ctrl = 100,
                             name = "KRT_features")
head(scRNA_KRT@meta.data)
colnames(scRNA_KRT@meta.data)[15] <- 'Stress_Score'
VlnPlot(scRNA_KRT,features = 'Stress_Score', 
        pt.size = 0,group.by = "seurat_clusters",col = cors)

KRT_celltype=c('KRT',
               'KRT','Stressed KRT','Stressed KRT','KRT','KRT',
               'KRT','KRT','KRT','KRT')

Idents(scRNA_KRT) <- scRNA_KRT@meta.data$seurat_clusters
names(KRT_celltype) <- levels(scRNA_KRT)
scRNA_KRT <- RenameIdents(scRNA_KRT, KRT_celltype)
scRNA_KRT$KRT_celltype <- Idents(scRNA_KRT)
# 删除unknown
scRNA_KRT <- subset(scRNA_KRT, KRT_celltype != 'Unknown')
scRNA_KRT$KRT_celltype <- Idents(scRNA_KRT) #完全去除Unknown
table(scRNA_KRT$KRT_celltype)
#自定义celltype顺序
cell = c('KRT','Stressed KRT') 
scRNA_KRT$KRT_celltype <- factor(scRNA_KRT$KRT_celltype,levels = cell)
Idents(scRNA_KRT) <- scRNA_KRT$KRT_celltype
saveRDS(scRNA_KRT, 'scRNA_KRT.rds')



#### 注释后可视化 ####
## tsne图
DimPlot(scRNA_KRT,reduction = "tsne",group.by = 'KRT_celltype',label = T,cols = cors) 

library(scRNAtoolVis)
pdf(file=paste0('scRNA_KRT_celltype_tsne.pdf'),width = 6,height = 4)
clusterCornerAxes(object = scRNA_KRT, reduction = 'tsne',
                  clusterCol = "KRT_celltype", #分类依据
                  #noSplit = F, groupFacet = 'group', #分组时使用
                  cellLabel = T, cellLabelSize = 3.5)+  #标签
  #addCircle = TRUE, cicAlpha = 0.1, nbin = 200) + #画圈
  scale_color_npg() + scale_fill_npg()
dev.off()


##绘制细胞marker气泡图
KRT_marker <- c('KRT16','KRT6A','KRT6B','S100A8','S100A9', #Stressed KRT 
                top10_markers$gene[2:5])            

pdf(file=paste0('scRNA_KRT_celltype_dotplot.pdf'),width = 4,height = 4)
DotPlot(scRNA_KRT, features = KRT_marker, assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust=1,vjust=1))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF')) #颜色
dev.off()
#美化气泡图  6*4
jjDotPlot(object = scRNA_KRT,gene = KRT_marker,
          id = 'KRT_celltype',  #分组设置
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
prop_df <- table(scRNA_KRT$KRT_celltype,scRNA_KRT$group) %>% melt()
colnames(prop_df) <- c("Cluster","Sample","Number")
prop_df$Cluster <- factor(prop_df$Cluster)
prop_df$Proportion <- ave(prop_df$Number, prop_df$Sample, 
                          FUN = function(x) x/sum(x))
#比例图1  尺寸6*3
pdf(file=paste0('scRNA_KRT_proportion_barplot.pdf'),width = 4,height = 3)
ggplot(data = prop_df, aes(x =Sample, y = Number, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8, position="fill")+
  scale_fill_manual(values=cors) + #配色
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45))
dev.off()


####角质细胞组间DEG####
## 1)pseudobulks差异分析####
bs = split(colnames(scRNA_KRT),scRNA_KRT$orig.ident)
ct = do.call(
  cbind,lapply(names(bs), function(x){ 
    kp =colnames(scRNA_KRT) %in% bs[[x]]
    rowSums( as.matrix(scRNA_KRT@assays$RNA@counts[, kp]  ))
  })
)
colnames(ct) <- names(bs)
phe = unique(scRNA_KRT@meta.data[,c('orig.ident','group')])#样本&信息，自行修改
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
write.csv(DEG_DESeq2,"scRNA_KRT_DEG_H_vs_V.csv")

## 2)FindMarkers组间DEG####
Idents(scRNA_KRT)='group'
all_markers <- FindAllMarkers(object = scRNA_KRT,
                              only.pos = T, #only.pos改为T则只输出高表达gene
                              min.pct = 0.25,logfc.threshold = 0.25)
top10_markers <- all_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top10_markers <- top10_markers[!duplicated(top10_markers$gene),]
# 各亚群平均表达量提取
genes <- unique(top10_markers$gene)
aver_dt <- AverageExpression(scRNA_KRT,
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

pdf(file="scRNA_KRT_top10marker_heatmap.pdf", width=4, height=6)
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


####3)GSVA####
library(GSVA)
# 各亚群平均表达量提取
genes <- unique(rownames(scRNA_KRT@assays$RNA@counts))
aver_dt <- AverageExpression(scRNA_KRT,
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

pdf('scRNA_KRT_group_GSVA_interest.pdf', width = 5, height = 4.5)
pheatmap::pheatmap(exp_gsva,
                   show_colnames = T, 
                   cluster_rows = F, cluster_cols = F, 
                   color = colorRampPalette(c("#3C5488FF", "white", "#DC0000FF"))(50),
                   annotation_col = annotation_col)
dev.off()


####4)取趋势相反DEG交集####
#变化趋势相反的DEG
Idents(scRNA_KRT)='group'
DEG_Viti <- FindMarkers(object = scRNA_KRT,ident.1 = 'Vitiligo',ident.2 = 'Healthy',
                        only.pos = F, #only.pos改为T则只输出高表达gene
                        min.pct = 0.25,logfc.threshold = 0.25)
DEG_Mel <- FindMarkers(object = scRNA_KRT,ident.1 = 'Melanoma',ident.2 = 'Healthy',
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
write.csv(inverse_DEG, 'scRNA_KRT_inverse_DEG.csv',row.names = F, col.names = F)

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
pdf(file=paste0('scRNA_KRT_inverseDEG_dotplot.pdf'),width = 4.5,height = 5)
DotPlot(scRNA_KRT, features = inverse_DEG, group.by = 'group',assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colors = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF')) #颜色
dev.off()
#关注ZNF683和ATP5I


####mCAF组间DEG####
library(ggVolcano)
library(ggpubr)
library(ggthemes)
scRNA_KRTrm = subset(scRNA_KRT,KRT_celltype %in% c('CD8+ Trm'))
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

pdf('scRNA_KRTrm_marker_volcano.pdf', width = 5, height = 4)
gradual_volcano(Trm_markers, x = "avg_log2FC", y = "p_val_adj",
                log2FC_cut = 0.5,
                label = "Label", label_number = nrow(Trm_markers), output = FALSE)
dev.off()


####亚群marker富集分析####
####0)提取亚群marker####
all_markers <- FindAllMarkers(object = scRNA_KRT,
                              only.pos = T, #only.pos改为T则只输出高表达gene
                              min.pct = 0.25,logfc.threshold = 0.25) #画火山图时logfc.threshold = 0
top10_markers <- all_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC) #输出每个组的前5个高的gene
top10_markers <- subset(top10_markers, !duplicated(gene))

##热图
#法一：平均热图
library(scRNAtoolVis)
AverageHeatmap(object = scRNA_KRT,
               markerGene = top10_markers$gene,
               htCol = c("#3C5488FF", "white", "#DC0000FF"))
#法二：原始热图
pdf(file=paste0('scRNA_KRT_DEG_heatmap.pdf'),width = 5,height = 4)
DoHeatmap(scRNA_KRT, top10_markers$gene, 
          group.by = 'KRT_celltype', group.colors = cors) +
  scale_fill_gradientn(colors = c("#3C5488FF", "white", "#DC0000FF"))
dev.off()

##火山图
mCAF_markers = all_markers[all_markers$cluster=='mCAF',]
mCAF_markers <- mCAF_markers %>% 
  filter(p_val_adj > 0) %>%  #把p等于0的删除掉
  mutate(change = if_else(p_val_adj > 0.05, "ns", #pvalue/padj
                          if_else(abs(avg_log2FC) < 0.25, "ns",
                                  if_else(avg_log2FC >= 0.25, "up", "down")))) 
mCAF_markers$logP <- -log10(mCAF_markers$p_val_adj)
#添加火山图的基因标签
mCAF_markers$Label = ""   #新加一列label
mCAF_markers <- mCAF_markers %>% arrange(avg_log2FC)  #对差异基因的logfc进行从小到大的排序
up.genes <- tail(mCAF_markers$gene, 5)   #高表达的基因
down.genes <- head(mCAF_markers$gene, 5) #低表达的基因
#将up.genes和down.genes合并，加入到Label中
DEG.top5.genes <- c(as.character(up.genes), as.character(down.genes))
mCAF_markers$Label[match(DEG.top5.genes, mCAF_markers$gene)] <- DEG.top5.genes
#绘图
pdf('scRNA_mCAF_marker_volcano.pdf', width = 5, height = 4)
gradual_volcano(mCAF_markers, x = "avg_log2FC", y = "p_val_adj",
                log2FC_cut = 0.25,
                label = "Label", label_number = nrow(mCAF_markers), output = FALSE)#+xlim(-4, 4)
dev.off()


####1)GO分析####
KRT_markers = all_markers[all_markers$cluster %in% c('KRT'),]$gene
Stress_markers = all_markers[all_markers$cluster %in% c('Stressed KRT'),]$gene

#产生ENTREZID
genelist <- bitr(KRT_markers, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')

ego1 <- enrichGO(gene = genelist$ENTREZID,
                 OrgDb = org.Hs.eg.db, 
                 ont = "BP",            #此处选择BP
                 pAdjustMethod = "BH",
                 minGSSize = 1,
                 pvalueCutoff =0.05, qvalueCutoff =0.2,
                 readable = TRUE)
ego1_res <- ego1@result
ego1_res$Group <- 'KRT'
ego1_res$LogP <- -log(ego1_res$p.adjust) #计算-logP
##ego2_res,ego3_res同理计算

#只提取与xxx有关的基因集   pathway|signal|immune|oxidat|interferon|inflam|ediper
ego1_res = ego1_res  %>%
  dplyr::filter(grepl('pathway|signal|immune|oxidat|interferon|inflam|ediper', Description)) %>%
  dplyr::arrange(dplyr::desc(LogP),dplyr::desc(Description)) %>%
  mutate(Description =forcats:: fct_inorder(Description))

## 气泡图
barplot(ego1, showCategory = 10)
# 个性化柱状图
pdf('scRNA_KRTem_GO_BP_barplot.pdf', width = 6, height = 5)
ggplot(ego1_res[1:10,],aes(Count,Description))+ #只展示前十条通路
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=LogP),stat='identity')+
  scale_fill_gradient(low="#F39B7FFF",high="#DC0000FF")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=0.5))
dev.off()

# 合并富集结果
ego_res = rbind(ego1_res0[c(1,3:6,9:13),],ego2_res0[c(1:4,7:11,15),])

# 分组柱状图
pdf('scRNA_KRT_GOBP_barplot.pdf', width = 6, height = 6)
ggplot(ego_res, aes(Count, Description)) +
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=LogP),stat='identity') +
  scale_fill_gradient(low="#F39B7FFF",high="#DC0000FF")+
  facet_wrap(~Group, scales="free", ncol=1)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=1))
dev.off()

# 分组气泡图
pdf('scRNA_KRT_GOBP_dotplot.pdf', width = 6.5, height = 4)
ggplot(ego_res, aes(Group, factor(Description, levels=ego_res$Description))) +
  geom_point(aes(color=LogP, size=Count))+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=1))+
  scale_color_gradient(low='#F39B7FFF',high='#DC0000FF')+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1))
dev.off()



####2)GSVA分析####
library(GSVA)
# 各亚群平均表达量提取
genes <- unique(rownames(scRNA_KRT@assays$RNA@counts))
aver_dt <- AverageExpression(scRNA_KRT,
                             features = genes,
                             group.by = 'KRT_celltype',
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
group_list <- levels(scRNA_KRT$KRT_celltype)
group_list = factor(group_list,levels = levels(scRNA_KRT$KRT_celltype))

##气泡热图
library(dplyr)
library(reshape)
data <- data.frame(exp_gsva) %>% 
  mutate(gene = rownames(data.frame(exp_gsva)))
data <- melt(data,id.vars='gene')
colnames(data) <- c('gene','group','value')
head(data)
## plot
library(ggplot2)
p1 <- ggplot(data,aes(x=group ,y= gene)) + 
  geom_point(aes(size=-log10(value),fill=value),shape=21,color="black") +
  scale_fill_gradient2(name = 'value\n(Expression)', #灵活调整颜色区间
                       limit = c(0,1.001),breaks = c(0,0.5,1.0),
                       low='#4DBBD5FF',high='#DC0000FF',mid="#F39B7FFF",midpoint = 0.5)+
  scale_size_continuous(name = 'value',
                        limit = c(-1,3.1), #c(-0.001,3.1)
                        breaks = c(0,1,2,3))+
  #geom_hline(yintercept=c(5.5, 10.5))+  #图中添加横线
  labs(x=NULL,y=NULL)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 10, color = "black"),
        axis.text.y = element_blank(),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1))
p1  

anotate <- data %>% distinct(gene,.keep_all = T)
p2 <- ggplot(anotate, aes(x = 0,y = gene,label= gene)) +
  geom_text() + theme_void()
p2
# patch
library(patchwork)
library(cowplot)
p2+p1#+plot_layout(nrow= 1,width = c(1, 2))
ggsave('KRT_celltype_GSVA_interest.pdf',width = 6,height = 5)



####AUCell关键通路评分####
###AUCell与GSVA方法结果可能有不一致！
#角质细胞关键geneSet: "IFN"

cells_AUC <- AUCell_run(scRNA_KRT@assays$RNA@data, my_geneSet)
#提取PATHWAY
geneSet <- "KEGG_JAK_STAT_SIGNALING_PATHWAY"   
AUCell_auc <- as.numeric(getAUC(cells_AUC)[geneSet, ])
scRNA_KRT$AUCell <- AUCell_auc #添加至metadata中
head(scRNA_KRT@meta.data)

##绘图可视化
#Seurat自带小提琴图
pdf(file='scRNA_KRT_IFNG_vlnplot.pdf',width = 5,height = 4)
VlnPlot(scRNA_KRT,features = 'AUCell', #features也可改为AUCell
        pt.size = 0,group.by = "KRT_celltype",col = cors) #按细胞类型分组
dev.off()

#箱式图
my_comparisons = list(c("Healthy","Melanoma"),c("Healthy","Vitiligo"),c("Melanoma","Vitiligo"))

pdf(file='scRNA_KRT_JAK_boxplot.pdf',width = 4,height = 3)
ggboxplot(scRNA_KRT@meta.data, x="group", y="AUCell", width = 0.6, #按group分组
          color = "black",#轮廓颜色
          fill="group",#填充
          palette = cors,
          xlab = F, #不显示x轴的标签
          bxp.errorbar=T,#显示误差条
          bxp.errorbar.width=0.5, #误差条大小
          size=0.5, #箱型图边线的粗细
          outlier.shape=NA, #不显示outlier
          legend = "right") + 
  #ylim(0,0.1) +
  #facet_wrap(~KRT_celltype,ncol=4,scales="free_y") + # 按照细胞类型分面
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  stat_compare_means(comparisons=my_comparisons,
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1),symbols=c("***", "**", "*", "ns")),label="p.signif") # 添加t检验
dev.off()

## tsne图
# 法二  https://zhuanlan.zhihu.com/p/482523999
library(ggrepel)
#提取tsne坐标数据
tsne <- data.frame(scRNA_KRT@meta.data, scRNA_KRT@reductions$tsne@cell.embeddings)
library(ggplot2)
scRNA_KRT$tsne_1 <- tsne$tSNE_1
scRNA_KRT$tsne_2 <- tsne$tSNE_2
mydata <- FetchData(scRNA_KRT,vars = c("tsne_1","tsne_2","AUCell"))
#绘图
pdf(file='scRNA_KRT_IFNG_tsne.pdf',width = 4,height = 3)
ggplot(mydata,aes(x = tsne_1,y =tsne_2,colour = AUCell))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('#333366',"#6666FF",'#CC3333','#FFCC33')) + 
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"))
dev.off()


###感兴趣通路基因dotplot
#趋化因子，炎症因子
genes <- unique(my_geneSet$BIOCARTA_IFNG_PATHWAY)

pdf(file=paste0('scRNA_KRT_IFNG_dotplot.pdf'),width = 5,height = 2.5)
DotPlot(scRNA_KRT, features = genes, group.by = 'KRT_celltype', assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF')) #颜色
dev.off()


#### KRT相关基因表达比较 ####
geneset = c('TNF','IL8','IL13','IL18','IL19','IL20','MCP1',
            'IFNG','CXCL9','CXCL10','JAK3','STAT3')#IFNG相关

DotPlot(scRNA_KRT, features = geneset, group.by = 'group', assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust =1,vjust=1))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF'))

#美化气泡图  9*7
jjDotPlot(object = scRNA_KRT,gene = geneset,
          id = 'group',  #分组设置
          ytree = F,
          dot.col = c('#3C5488FF','white','#E64B35FF'),
          rescale.min = 0,rescale.max = 2,midpoint = 0) #+ coord_flip()
ggsave('scRNA_KRT_related_dotplot.pdf',width = 8,height = 6)

##提取stressKRT组间分析
stressKRT <- subset(scRNA_KRT, KRT_celltype == 'Stressed KRT')
#组间GSVA、组间看IFNG