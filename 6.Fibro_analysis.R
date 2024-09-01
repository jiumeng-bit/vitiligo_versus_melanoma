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
setwd("E:/【科研学习】/【皮肤科】/卡波西肉瘤课题")
scRNA_Fib=readRDS('./scRNA_Fib.rds')
table(scRNA_Fib$Fib_celltype,scRNA_Fib$group)

#### 提取成纤维细胞 ####
#https://zhuanlan.zhihu.com/p/375318689 朗格汉斯细胞只是DC的一个亚群！
scRNA_Fib=subset(scRNA,celltype=='Fibroblast') 
# 提细胞亚群,重新降维聚类
scRNA_Fib <- FindVariableFeatures(scRNA_Fib, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNA_Fib)
scRNA_Fib <- ScaleData(scRNA_Fib, features = scale.genes)
scRNA_Fib<- RunPCA(scRNA_Fib, features = VariableFeatures(scRNA_Fib))
DimPlot(scRNA_Fib, reduction = "pca", group.by = "orig.ident")
ElbowPlot(scRNA_Fib)

## 重新harmony
library(harmony)
set.seed(1000)
scRNA_Fib <- RunHarmony(scRNA_Fib, group.by.vars = "orig.ident")
DimPlot(scRNA_Fib, reduction = "harmony", group.by = "orig.ident")
ElbowPlot(scRNA_Fib,reduction = 'harmony')

scRNA_Fib <- FindNeighbors(scRNA_Fib, reduction = 'harmony',dims = 1:10)
sce_res <- scRNA_Fib
for (i in c(0.01, 0.05, 0.1, 0.15, 0.2, 0.3,0.4, 0.5,0.8,1)){
  sce_res <- FindClusters(sce_res,resolution = i)
}
scRNA_Fib <- FindClusters(scRNA_Fib, resolution = 0.8)
scRNA_Fib <- RunUMAP(scRNA_Fib, reduction = 'harmony',dims = 1:10)
scRNA_Fib <- RunTSNE(scRNA_Fib, reduction = 'harmony',dims = 1:10)
DimPlot(scRNA_Fib,reduction = 'tsne',group.by = 'group')
DimPlot(scRNA_Fib,reduction = "tsne",group.by = "seurat_clusters",label = T) 

## 成纤维细胞亚群注释参考
# marker：https://mp.weixin.qq.com/s/HPkhAo1r7PT4zJW1AOjPvA
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9882988/
# 通过group分组确定哪些cluster是正常Fib还是CAF
Idents(scRNA_Fib)=scRNA_Fib$seurat_clusters
Fib_marker <- c("WISP2","SLPI", "CTHRC1",'MFAP5','TSPAN8', #Secretory-reticular
                  "APCDD1", "ID1", "WIF1",'COL18A1','PTGDS', #Secretory-papillary
                  "CCL19", "APOE","CXCL2",'CXCL3','EFEMP1', #Pro-inflammatory
                  'ASPN',"POSTN",'GPC3','TNN','SFRP1') #Mesenchymal 

Fib_marker <- c("COL1A1","SFRP2","TWIST2","PDGFRB","ACTA2", #Fib
                'PDGFRA','DCN','LUM','COL6A3','POSTN', #mCAF
                'MCAM','RGS5','DSTN','ADIRF', #vCAF
                'APOC1','APOC2','APOC3','APOA1',  #lpCAF
                'CD74','CCL5','HLADRB1') #apCAF

DotPlot(scRNA_Fib, 
        features = Fib_marker,
        group.by = "seurat_clusters") + coord_flip()

Fib_celltype=c('Fib1',
               'Fib2','mCAF','Fib3','Fib1','Fib3',
               'vCAF','Fib1','Fib4','Fib5','mCAF')

Idents(scRNA_Fib) <- scRNA_Fib@meta.data$seurat_clusters
names(Fib_celltype) <- levels(scRNA_Fib)
scRNA_Fib <- RenameIdents(scRNA_Fib, Fib_celltype)
scRNA_Fib$Fib_celltype <- Idents(scRNA_Fib)
# 删除unknown
scRNA_Fib <- subset(scRNA_Fib, Fib_celltype != 'Unknown')
scRNA_Fib$Fib_celltype <- Idents(scRNA_Fib) #完全去除Unknown
table(scRNA_Fib$Fib_celltype)
#自定义celltype顺序
cell = c('Fib1','Fib2','Fib3','Fib4','Fib5','mCAF','vCAF') 
scRNA_Fib$Fib_celltype <- factor(scRNA_Fib$Fib_celltype,levels = cell)
Idents(scRNA_Fib) <- scRNA_Fib$Fib_celltype
saveRDS(scRNA_Fib, 'scRNA_Fib.rds')



#### 注释后可视化 ####
## tsne图
DimPlot(scRNA_Fib,reduction = "tsne",group.by = 'Fib_celltype',label = T,cols = cors) 

library(scRNAtoolVis)
pdf(file=paste0('scRNA_Fib_celltype_tsne.pdf'),width = 6,height = 4)
clusterCornerAxes(object = scRNA_Fib, reduction = 'tsne',
                  clusterCol = "Fib_celltype", #分类依据
                  #noSplit = F, groupFacet = 'group', #分组时使用
                  cellLabel = T, cellLabelSize = 3.5)+  #标签
  #addCircle = TRUE, cicAlpha = 0.1, nbin = 200) + #画圈
  scale_color_npg() + scale_fill_npg()
dev.off()


##绘制细胞marker气泡图
Fib_marker <- c("COL1A1","SFRP2","TWIST2","PDGFRB", #Fib
                'PDGFRA','LUM','COL6A3','POSTN', #mCAF
                'DSTN','ADIRF') #vCAF
Fib_marker <- top5_markers$gene   #先计算all_markers             

pdf(file=paste0('scRNA_Fib_celltype_dotplot.pdf'),width = 6,height = 8)
DotPlot(scRNA_Fib, features = Fib_marker, assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 60, hjust=1,vjust=1))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF')) #颜色
dev.off()
#美化气泡图  10*8
jjDotPlot(object = scRNA_Fib,gene = Fib_marker,
          id = 'Fib_celltype',  #分组设置
          ytree = F,
          dot.col = c('#3C5488FF','white','#E64B35FF'),
          rescale.min = 0,rescale.max = 2,midpoint = 0) #+ coord_flip()


####细胞比例####
library(reshape2)
library(ggplot2)
#define the color
library(ggsci)
cors <- pal_npg()(10) #定义颜色

#准备绘图数据
prop_df <- table(scRNA_Fib$Fib_celltype,scRNA_Fib$group) %>% melt()
colnames(prop_df) <- c("Cluster","Sample","Number")
prop_df$Cluster <- factor(prop_df$Cluster)
prop_df$Proportion <- ave(prop_df$Number, prop_df$Sample, 
                          FUN = function(x) x/sum(x))
#比例图1  尺寸6*3
pdf(file=paste0('scRNA_Fib_proportion_barplot.pdf'),width = 6,height = 3)
ggplot(data = prop_df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=cors) + #配色
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Proportion")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45))
dev.off()


####成纤维细胞组间DEG####
## 1)pseudobulks差异分析####
bs = split(colnames(scRNA_Fib),scRNA_Fib$orig.ident)
ct = do.call(
  cbind,lapply(names(bs), function(x){ 
    kp =colnames(scRNA_Fib) %in% bs[[x]]
    rowSums( as.matrix(scRNA_Fib@assays$RNA@counts[, kp]  ))
  })
)
colnames(ct) <- names(bs)
phe = unique(scRNA_Fib@meta.data[,c('orig.ident','group')])#样本&信息，自行修改
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
write.csv(DEG_DESeq2,"scRNA_Fib_DEG_H_vs_V.csv")

## 2)FindMarkers组间DEG####
Idents(scRNA_Fib)='group'
all_markers <- FindAllMarkers(object = scRNA_Fib,
                              only.pos = T, #only.pos改为T则只输出高表达gene
                              min.pct = 0.25,logfc.threshold = 0.25)
top10_markers <- all_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top10_markers <- top10_markers[!duplicated(top10_markers$gene),]
# 各亚群平均表达量提取
genes <- unique(top10_markers$gene)
aver_dt <- AverageExpression(scRNA_Fib,
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

pdf(file="scRNA_Fib_top10marker_heatmap.pdf", width=4, height=6)
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
genes <- unique(rownames(scRNA_Fib@assays$RNA@counts))
aver_dt <- AverageExpression(scRNA_Fib,
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

pdf('scRNA_Fib_group_GSVA_interest.pdf', width = 5, height = 4.5)
pheatmap::pheatmap(exp_gsva,
                   show_colnames = T, 
                   cluster_rows = F, cluster_cols = F, 
                   color = colorRampPalette(c("#3C5488FF", "white", "#DC0000FF"))(50),
                   annotation_col = annotation_col)
dev.off()


####4)GSEA####
##准备分析用的geneSet
geneSet = msigdbr(species = "Homo sapiens",category = "C2") #category = "C2", subcategory = 'KEGG'
geneSet %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
geneSet = geneSet %>% dplyr::select(gs_name, gene_symbol)
# 获取GeneList
Idents(scRNA_Fib)='group'
all_markers <- FindAllMarkers(object = scRNA_Fib,
                              only.pos = F, #only.pos改为T则只输出高表达gene
                              min.pct = 0.25,logfc.threshold = 0)
Fib_markers = all_markers[all_markers$cluster=='Melanoma',]
geneList <- Fib_markers$avg_log2FC            
names(geneList) <- Fib_markers$gene         # 对GeneList命名
geneList <- sort(geneList, decreasing = T)   # 从高到低排序
# 开始GSEA富集分析
GSEA_enrichment <- GSEA(geneList,                 # 排序后的gene
                        TERM2GENE = geneSet,      # 基因集
                        pvalueCutoff = 0.1,      # P值阈值
                        pAdjustMethod = "BH")     # 校正P值的计算方法
result <- data.frame(GSEA_enrichment)
#只提取与xxx有关的基因集   HYPOXIA|WNT
result1 = result  %>%
  dplyr::filter(grepl('HYPO', Description))
# plot
library(enrichplot)
geneset_plot <- c("GROSS_HYPOXIA_VIA_ELK3_AND_HIF1A_DN")
pdf('scRNA_CAF_HYPOXIA_GSEA.pdf', width = 5, height = 4)
gseaplot2(GSEA_enrichment,
          geneSetID = geneset_plot,
          color = cors,
          title = 'Hypoxia',
          rel_heights = c(1.3, 0.3, 0.6),
          pvalue_table = F)
dev.off()


####mCAF组间DEG####
library(ggVolcano)
library(ggpubr)
library(ggthemes)
Idents(scRNA_Fib)='group'
all_markers <- FindAllMarkers(object = scRNA_Fib,
                              only.pos = F, #only.pos改为T则只输出高表达gene
                              min.pct = 0.25,logfc.threshold = 0)
Fib_markers = all_markers[all_markers$cluster=='mCAF',] #组间比较的all_markers
Fib_markers <- Fib_markers %>%
  filter(p_val_adj > 0) %>%  #把p等于0的删除掉
  mutate(change = if_else(p_val_adj > 0.05, "ns", #pvalue/padj
                          if_else(abs(avg_log2FC) < 0.25, "ns",
                                  if_else(avg_log2FC >= 0.25, "up", "down"))))

Fib_markers$logP <- -log10(Fib_markers$p_val_adj)
#添加火山图的基因标签
Fib_markers$Label = ""   #新加一列label
Fib_markers <- Fib_markers[order(Fib_markers$p_val_adj), ]  #对差异基因的p值进行从小到大的排序
Fib_markers$Gene <- rownames(Fib_markers)
up.genes <- head(Fib_markers$Gene[which(Fib_markers$change == "up")], 5) #高表达的基因中选择fdr值最小的5个
down.genes <- head(Fib_markers$Gene[which(Fib_markers$change == "down")], 5) #低表达的基因中选择fdr值最小的5个
#将up.genes和down.genes合并，加入到Label中
DEG.top5.genes <- c(as.character(up.genes), as.character(down.genes))
Fib_markers$Label[match(DEG.top5.genes, Fib_markers$Gene)] <- DEG.top5.genes

pdf('scRNA_Fibrm_marker_volcano.pdf', width = 5, height = 4)
gradual_volcano(Fib_markers, x = "avg_log2FC", y = "p_val_adj",
                log2FC_cut = 0.5,
                label = "Label", label_number = nrow(Fib_markers), output = FALSE)
dev.off()


####亚群marker富集分析####
####0)提取亚群marker####
all_markers <- FindAllMarkers(object = scRNA_Fib,
                              only.pos = T, #only.pos改为T则只输出高表达gene
                              min.pct = 0.25,logfc.threshold = 0.25) #画火山图时logfc.threshold = 0
write.csv(all_markers, 'Fib_markers.csv')
top5_markers <- all_markers %>% group_by(cluster) %>% top_n(5, avg_log2FC) #输出每个组的前5个高的gene
top5_markers <- subset(top5_markers, !duplicated(gene))
# 各亚群平均表达量提取
genes <- unique(top5_markers$gene)
aver_dt <- AverageExpression(scRNA_Fib,
                             features = genes,
                             group.by = 'Fib_celltype',
                             slot = 'data')
aver_dt <- as.data.frame(aver_dt$RNA)

##热图
library(ComplexHeatmap)
library(cols4all)
#热图配色自定义
mycol <- colorRampPalette(c("#4DBBD5FF", "white", "#E64B35FF"))(50)
#行列注释配色自定义
celltype_col <- c4a('10', 7) #此处改为细胞类型数目
gene_anno <- data.frame(gene_anno = top5_markers$cluster,#行注释：Top5marker对应celltype
                        row.names = top5_markers$gene)
cell_anno <- data.frame(cell_anno = colnames(aver_dt),#列注释：celltype
                        row.names = colnames(aver_dt))
names(celltype_col) <- cell_anno$cell_anno
anno_col <- list(cell_anno = celltype_col,gene_anno = celltype_col)

pdf(file="scRNA_Fib_celltype_top5marker_heatmap.pdf", height=6, width=4.5)
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
all_markers <- FindAllMarkers(object = scRNA_Fib,
                              only.pos = F, #only.pos改为T则只输出高表达gene
                              min.pct = 0.25,logfc.threshold = 0)
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
Fib_markers = all_markers[!all_markers$cluster %in% c('mCAF','vCAF'),]$gene
mCAF_markers = all_markers[all_markers$cluster=='mCAF',]$gene
vCAF_markers = all_markers[all_markers$cluster=='vCAF',]$gene
#产生ENTREZID
genelist <- bitr(mCAF_markers, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')

ego1 <- enrichGO(gene = genelist$ENTREZID,
                 OrgDb = org.Hs.eg.db, 
                 ont = "BP",            #此处选择BP
                 pAdjustMethod = "BH",
                 minGSSize = 1,
                 pvalueCutoff =0.05, qvalueCutoff =0.2,
                 readable = TRUE)
ego1_res <- ego1@result
ego1_res$Group <- 'mCAF'
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
pdf('scRNA_Fibem_GO_BP_barplot.pdf', width = 6, height = 5)
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
pdf('scRNA_Fib_GOBP_barplot.pdf', width = 6, height = 6)
ggplot(ego_res, aes(Count, Description)) +
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=LogP),stat='identity') +
  scale_fill_gradient(low="#F39B7FFF",high="#DC0000FF")+
  facet_wrap(~Group, scales="free", ncol=1)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=1))
dev.off()

# 分组气泡图
pdf('scRNA_Fib_GOBP_dotplot.pdf', width = 6.5, height = 4)
ggplot(ego_res, aes(Group, factor(Description, levels=ego_res$Description))) +
  geom_point(aes(color=LogP, size=Count))+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=1))+
  scale_color_gradient(low='#F39B7FFF',high='#DC0000FF')+
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
genes <- unique(rownames(scRNA_Fib@assays$RNA@counts))
aver_dt <- AverageExpression(scRNA_Fib,
                             features = genes,
                             group.by = 'Fib_celltype_stress',
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
group_list <- levels(scRNA_Fib$Fib_celltype_stress)
group_list = factor(group_list,levels = levels(scRNA_Fib$Fib_celltype_stress))

##热图
library(pheatmap) 
annotation_col = data.frame(group=group_list)
rownames(annotation_col) <- colnames(exp_gsva)

pdf('stressFib_celltype_GSVA_interest.pdf', width = 5.5, height = 4.5)
pheatmap::pheatmap(exp_gsva,
                   show_colnames = T, 
                   cluster_rows = F, cluster_cols = F, 
                   color = colorRampPalette(c("#4DBBD5FF", "white", "#E64B35FF"))(50),
                   annotation_col = annotation_col)
dev.off()


####AUCell关键通路评分####
###AUCell与GSVA方法结果可能有不一致！
#成纤维细胞关键geneSet: "IFN"/WP_OXIDATIVE_STRESS_RESPONSE
geneSet_name = c('WP_OXIDATIVE_STRESS_RESPONSE',
                 'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_ENDOGENOUS_ANTIGEN')

cells_AUC <- AUCell_run(scRNA_Fib@assays$RNA@data, my_geneSet)
#提取PATHWAY
geneSet <- "REACTOME_INTERFERON_GAMMA_SIGNALING"   
AUCell_auc <- as.numeric(getAUC(cells_AUC)[geneSet, ])
scRNA_Fib$AUCell <- AUCell_auc #添加至metadata中
head(scRNA_Fib@meta.data)

##绘图可视化
#Seurat自带小提琴图
pdf(file='scRNA_Fib_Stress_vlnplot.pdf',width = 5,height = 4)
VlnPlot(scRNA_Fib,features = 'AUCell', #features也可改为AUCell
        pt.size = 0,group.by = "Fib_celltype_stress",col = cors) #按细胞类型分组
dev.off()

#箱式图
my_comparisons = list(c("Healthy","Melanoma"),c("Healthy","Vitiligo"),c("Melanoma","Vitiligo"))

pdf(file='scRNA_Fib_Stress_boxplot.pdf',width = 7,height = 3)
ggboxplot(scRNA_Fib@meta.data, x="group", y="AUCell", width = 0.6, #按group分组
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
  facet_wrap(~Fib_celltype,ncol=4,scales="free_y") + # 按照细胞类型分面
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  stat_compare_means(comparisons=my_comparisons,
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1),symbols=c("***", "**", "*", "ns")),label="p.signif") # 添加t检验
dev.off()

## tsne图
# 法二  https://zhuanlan.zhihu.com/p/482523999
library(ggrepel)
#提取tsne坐标数据
tsne <- data.frame(scRNA_Fib@meta.data, scRNA_Fib@reductions$tsne@cell.embeddings)
library(ggplot2)
scRNA_Fib$tsne_1 <- tsne$tSNE_1
scRNA_Fib$tsne_2 <- tsne$tSNE_2
mydata <- FetchData(scRNA_Fib,vars = c("tsne_1","tsne_2","AUCell"))
#绘图
pdf(file='scRNA_Fib_Stress_tsne.pdf',width = 4,height = 3)
ggplot(mydata,aes(x = tsne_1,y =tsne_2,colour = AUCell))+
  geom_point(size = 1)+
  scale_color_gradientn(values = seq(0,1,0.2),
                        colours = c('#333366',"#6666FF",'#CC3333','#FFCC33'))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"))
dev.off()


###感兴趣通路基因dotplot
#趋化因子，炎症因子
genes <- unique(my_geneSet$KEGG_CHEMOKINE_SIGNALING_PATHWAY) %>% 
  data.frame() %>% 
  dplyr::filter(grepl('CCL|CXCL', genes$gene)) #CCR|CXCR
genes = c('CCL2','CCL4','CCL5','CCL8','CCL19','CCL27',
          'CXCL1','CXCL2','CXCL3','CXCL12','CXCL13','CXCL14')

pdf(file=paste0('scRNA_Fib_chemokine_dotplot.pdf'),width = 4,height = 3.5)
DotPlot(scRNA_Fib, features = genes, group.by = 'Fib_celltype_stress', assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF')) #颜色
dev.off()


#### Fib相关基因表达比较 ####
geneset = c('TGFB1','FGF2','EGF','CTGF', #生长因子
            'VEGFA', 'IL6','CXCL12', #免疫抑制
            'IFNG','CXCL9','CXCL10','JAK1','JAK2','STAT1')#IFNG相关

DotPlot(scRNA_Fib, features = geneset, group.by = 'Fib_celltype_stress', assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust =1,vjust=1))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF'))

#美化气泡图  9*7
jjDotPlot(object = scRNA_Fib,gene = geneset,
          id = 'Fib_celltype_stress',  #分组设置
          ytree = F,
          dot.col = c('#3C5488FF','white','#E64B35FF'),
          rescale.min = 0,rescale.max = 2,midpoint = 0) #+ coord_flip()



#### 氧化应激评分分群再分析 ####
   ##原因：之前的传统分群找不出白癜风和健康组的差异，因此换个角度分群再分析！
#准备基因集
geneSet = msigdbr(species = "Homo sapiens")  #category = "C2"
geneSet %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
geneSet_name = 'WP_OXIDATIVE_STRESS_RESPONSE'
my_geneSet = subset(geneSet, gs_name %in% geneSet_name)
my_geneSet = my_geneSet %>% split(x = .$gene_symbol, f = .$gs_name)
#AUCell评分
cells_AUC <- AUCell_run(scRNA_Fib@assays$RNA@data, my_geneSet)
geneSet <- "WP_OXIDATIVE_STRESS_RESPONSE"   
AUCell_auc <- as.numeric(getAUC(cells_AUC)[geneSet, ])
scRNA_Fib$AUCell <- AUCell_auc 
##绘图可视化
#Seurat自带小提琴图、箱式图、tsne图同理绘制

##根据氧化应激评分分群(选择分高的Fib1,Fib2,Fib4)
Fib_celltype=c('Stressed Fib','Stressed Fib','Fib','Stressed Fib','Fib',
               'mCAF','vCAF')
Idents(scRNA_Fib) <- scRNA_Fib$Fib_celltype
names(Fib_celltype) <- levels(scRNA_Fib)
scRNA_Fib <- RenameIdents(scRNA_Fib, Fib_celltype)
scRNA_Fib$Fib_celltype_stress <- Idents(scRNA_Fib)
#自定义celltype顺序
cell = c('Stressed Fib','Fib','mCAF','vCAF') 
scRNA_Fib$Fib_celltype_stress <- factor(scRNA_Fib$Fib_celltype_stress,levels = cell)
Idents(scRNA_Fib) <- scRNA_Fib$Fib_celltype_stress
#saveRDS(scRNA_Fib, 'scRNA_Fib.rds')

##再次可视化

##再次绘制细胞比例

##再次亚群marker&富集分析

##提取stressFib组间分析
stressFib <- subset(scRNA_Fib, Fib_celltype_stress == 'Stressed Fib')
 #组间GSVA、组间看IFNG