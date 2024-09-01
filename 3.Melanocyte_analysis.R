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
scRNA_Mela=readRDS('./scRNA_Mela.rds')
table(scRNA_Mela$Mela_celltype,scRNA_Mela$group)


#### 提取黑素细胞亚群 ####
scRNA_Mela=subset(scRNA,celltype=='Melanocyte')

# 提细胞亚群,重新标准化、降维聚类
scRNA_Mela<- FindVariableFeatures(scRNA_Mela, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNA_Mela)
scRNA_Mela <- ScaleData(scRNA_Mela, features = scale.genes)
scRNA_Mela<- RunPCA(scRNA_Mela, features = VariableFeatures(scRNA_Mela))
DimPlot(scRNA_Mela, reduction = "pca", group.by = "orig.ident")
ElbowPlot(scRNA_Mela)

## 重新harmony
library(harmony)
set.seed(1000)
scRNA_Mela <- RunHarmony(scRNA_Mela, group.by.vars = "orig.ident")
DimPlot(scRNA_Mela, reduction = "harmony", group.by = "orig.ident")
ElbowPlot(scRNA_Mela,reduction = 'harmony')

scRNA_Mela <- FindNeighbors(scRNA_Mela, reduction = 'harmony',dims = 1:10)
scRNA_Mela <- FindClusters(scRNA_Mela, resolution = 0.8)
scRNA_Mela <- RunUMAP(scRNA_Mela, reduction = 'harmony',dims = 1:10)
scRNA_Mela <- RunTSNE(scRNA_Mela, reduction = 'harmony',dims = 1:10)
DimPlot(scRNA_Mela,reduction = "tsne",label = T,split.by = 'group')  #此处tsne比umap好看
dev.off()

##如果文献中提供的marker效果不好，则直接按各cluster的top marker编号
all_markers <- FindAllMarkers(object = scRNA_Mela,
                              only.pos = T, #only.pos改为T则只输出高表达gene
                              min.pct = 0.25,logfc.threshold = 0.5)
# 提取各亚群Top5Marker
top5_markers <- all_markers %>% group_by(cluster) %>% top_n(5, avg_log2FC) #输出每个组的前5个高的gene
# 先用编号注释过渡，后面根据Top5Marker生物学意义注释
Mela_celltype=c('Mela1','Mela2','Mela3','Mela4','Mela5','Mela6') 
Idents(scRNA_Mela) <- scRNA_Mela@meta.data$seurat_clusters
names(Mela_celltype) <- levels(scRNA_Mela)
scRNA_Mela <- RenameIdents(scRNA_Mela, Mela_celltype)
scRNA_Mela@meta.data$Mela_celltype0 <- Idents(scRNA_Mela)


###功能注释教程：https://mp.weixin.qq.com/s/-9V0WQG5v1L9Fjy9d7PFKg
### 法一：根据功能基因集进行注释
geneSet = msigdbr(species = "Homo sapiens")  #category = "C2"
geneSet %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
#挑选感兴趣的基因集
geneSet_name = c('BIOCARTA_MELANOCYTE_PATHWAY', #melanocytic
                 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION')   #mesenchymal-like
geneSet = subset(geneSet, gs_name %in% geneSet_name)
geneSet = geneSet %>% split(x = .$gene_symbol, f = .$gs_name)
cells_AUC <- AUCell_run(scRNA_Mela@assays$RNA@data, geneSet) #计算基因集评分
#提取PATHWAY
my_geneSet <- "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"   
AUCell_auc <- as.numeric(getAUC(cells_AUC)[my_geneSet, ])
scRNA_Mela$AUCell <- AUCell_auc #添加至metadata中
head(scRNA_Mela@meta.data)
#Seurat自带小提琴图
VlnPlot(scRNA_Mela,features = 'AUCell', 
        pt.size = 0,group.by = "Mela_celltype0",col = cors) #按细胞类型分组
VlnPlot(scRNA_Mela,features = 'SERPINE1', #TYR,MITF代表色素沉着;SERPINE2,TIMP1代表EMT
        pt.size = 0,group.by = "Mela_celltype0",col = cors)
### 法二：根据标记基因进行注释
Idents(scRNA_Mela)=scRNA_Mela$seurat_clusters
Mela_marker <- c(unique(geneSet$BIOCARTA_MELANOCYTE_PATHWAY), #Melanocytic
                 unique(geneSet$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION))  #Mesenchymal-like 
DotPlot(scRNA_Mela, 
        features = Mela_marker,
        group.by = "seurat_clusters") + coord_flip()
#注释
Mela_celltype=c('Mesenchymal-like','Melanocytic','Melanocytic','Melanocytic','Melanocytic')
Idents(scRNA_Mela) <- scRNA_Mela@meta.data$seurat_clusters
names(Mela_celltype) <- levels(scRNA_Mela)
scRNA_Mela <- RenameIdents(scRNA_Mela, Mela_celltype)
scRNA_Mela@meta.data$Mela_celltype <- Idents(scRNA_Mela)
#自定义celltype顺序
cell = c('Melanocytic','Mesenchymal-like') 
scRNA_Mela$Mela_celltype <- factor(scRNA_Mela$Mela_celltype,levels = cell)
Idents(scRNA_Mela) <- scRNA_Mela$Mela_celltype
#保存
saveRDS(scRNA_Mela, 'scRNA_Mela.rds')

### 可视化
## umap图
DimPlot(scRNA_Mela,reduction = "tsne",group.by = 'Mela_celltype0',label = T,cols = cors) 
DimPlot(scRNA_Mela,reduction = "tsne",group.by = 'group',label = F,cols = cors) 
library(scRNAtoolVis)
pdf(file=paste0('scRNA_Mela_celltype_umap.pdf'),width = 6,height = 4)
clusterCornerAxes(object = scRNA_Mela, reduction = 'tsne',
                  clusterCol = "Mela_celltype0", #分类依据
                  #noSplit = F, groupFacet = 'group', #分组时使用
                  cellLabel = T, cellLabelSize = 3.5)+  #标签
  #addCircle = TRUE, cicAlpha = 0.1, nbin = 200) + #画圈
  scale_color_npg() + scale_fill_npg()
dev.off()

##marker气泡图
pdf(file=paste0('scRNA_Mela_celltype_dotplot.pdf'),width = 5,height = 6)
DotPlot(scRNA_Mela, features = top5_markers$gene, group.by='Mela_celltype0',assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF')) #颜色
dev.off()




####细胞比例####
library(reshape2)
library(ggplot2)
#define the color
library(ggsci)
cors <- pal_npg()(10) #定义颜色
#准备绘图数据
prop_df <- table(scRNA_Mela$Mela_celltype0,scRNA_Mela$group) %>% melt()
colnames(prop_df) <- c("Cluster","Sample","Number")
prop_df$Cluster <- factor(prop_df$Cluster)
prop_df$Proportion <- ave(prop_df$Number, prop_df$Sample, 
                          FUN = function(x) x/sum(x))

#比例图1  尺寸6*3
pdf(file=paste0('scRNA_Mela_proportion_barplot.pdf'),width = 4.5,height = 2)
ggplot(data = prop_df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=cors) + #配色
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
dev.off()


####黑素细胞infercnv####
#服务器上完成


####黑素细胞组间DEG####
## 1)FindMarkers组间DEG####
Idents(scRNA_Mela)='group'
all_markers <- FindAllMarkers(object = scRNA_Mela,
                              only.pos = T, #only.pos改为T则只输出高表达gene
                              min.pct = 0.25,logfc.threshold = 0.25)
top10_markers <- all_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top10_markers <- top10_markers[!duplicated(top10_markers$gene),]
# 各亚群平均表达量提取
genes <- unique(top10_markers$gene)
aver_dt <- AverageExpression(scRNA_Mela,
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

pdf(file="scRNA_Mela_top10marker_heatmap.pdf", width=4, height=6)
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


####2)GO/KEGG####
V_markers=read.csv('DEG_V vs H/Melanocyte.csv') %>% filter(p_val<0.05 & avg_log2FC>0.5)
M_markers=read.csv('DEG_M vs H/Melanocyte.csv') %>% filter(p_val<0.05 & avg_log2FC>0.5)

genelist <- bitr(M_markers$X, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')

##GO
ego1 <- enrichGO(gene = genelist$ENTREZID,
                 OrgDb = org.Hs.eg.db, 
                 ont = "BP",            #此处选择BP
                 pAdjustMethod = "BH",
                 minGSSize = 1,
                 pvalueCutoff =0.05, qvalueCutoff =0.2,
                 readable = TRUE)
ego1_res <- ego1@result
ego1_res$Group <- 'V vs H'
ego1_res$LogP <- -log(ego1_res$p.adjust) #计算-logP

ego2 <- enrichGO(gene = genelist$ENTREZID,
                 OrgDb = org.Hs.eg.db, 
                 ont = "BP",            #此处选择BP
                 pAdjustMethod = "BH",
                 minGSSize = 1,
                 pvalueCutoff =0.05, qvalueCutoff =0.2,
                 readable = TRUE)
ego2_res <- ego2@result
ego2_res$Group <- 'M vs H'
ego2_res$LogP <- -log(ego2_res$p.adjust) #计算-logP

# 合并富集结果
ego_res = rbind(ego1_res[1:10,],ego2_res[1:10,])
# 分组柱状图
pdf('scRNA_Mela_KEGG_barplot.pdf', width = 6, height = 6)
ggplot(ego_res, aes(Count, Description)) +
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=LogP),stat='identity') +
  scale_fill_gradient(low="#F39B7FFF",high="#DC0000FF")+
  facet_wrap(~Group, scales="free", ncol=1)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=1))
dev.off()
# 分组气泡图
pdf('scRNA_Mela_GOBP_dotplot.pdf', width = 6.5, height = 4)
ggplot(ego_res, aes(Group, Description)) +
  geom_point(aes(color=LogP, size=Count))+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=0.5))+
  scale_color_gradient(low='#4DBBD5FF',high='#DC0000FF')+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1))
dev.off()


####3)GSEA/GSVA####
##准备分析用的geneSet
geneSet = msigdbr(species = "Homo sapiens") #category = "C2", subcategory = 'KEGG'
geneSet %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
geneSet = geneSet %>% dplyr::select(gs_name, gene_symbol)
# 获取GeneList
geneList <- all_markers$avg_log2FC            
names(geneList) <- all_markers$gene         # 对GeneList命名
geneList <- sort(geneList, decreasing = T)   # 从高到低排序
# 开始GSEA富集分析
GSEA_enrichment <- GSEA(geneList,                 # 排序后的gene
                        TERM2GENE = geneSet,      # 基因集
                        pvalueCutoff = 0.1,      # P值阈值
                        pAdjustMethod = "BH")     # 校正P值的计算方法
result <- data.frame(GSEA_enrichment)
# plot
library(enrichplot)
geneset_plot <- c("KEGG_OXIDATIVE_PHOSPHORYLATION")
gseaplot2(GSEA_enrichment,
          geneSetID = geneset_plot,
          color = cors,
          title = 'Specific GO_BP',
          rel_heights = c(1.3, 0.3, 0.6),
          pvalue_table = F)

###GSVA
library(GSVA)
# 各亚群平均表达量提取
genes <- unique(rownames(scRNA_Mela@assays$RNA@counts))
aver_dt <- AverageExpression(scRNA_Mela,
                             features = genes,
                             group.by = 'group',
                             slot = 'data')
aver_dt <- as.data.frame(aver_dt$RNA)
aver_dt <- aver_dt[rowSums(aver_dt)>0,]  #过滤细胞表达量全为零的基因
# 余下操作参考后面GSVA

####3)取趋势相反DEG交集####
#变化趋势相反的DEG
DEG_Viti = read.csv("CD8T_DEG_H_vs_V.csv")
DEG_Mel = read.csv("CD8T_DEG_H_vs_M.csv")

DEG_Vitiligo_up = subset(DEG_Viti, Type == 'up')
DEG_Vitiligo_down = subset(DEG_Viti, Type == 'down')
DEG_Melanoma_up = subset(DEG_Mel, Type == 'up')
DEG_Melanoma_down = subset(DEG_Mel, Type == 'down')

inverse_DEG = union(intersect(DEG_Vitiligo_up$Gene_Symbol, DEG_Melanoma_down$Gene_Symbol), 
                    intersect(DEG_Vitiligo_down$Gene_Symbol, DEG_Melanoma_up$Gene_Symbol))
write.csv(inverse_DEG, 'T_inverse_DEG.csv',row.names = F, col.names = F)

##韦恩图
library(ggvenn)
venn_dat <- list(
  "DEG_Vitiligo_up" = DEG_Vitiligo_up$Gene_Symbol,
  "DEG_Melanoma_down" = DEG_Melanoma_down$Gene_Symbol,
  "DEG_Vitiligo_down" = DEG_Vitiligo_down$Gene_Symbol,
  "DEG_Melanoma_up" = DEG_Melanoma_up$Gene_Symbol)

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
pdf(file=paste0('CD8T_inverseDEG_dotplot.pdf'),width = 5,height = 4.5)
DotPlot(scRNA_Mela, features = inverse_DEG, group.by = 'group',assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF')) #颜色
dev.off()

####4)趋势相同/反基因散点图####
##以黑素细胞为例，提取两组差异基因
Mela_Viti = read.csv('Melanocye.csv')
Mela_Mel = read.csv('Melanocye.csv')
#合并数据
ids = intersect(rownames(Mela_Viti), rownames(Mela_Mel))
df = data.frame(Mela_Viti = Mela_Viti[ids,'avg_log2FC'],
                Mela_Mel = Mela_Mel[ids,'avg_log2FC'])
rownames(df) = ids
df$change = ifelse(df$Mela_Viti>0.5 & df$Mela_Mel< -0.5, 'inverse',
                   ifelse(df$Mela_Viti< -0.5 & df$Mela_Mel>0.5, 'inverse','not'))
#添加基因标签
df$Label = ""   #新加一列label
df_inverse = subset(df, change == 'inverse')
df_inverse <- df_inverse[order(df_inverse$Mela_Viti-df_inverse$Mela_Mel), ]  #对差异基因的p值进行从小到大的排序
df_inverse$Gene <- rownames(df_inverse)
up.genes <- tail(df_inverse$Gene[which(df_inverse$change == "inverse")], 5) #高表达的基因中选择fdr值最小的5个
down.genes <- head(df_inverse$Gene[which(df_inverse$change == "inverse")], 5) #低表达的基因中选择fdr值最小的5个
#将up.genes和down.genes合并，加入到Label中
df.top5.genes <- c(as.character(up.genes), as.character(down.genes))
df$Label[match(df.top5.genes, df$Gene)] <- df.top5.genes

##对角散点图
ggscatter(df, x = "Mela_Viti", y = "Mela_Mel",
          color = "change",
          palette = c("#E41A1C","black","#E41A1C"),
          size = 1.0,
          label = df$Label,
          font.label = 10,
          repel = T,
          xlab = "V vs H",
          ylab = "M vs H") +
  geom_hline(yintercept = c(0), linetype = "dashed") +
  geom_vline(xintercept = c(0), linetype = "dashed") #改logfc阈值



####亚群marker富集分析####
####0)提取亚群marker####
all_markers <- FindAllMarkers(object = scRNA_Mela,
                              only.pos = T, #only.pos改为T则只输出高表达gene
                              min.pct = 0.25,logfc.threshold = 0.25)
top5_markers <- all_markers %>% group_by(cluster) %>% top_n(5, avg_log2FC) #输出每个组的前5个高的gene
top10_markers <- all_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
# 各亚群平均表达量提取
genes <- unique(top10_markers$gene)
aver_dt <- AverageExpression(scRNA_Mela,
                             features = genes,
                             group.by = 'Mela_celltype',
                             slot = 'data')
aver_dt <- as.data.frame(aver_dt$RNA)

##热图
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

pdf(file="CD8T_markertop10_heatmap.pdf", height=6, width=5)
pheatmap(as.matrix(aver_dt),
         scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col = cell_anno,
         annotation_row = gene_anno,
         annotation_colors = anno_col, #注释配色
         color = mycol, #热图配色
         border_color = 'white',
         gaps_row = c(10,20)) #描边颜色
dev.off()

#产生ENTREZID（富集分析用）
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)

Tem_markers = all_markers[all_markers$cluster=='CD8+ Tem',]$gene
Tcm_markers = all_markers[all_markers$cluster=='CD8+ Tcm',]$gene
Trm_markers = all_markers[all_markers$cluster=='CD8+ Trm',]$gene

genelist <- bitr(Tem_markers, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')


####1)GO分析####
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
  scale_fill_gradient(low="#F39B7FFF",high="#DC0000FF")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=0.5))
dev.off()
# 分组气泡图
ego_res = rbind(ego1_res[1:10,],ego2_res[1:10,],ego3_res[1:10,])

pdf('CD8T_GO_BP_dotplot.pdf', width = 7, height = 5)
ggplot(ego_res, aes(Group, Description)) +
  geom_point(aes(color=LogP, size=Count))+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=0.5))+
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

# 分组气泡图
kk_res = rbind(kk1_res[1:10,],kk2_res[1:10,],kk3_res[1:10,])

pdf('CD8T_KEGG_dotplot.pdf', width = 6, height = 5)
ggplot(kk_res, aes(Group, Description)) +
  geom_point(aes(color=LogP, size=Count))+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust=0.5,vjust=0.5))+
  scale_color_gradient(low='#4DBBD5FF',high='#DC0000FF')+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1))
dev.off()


####3)GSVA分析####
library(GSVA)
# 各亚群平均表达量提取
genes <- unique(rownames(scRNA_Mela@assays$RNA@counts))
aver_dt <- AverageExpression(scRNA_Mela,
                             features = genes,
                             group.by = 'Mela_celltype0',
                             slot = 'data')
aver_dt <- as.data.frame(aver_dt$RNA)
aver_dt <- aver_dt[rowSums(aver_dt)>0,]  #过滤细胞表达量全为零的基因

##法一：选择msigdb的KEGG基因集
#msigdbr包教程：https://www.jianshu.com/p/f2febb3123d8
geneSet = msigdbr(species = "Homo sapiens")  #category = "C2"
geneSet %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
#挑选感兴趣的基因集  死亡通路
geneSet_name = c('KEGG_APOPTOSIS',
                 'REACTOME_CELLULAR_SENESCENCE',
                 'GOBP_RESPONSE_TO_OXIDATIVE_STRESS',
                 'GOBP_AUTOPHAGY_OF_MITOCHONDRION',
                 'GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS',
                 'HP_ABNORMAL_ACTIVITY_OF_MITOCHONDRIAL_RESPIRATORY_CHAIN',
                 'WP_FERROPTOSIS',
                 'REACTOME_REGULATED_NECROSIS',
                 'HP_AVASCULAR_NECROSIS',
                 'GOBP_NECROPTOTIC_SIGNALING_PATHWAY',
                 'REACTOME_PYROPTOSIS',
                 'BIOCARTA_P53_PATHWAY',
                 'GOBP_POSITIVE_REGULATION_OF_DNA_DAMAGE_RESPONSE_SIGNAL_TRANSDUCTION_BY_P53_CLASS_MEDIATOR_RESULTING_IN_TRANSCRIPTION_OF_P21_CLASS_MEDIATOR'
                 )
geneSet = subset(geneSet, gs_name %in% geneSet_name)
geneSet = geneSet %>% split(x = .$gene_symbol, f = .$gs_name)#基因集是list
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
exp_gsva <- gsva(as.matrix(aver_dt), geneSet, min.sz > 1) #表达矩阵需转换为matrix格式
#method=c("gsva", "ssgsea", "zscore", "plage")
#exp_gsva <- MinMax(scale(exp_gsva), -2, 2) #标准化
#调整exp_gsva行顺序
rownames(exp_gsva)
exp_gsva <- exp_gsva[c(5,3,2,6,1,9,8,4,7,11,10,12),] 
rownames(exp_gsva) <- c('Oxidative stress',
                        'Endoplasmic reticulum stress',
                        'Autophagy of mitochondrion',
                        'Mitochondrion dysfunction',
                        'p53 pathway',
                        'Cellular senescence',
                        'Apoptosis',
                        'Necroptosis',
                        'Avascular necrosis',
                        'Regulated necrosis',
                        'Pyroptosis',
                        'Ferroptosis')

#设置参考水平
group_list <- c('Mela1','Mela2','Mela3','Mela4','Mela5','Mela6')
group_list = factor(group_list,levels = c('Mela1','Mela2','Mela3','Mela4','Mela5','Mela6'))
#选择相应表达水平中位绝对偏差中位数排名前20位的通路（看情况运行）
mad_scores <- apply(exp_gsva, 1, mad)
top_genes <- order(mad_scores, decreasing = TRUE)[1:20]
exp_gsva_top <- exp_gsva[top_genes, ]

##热图
library(pheatmap) 
annotation_col = data.frame(group=group_list)
rownames(annotation_col) <- colnames(exp_gsva)

pdf('scRNA_Mela_celltype_GSVA_death.pdf', width = 5, height = 3.5)
pheatmap::pheatmap(exp_gsva,
                   scale = "row",
                   show_colnames = T, 
                   cluster_rows = F, cluster_cols = F, 
                   color = colorRampPalette(c("#3C5488FF", "white", "#DC0000FF"))(50),
                   annotation_col = annotation_col)
dev.off()



####AUCell关键通路评分####
###AUCell与GSVA方法结果可能有不一致！
##通过R获取msigdbr数据集
#geneSet和GSVA一样
geneSet_name = c('WP_EMBRYONIC_STEM_CELL_PLURIPOTENCY_PATHWAYS',
                 'BHATTACHARYA_EMBRYONIC_STEM_CELL')
cells_AUC <- AUCell_run(scRNA_Mela@assays$RNA@data, my_geneSet)
#提取PATHWAY
geneSet <- 'REACTOME_PYROPTOSIS'  
AUCell_auc <- as.numeric(getAUC(cells_AUC)[geneSet, ])
scRNA_Mela$AUCell <- AUCell_auc #添加至metadata中
head(scRNA_Mela@meta.data)

##绘图可视化
#Seurat自带小提琴图
pdf(file='Mela_REACTOME_PYROPTOSIS_vlnplot.pdf',width = 5,height = 4)
VlnPlot(scRNA_Mela, features = 'AUCell', 
        group.by = "group", pt.size = 0, cols = cors) + 
  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(title = "REACTOME_PYROPTOSIS", x="") +  #title命名
  theme(legend.position="right") +
  stat_compare_means(method = "wilcox.test", label.x.npc ="center", size = 5) + 
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black") #添加均值
dev.off()

#箱式图
my_comparisons = list(c("Healthy","Melanoma"),c("Healthy","Vitiligo"),c("Melanoma","Vitiligo"))

pdf(file='CD8T_CYTOTOXICITY_boxplot.pdf',width = 4,height = 3)
ggboxplot(scRNA_Mela@meta.data, x="group", y="AUCell", width = 0.6, #按group分组
          color = "black",#轮廓颜色
          fill="group",#填充
          palette = cors,
          xlab = F, #不显示x轴的标签
          bxp.errorbar=T,#显示误差条
          bxp.errorbar.width=0.5, #误差条大小
          size=0.5, #箱型图边线的粗细
          outlier.shape=NA, #不显示outlier
          legend = "right") + 
  ylim(0,0.1) +
  facet_wrap(~Mela_celltype) + # 按照细胞类型分面
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test") # 添加t检验
dev.off()

## tsne图
# 法二  https://zhuanlan.zhihu.com/p/482523999
library(ggrepel)
#提取tsne坐标数据
tsne <- data.frame(scRNA_Mela@meta.data, scRNA_Mela@reductions$tsne@cell.embeddings)
library(ggplot2)
scRNA_Mela$tsne_1 <- tsne$tSNE_1
scRNA_Mela$tsne_2 <- tsne$tSNE_2
mydata <- FetchData(scRNA_Mela,vars = c("tsne_1","tsne_2","AUCell"))
#绘图
pdf(file='scRNA_Mela_stemness_tsne.pdf',width = 4,height = 3)
ggplot(mydata,aes(x = tsne_1,y =tsne_2,colour = AUCell))+
  geom_point(size = 1)+
  scale_color_gradientn(values = seq(0,1,0.2),
                        colours = c('#333366',"#6666FF",'#CC3333','#FFCC33'))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"))
dev.off()



###感兴趣通路基因dotplot
geneSet = msigdbr(species = "Homo sapiens")  #category = "C2"
geneSet %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
#挑选感兴趣的基因集  死亡通路
geneSet_name = c('REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES')
geneSet = subset(geneSet, gs_name %in% geneSet_name)
geneSet = geneSet %>% split(x = .$gene_symbol, f = .$gs_name)#基因集是list
genes = unique(geneSet$REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES) %>% data.frame()
colnames(genes) = 'gene'
genes = genes %>% dplyr::filter(grepl('CCR|CXCR', genes$gene)) #CCR|CXCR
genes = genes %>% dplyr::filter(grepl('CCL|CXCL', genes$gene)) #CCR|CXCR
#Melanocyte
genes = c('IFNAR1','IFNAR2','IFNGR1','IFNGR2',
         'JAK1','JAK2','STAT1','STAT2','STAT3')

pdf(file=paste0('Mela_chemokine_dotplot.pdf'),width = 5,height = 4)
DotPlot(scRNA_Mela, features = genes, group.by = 'Mela_celltype0',assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF')) #颜色
dev.off()



#### 氧化应激评分分群再分析 ####
##原因：之前的传统分群找不出白癜风和健康组的差异，因此换个角度分群再分析！
#准备基因集
geneSet = msigdbr(species = "Homo sapiens")  #category = "C2"
geneSet %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
geneSet_name = 'WP_OXIDATIVE_STRESS_RESPONSE'
my_geneSet = subset(geneSet, gs_name %in% geneSet_name)
my_geneSet = my_geneSet %>% split(x = .$gene_symbol, f = .$gs_name)
#AUCell评分
cells_AUC <- AUCell_run(scRNA_Mela@assays$RNA@data, my_geneSet)
geneSet <- "WP_OXIDATIVE_STRESS_RESPONSE"   
AUCell_auc <- as.numeric(getAUC(cells_AUC)[geneSet, ])
scRNA_Mela$AUCell <- AUCell_auc 
##绘图可视化
#Seurat自带小提琴图、箱式图、tsne图同理绘制

##根据氧化应激评分分群(选择分高的Fib1,Fib2,Fib4)
Mela_celltype=c('Malignant Mel','Mel','Malignant Mel','Malignant Mel','Stressed Mel','Mel')
Idents(scRNA_Mela) <- scRNA_Mela$Mela_celltype0
names(Mela_celltype) <- levels(scRNA_Mela)
scRNA_Mela <- RenameIdents(scRNA_Mela, Mela_celltype)
scRNA_Mela$Mela_celltype1 <- Idents(scRNA_Mela)
#自定义celltype顺序
cell = c('Malignant Mel','Stressed Mel','Mel') 
scRNA_Mela$Mela_celltype1 <- factor(scRNA_Mela$Mela_celltype1,levels = cell)
Idents(scRNA_Mela) <- scRNA_Mela$Mela_celltype1
saveRDS(scRNA_Mela, 'scRNA_Mela.rds')

##再次可视化

##再次绘制细胞比例

##再次亚群marker&富集分析
Mf_stress_markers = FindMarkers(scRNA_Mela,ident.1 = 'Stressed Mφ', 
                                verbose = FALSE, test.use = 'wilcox',
                                min.pct = 0.25, logfc.threshold = 0)
