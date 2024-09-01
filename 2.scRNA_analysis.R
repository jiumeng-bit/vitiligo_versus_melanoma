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

####细胞比例####
library(reshape2)
library(ggplot2)
#define the color
library(ggsci)
cors <- pal_npg()(10) #定义颜色
#准备绘图数据
prop_df <- table(scRNA$celltype,scRNA$orig.ident) %>% melt()
colnames(prop_df) <- c("Cluster","Sample","Number")
prop_df$Cluster <- factor(prop_df$Cluster)
prop_df$Proportion <- ave(prop_df$Number, prop_df$Sample, 
                          FUN = function(x) x/sum(x))

#比例图1  尺寸6*3
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
#比例图2  尺寸8*4
ggplot(prop_df, aes(x = Cluster, y = Proportion, fill = Sample)) +
  geom_col(position=position_dodge(0.8), width = 0.8)+
  scale_fill_manual(values=c('#3C5488FF','#F39B7FFF','#4DBBD5FF','#E64B35FF')) + #配色
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45))

##批量散点箱式图 法一
#计算各组中各样本不同细胞群比例
Cellratio <- prop.table(table(Idents(scRNA), scRNA$orig.ident), margin = 2)
Cellratio <- data.frame(Cellratio)
colnames(Cellratio) <- c('celltype','Sample','Freq')
Cellratio$group <- ifelse(grepl('H',Cellratio$Sample),'Healthy',
                          ifelse(grepl('M',Cellratio$Sample),'Melanoma','Vitiligo'))

my_comparisons <- list(c("Healthy","Melanoma"),c("Healthy","Vitiligo"),c("Melanoma","Vitiligo"))

pdf(file="1.2.scRNA_proportion_boxplot.pdf", width=8, height=3.5)
ggboxplot(Cellratio, x="group", y="Freq",  #按group分组
          fill="group",color = "black",  #填充/轮廓颜色
          width = 0.6,
          palette = cors,
          xlab = F, #不显示x轴的标签
          bxp.errorbar=T,#显示误差条
          bxp.errorbar.width=0.5, #误差条大小
          size=0.5, #箱型图边线的粗细
          outlier.shape=NA, #不显示outlier
          legend = "right") + 
  facet_wrap(~celltype, ncol=4, scales = "free_y") + # 按照细胞类型分面，ncol为列数
  geom_jitter(size=1) + #添加散点
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.border = element_rect(colour = "black", fill=NA)) #+
  #stat_compare_means(comparisons = my_comparisons, method = "t.test") # 添加t检验
dev.off()

##批量散点箱式图 法二
#计算各组中各样本不同细胞群比例
Cellratio <- prop.table(table(Idents(scRNA), scRNA$orig.ident), margin = 2)
Cellratio <- data.frame(Cellratio)
library(reshape2)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#长数据转为宽数据
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]
###添加分组信息
sample <- unique(scRNA$orig.ident)
group <- c(rep('Healthy',5),rep('Vitiligo',10),rep('Melanoma',3))
samples <- data.frame(sample, group)#创建数据框
rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']#R添加列
cellper$group <- samples[rownames(cellper),'group']#R添加列

###作图展示
pplist = list()
sce_groups = c('Keratinocyte',"Melanocyte","Fibroblast","Endothelial cell", "Smooth muscle", 
               "T & NK","Langerhans cell","Mononulear phagocyte") 
for(group_ in sce_groups){
  cellper_  = cellper %>% dplyr::select(one_of(c('sample','group',group_)))#选择一组数据
  colnames(cellper_) = c('sample','group','percent')#对选择数据列命名
  cellper_$percent = as.numeric(cellper_$percent)#数值型数据
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                      lower = quantile(percent, 0.25),
                                                      mean = mean(percent),
                                                      median = median(percent))#上下分位数
  print(group_)
  print(cellper_$median)
  
  pp1 = ggplot(cellper_,aes(x=group,y=percent)) + #ggplot作图
    geom_jitter(shape = 21,aes(fill=group),width = 0.25) + 
    stat_summary(fun=mean, geom="point", color="grey60") +
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Percentage') +
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1) + scale_color_npg() + scale_fill_npg()
  ###组间t检验分析
  labely = max(cellper_$percent)
  compare_means(percent ~ group,  data = cellper_)
  my_comparisons <- list(c("Healthy","Melanoma"),c("Healthy","Vitiligo"),c("Melanoma","Vitiligo"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons, size = 3, method = "t.test")
  pplist[[group_]] = pp1
}
plot_grid(pplist[["T cell"]],  #自己选择想展示的celltype
          pplist[["NK cell"]],
          pplist[["Langerhans cell"]],
          pplist[["Mononulear phagocyte"]])



####细胞类群差异分析####
## 计算细胞marker
 #教程：https://www.jianshu.com/p/f5c8f9ea84af
scRNA_Vitiligo <- subset(scRNA, group == 'Vitiligo')
scRNA_Melanoma <- subset(scRNA, group == 'Melanoma')
all_markers <- FindAllMarkers(object = scRNA_Melanoma,
                              only.pos = F, #only.pos改为T则只输出高表达gene
                              min.pct = 0.25,logfc.threshold = 0.5)
write.csv(all_markers,file = "2.1.all_markers_Melanoma.csv")
# 提取各亚群Top5Marker
top5_markers <- all_markers %>% group_by(cluster) %>% top_n(5, avg_log2FC) #输出每个组的前5个高的gene
top5_markers <- top5_markers[!duplicated(top5_markers$gene),]
# 各亚群平均表达量提取
genes <- unique(top5_markers$gene)
aver_dt <- AverageExpression(scRNA,
                             features = genes,
                             group.by = 'celltype',
                             slot = 'data')
aver_dt <- as.data.frame(aver_dt$RNA)


### 绘制多组火山图
#方法一
source('D:/R语言/【单细胞分析】/【最全教程(闲鱼)】/scRNA_code/scVolcano_functions.R')
scVolcano(all_markers)  #10*8
#方法二
library(scRNAtoolVis)
pdf(file="2.2.scRNA_marker_multivolcano.pdf", width=10, height=8)
jjVolcano(diffData = all_markers,
          log2FC.cutoff = 0.5, topGeneN = 5,
          aesCol = c('#4DBBD5FF','#E64B35FF'))
dev.off()
#极坐标图  #8*8
jjVolcano(diffData = all_markers,
          log2FC.cutoff = 0.5, topGeneN = 5, sisze  = 3,
          aesCol = c('#4DBBD5FF','#E64B35FF'),
          tile.col = corrplot::COL2('RdBu', 15)[4:12],
          fontface = 'italic', polar = T)

### 绘制热图
library(ComplexHeatmap)
library(cols4all)
#热图配色自定义
mycol <- colorRampPalette(c("#4DBBD5FF", "white", "#E64B35FF"))(50)
#行列注释配色自定义
celltype_col <- c4a('10', 8) #此处改为细胞类型数目

gene_anno <- data.frame(gene_anno = top5_markers$cluster,#行注释：Top5marker对应celltype
                        row.names = top5_markers$gene)
cell_anno <- data.frame(cell_anno = colnames(aver_dt),#列注释：celltype
                        row.names = colnames(aver_dt))
names(celltype_col) <- cell_anno$cell_anno
anno_col <- list(cell_anno = celltype_col,gene_anno = celltype_col)

pdf(file="2.3.scRNA_top5marker_heatmap.pdf", width=6, height=7)
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

##更便捷方法
#法一：平均热图
library(scRNAtoolVis)
AverageHeatmap(object = scRNA,
               markerGene = top5_markers$gene,
               clusterAnnoName = F,
               htCol = c("#4DBBD5FF", "white", "#DC0000FF"), #热图方块颜色
               #showRowNames = F,markGenes = annoGene, #标记展示基因，annoGene是基因集合
               annoCol = TRUE, myanCol = cors[1:8])
#法二：原始热图
pdf(file=paste0('scRNA_KRT_DEG_heatmap.pdf'),width = 5,height = 4)
DoHeatmap(scRNA_KRT, top10_markers$gene, 
          group.by = 'KRT_celltype', group.colors = cors) +
  scale_fill_gradientn(colors = c("#3C5488FF", "white", "#DC0000FF"))
dev.off()


####组间单细胞差异分析####
## 1)pseudobulks差异分析####
bs = split(colnames(scRNA),scRNA$orig.ident)
ct = do.call(
  cbind,lapply(names(bs), function(x){ 
    # x=names(bs)[[1]]
    kp =colnames(scRNA) %in% bs[[x]]
    rowSums( as.matrix(scRNA@assays$RNA@counts[, kp]  ))
  })
)
colnames(ct) <- names(bs)
phe = unique(scRNA@meta.data[,c('orig.ident','group')])#样本&信息，自行修改
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
  mutate(Type = if_else(padj > 0.05, "ns",
                        if_else(abs(log2FoldChange) < 0.5, "ns",
                                if_else(log2FoldChange >= 0.5, "up", "down")))) %>%
  arrange(desc(abs(log2FoldChange))) %>% rownames_to_column("Gene_Symbol")
table(DEG_DESeq2$Type)
write.csv(DEG_DESeq2,"DEG_H vs V.csv")

## 2)FindMarkers差异分析####
names(scRNA@meta.data)
unique(scRNA$group)

scRNA$celltype.group <- paste(scRNA$celltype, scRNA$group, sep = "_")
Idents(scRNA) <- "celltype.group"

##选择某细胞亚群差异分析（以T cell为例）
  #ident.1是要去比的组(Disease)，ident.2是被比较的组(Control)
CELLDEG <- FindMarkers(scRNA,ident.1 = 'T cell_Healthy',ident.2 = 'T cell_Vitiligo', 
                    verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
##批量操作
cellfordeg <- levels(scRNA$celltype)

for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(scRNA, ident.1 = paste0(cellfordeg[i],'_Vitiligo'), ident.2 = paste0(cellfordeg[i],'_Healthy'), 
                         verbose = FALSE)
  write.csv(CELLDEG,paste0(cellfordeg[i],".csv"))
}

##差异基因可视化
library(dplyr)
top10 <- CELLDEG  %>% top_n(n = 10, wt = avg_log2FC) %>% row.names()
Top10 <- as.data.frame(top10)
#scRNA <- ScaleData(scRNA, features =  rownames(scRNA))
DoHeatmap(scRNA, features = top10, size = 3)


##小提琴图（以T细胞为例）
Idents(scRNA) <- "celltype"
VlnPlot(scRNA,features = top10,split.by = 'group',idents = 'T & NK')
##umap图
FeaturePlot(scRNA,features = top10,split.by = 'group',idents = 'T & NK')
##气泡图
DotPlot(scRNA,features = top10,split.by ='group',idents = 'T & NK',cols = cors)



####趋势相同/反基因####
##以黑素细胞为例，提取两组差异基因
Mela_Viti = read.csv('Melanocye.csv')
Mela_Mel = read.csv('Melanocye.csv')
#合并数据
ids = intersect(rownames(Mela_Viti), rownames(Mela_Mel))
df = data.frame(Mela_Viti = Mela_Viti[ids,'avg_log2FC'],
                Mela_Mel = Mela_Mel[ids,'avg_log2FC'])
rownames(df) = ids
df$change = ifelse(df$Mela_Viti>0.5 & df$Mela_Mel>0.5, 'up',
                   ifelse(df$Mela_Viti< -0.5 & df$Mela_Mel< -0.5, 'down','not'))
#添加基因标签
df$Label = ""   #新加一列label
df <- df[order(df$Mela_Viti+df$Mela_Mel), ]  #对差异基因的p值进行从小到大的排序
df$Gene <- rownames(df)
up.genes <- tail(df$Gene[which(df$change == "up")], 5) #高表达的基因中选择fdr值最小的5个
down.genes <- head(df$Gene[which(df$change == "down")], 5) #低表达的基因中选择fdr值最小的5个
#将up.genes和down.genes合并，加入到Label中
df.top5.genes <- c(as.character(up.genes), as.character(down.genes))
df$Label[match(df.top5.genes, df$Gene)] <- df.top5.genes

##对角散点图
ggscatter(df, x = "Mela_Viti", y = "Mela_Mel",
          color = "change",
          palette = c("#377EB8","black","#E41A1C"),
          size = 1.5,
          label = df$Label,
          font.label = 10,
          repel = T,
          xlab = "V vs H",
          ylab = "M vs H") +
  geom_hline(yintercept = c(0), linetype = "dashed") +
  geom_vline(xintercept = c(0), linetype = "dashed") #改logfc阈值


####组间DEG富集分析####
###1)GO####
genelist <- bitr(DEG_DESeq2$Gene_Symbol, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')

ego1 <- enrichGO(gene = genelist$ENTREZID,
                 OrgDb = org.Hs.eg.db, 
                 ont = "BP",            #此处选择BP
                 pAdjustMethod = "BH",
                 minGSSize = 1,
                 pvalueCutoff =0.05, qvalueCutoff =0.2,
                 readable = TRUE)
ego1_res <- ego1@result
ego1_res$Group <- 'M vs H'
ego1_res$LogP <- -log(ego1_res$p.adjust) #计算-logP
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


###2)KEGG####
kk1 <- enrichKEGG(gene = genelist$ENTREZID,
                  keyType = 'kegg',
                  organism = 'hsa',
                  pvalueCutoff = 0.1,
                  qvalueCutoff =0.1)
kk1_res <- kk1@result
kk1_res$Group <- 'CD8+ Tem'
kk1_res$LogP <- -log(kk1_res$p.adjust)
##kk2_res,kk3_res同理计算


###3)GSEA####
#教程：https://zhuanlan.zhihu.com/p/667669372
#msigdbr包教程：https://www.jianshu.com/p/f2febb3123d8
library(clusterProfiler)

# 使用logFC进行基因排序
head(DEG_DESeq2)
# 加载基因集，基因集介绍往下滑
geneSet = msigdbr(species = "Homo sapiens",category = "C2") #category = "C2", subcategory = 'KEGG'
geneSet %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
geneSet = geneSet %>% dplyr::select(gs_name, gene_symbol)
# 接下来我们进行基因排序
geneList <- DEG_DESeq2$log2FoldChange               # 获取GeneList
names(geneList) <- DEG_DESeq2$Gene_Symbol      # 对GeneList命名
geneList <- sort(geneList, decreasing = T)  # 从高到低排序

# 开始GSEA富集分析
GSEA_enrichment <- GSEA(geneList,                 # 排序后的gene
                        TERM2GENE = geneSet,      # 基因集
                        pvalueCutoff = 1,      # P值阈值
                        minGSSize = 20,           # 最小基因数量
                        maxGSSize = 1000,         # 最大基因数量
                        eps = 0,                  # P值边界
                        pAdjustMethod = "BH")     # 校正P值的计算方法
result <- data.frame(GSEA_enrichment)
#只提取与xxx有关的基因集   HYPOXIA|WNT
result1 = result  %>%
  dplyr::filter(grepl('HYPO', Description))
# 特定通路绘图
library(enrichplot) # 富集结果可视化
#gseaplot2(GSEA_enrichment, c(2), color = "red3", pvalue_table = T)
pdf('mela_GSEA_plot_C7.pdf', width = 6, height = 5)
gseaplot2(GSEA_enrichment,'REACTOME_INTERFERON_SIGNALING', 
          color = "red3", pvalue_table = T)
dev.off()


# 将通路分为激活和抑制两个部分进行展示
library(ggplot2)     # 画图图
pdf('mela_GSEA_dotplot_C7.pdf', width = 12, height = 10)
dotplot(GSEA_enrichment, showCategory = 10, split = ".sign") + facet_grid(~.sign) +
  theme(plot.title = element_text(size = 10, color = "black", hjust = 0.5),
        axis.title = element_text(size = 10,color = "black"), 
        axis.text = element_text(size = 10,color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 1 ),
        legend.position = "right",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))
dev.off()


####关键通路富集分析####
##1)GSVA####
library(GSVA)
# 各分组各亚群平均表达量提取
scRNA$celltype.group <- paste(scRNA$celltype, scRNA$group, sep = "_")
Idents(scRNA) <- "celltype.group"

genes <- unique(rownames(scRNA@assays$RNA@counts))
aver_dt <- AverageExpression(scRNA,
                             features = genes,
                             group.by = 'celltype',
                             slot = 'data')
aver_dt <- as.data.frame(aver_dt$RNA)
aver_dt <- aver_dt[rowSums(aver_dt)>0,]  #过滤细胞表达量全为零的基因


##选择msigdb的KEGG基因集
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
geneSet = subset(geneSet, gs_name %in% geneSet_name)
geneSet = geneSet %>% split(x = .$gene_symbol, f = .$gs_name)#基因集是list

##GSVA分析 
exp_gsva <- gsva(as.matrix(aver_dt),  #表达矩阵需转换为matrix格式
                 gset.idx.list = geneSet) #method=c("gsva","ssgsea","zscore","plage")

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
  

##热图
library(pheatmap) 
#添加列注释
annotation_col = data.frame(group = rep(c('Healthy','Melanoma','Vitiligo'),8),
                            celltype = rep(sort(levels(scRNA$celltype)),each = 3))
rownames(annotation_col) <- colnames(exp_gsva)

pdf('scRNA_GSVA_interest.pdf', width = 10, height = 4)
pheatmap::pheatmap(exp_gsva,
                   show_colnames = F, 
                   cluster_rows = F, cluster_cols = F, 
                   #gaps_row = c(3,5,14),  #添加行分隔
                   gaps_col = c(3,6,9,12,15,18,21),  #添加列分隔
                   color = colorRampPalette(c("#3C5488FF", "white", "#DC0000FF"))(50),
                   annotation_col = annotation_col)
dev.off()



##2)AUCell####
#geneSet和GSVA一样

cells_AUC <- AUCell_run(scRNA@assays$RNA@data, geneSet)
#提取PATHWAY
my_geneSet <- "KEGG_CHEMOKINE_SIGNALING_PATHWAY"   
AUCell_auc <- as.numeric(getAUC(cells_AUC)[my_geneSet, ])
scRNA$AUCell <- AUCell_auc #添加至metadata中
head(scRNA@meta.data)

##绘图可视化
#Seurat自带小提琴图
pdf(file='scRNA_chemokine_vlnplot.pdf',width = 5,height = 4)
VlnPlot(scRNA,features = 'AUCell', #features也可改为AUCell
        pt.size = 0,group.by = "celltype",col = cors) #按细胞类型分组
dev.off()

#箱式图
my_comparisons = list(c("Healthy","Melanoma"),c("Healthy","Vitiligo"),c("Melanoma","Vitiligo"))

pdf(file='scRNA_chemokine_boxplot.pdf',width = 8,height = 4)
ggboxplot(scRNA@meta.data, x="group", y="AUCell", width = 0.6, #按group分组
          color = "black",#轮廓颜色
          fill="group",#填充
          palette = cors,
          xlab = F, #不显示x轴的标签
          bxp.errorbar=T,#显示误差条
          bxp.errorbar.width=0.5, #误差条大小
          size=0.5, #箱型图边线的粗细
          outlier.shape=NA, #不显示outlier
          legend = "right") + 
  #ylim(0,0.15) + 
  facet_wrap(~celltype,ncol=4,scales="free_y") + # 按照细胞类型分面
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test") # 添加t检验
dev.off()

## umap图
# 法二  https://zhuanlan.zhihu.com/p/482523999
library(ggrepel)
#提取umap坐标数据
umap <- data.frame(scRNA@meta.data, scRNA@reductions$umap@cell.embeddings)
library(ggplot2)
scRNA$UMAP_1 <- umap$UMAP_1
scRNA$UMAP_2 <- umap$UMAP_2
mydata <- FetchData(scRNA,vars = c("UMAP_1","UMAP_2","AUCell"))
#绘图
pdf(file='CD8T_T_ANTIGEN_umap.pdf',width = 5,height = 4)
ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = AUCell))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('#333366',"#6666FF",'#CC3333','#FFCC33')) + 
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"))
dev.off()


###感兴趣通路基因dotplot
# 各亚群平均表达量提取
genes = unique(geneSet$KEGG_CHEMOKINE_SIGNALING_PATHWAY) %>% data.frame()
colnames(genes) = 'gene'
genes = genes %>% dplyr::filter(grepl('CCL|CXCL', genes$gene)) #CCR|CXCR


pdf(file=paste0('CD8T_NOTCH_dotplot.pdf'),width = 5,height = 4.5)
DotPlot(scRNA, features = genes$gene, group.by = 'celltype', assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust=1,vjust=1))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF'))
dev.off()
