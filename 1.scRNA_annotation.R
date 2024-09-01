rm(list=ls())
#以后seurat包可根据安装路径指定版本
.libPaths(c("D:/R语言/Rlibrary/win-library/4.2","~/seurat_v5","C:/Program Files/R/R-4.3.1/library"))
library(Seurat)
packageVersion("Seurat")

library(stringr)
library(dplyr)
library(future)
library(future.apply)
library(dplyr)
library(msigdbr)
library(clusterProfiler)
library(devtools)
library(harmony)
library(clustree)
library(ggplot2)
library(reshape2)

options(future.globals.maxSize = 60000 * 1024^2)
getwd()

#define the color
library(ggsci)
cors <- pal_npg()(10) #定义颜色

#### 读取单细胞数据 ####
### 读取白癜风样本数据
## 循环读取
setwd("E:/【科研学习】/【皮肤科】/白癜风课题/MR+单细胞")

filename <- paste('processed_data_new/',list.files('processed_data_new/'),sep = '')
sceList <- lapply(filename, function(x){
  obj <- CreateSeuratObject(counts = Read10X(x),
                            min.cells = 10, 
                            min.features = 200,
                            project = str_split(x,'/')[[1]][2])
})
names(sceList) <- list.files('processed_data_new/')
scVitiligo <- merge(sceList[[1]],sceList[-1],add.cell.ids = names(sceList),project='OC')
dim(scVitiligo) # 18924个gene， 52845个cell
scVitiligo@meta.data[1:5,] #简单看一下
table(scVitiligo$orig.ident)
# 添加分组信息
scVitiligo@meta.data$group <- ifelse(substr(rownames(scVitiligo@meta.data),1,1)=='H',
                                       'Healthy','Vitiligo') 
scVitiligo@meta.data[1:5,]
table(scVitiligo$group)
# 保存数据集
saveRDS(scVitiligo,file = 'scRNA.rds')


### 读取黑色素瘤样本数据
scMelanoma1 = Read10X_h5('Melanoma_GSE215120/GSM6622299_CM1_filtered_feature_bc_matrix.h5')
scMelanoma1 = CreateSeuratObject(counts = scMelanoma1, project = "M01", 
                                 min.cells=10, min.features=200, names.delim = "_")
scMelanoma2 = Read10X_h5('Melanoma_GSE215120/GSM6622300_CM2_filtered_feature_bc_matrix.h5')
scMelanoma2 = CreateSeuratObject(counts = scMelanoma2, project = "M02", 
                                 min.cells=10, min.features=200, names.delim = "_")
scMelanoma3 = Read10X_h5('Melanoma_GSE215120/GSM6622301_CM3_filtered_feature_bc_matrix.h5')
scMelanoma3 = CreateSeuratObject(counts = scMelanoma3, project = "M03", 
                                 min.cells=10, min.features=200, names.delim = "_")
#合并数据集
scMelanoma <- merge(scMelanoma1,y=c(scMelanoma2,scMelanoma3))
rm(scMelanoma1,scMelanoma2,scMelanoma3,hd5,mat)
gc()
#sample_name <- str_split_fixed(rownames(scMelanoma@meta.data),'_',n=2) %>% data.frame()
#scMelanoma$orig.ident <- paste0('M0',sample_name$X2)
scMelanoma$group <- 'Melanoma'
table(scMelanoma$orig.ident)


### 数据合并
## 合并单细胞数据集
scRNA <- merge(scVitiligo,y=scMelanoma)
table(scRNA$group)
rm(scMelanoma,scVitiligo)
gc()

#### 质量控制 ####
##计算质控指标
#计算细胞中线粒体核糖体基因比例
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
hist(scRNA[["percent.mt"]]$percent.mt) #展示线粒体基因百分比的直方图

col.num <- length(levels(scRNA@active.ident))
# 绘制基因特征的小提琴图
VlnPlot(scRNA,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        cols =rainbow(col.num), 
        pt.size = 0, #不需要显示点，可以设置pt.size = 0
        ncol = 3) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
# 绘制测序深度的相关性图
plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

## 数据质控
# 设置质控标准
print(c("请输入允许基因数和核糖体比例，示例如下：", "minGene=500", "maxGene=4000", "pctMT=20"))
minCount=500
minGene=200
maxGene=5000
pctMT=20
scRNA <- subset(scRNA, subset = nFeature_RNA> minCount & 
                  nFeature_RNA > minGene & nFeature_RNA < maxGene & 
                  percent.mt < pctMT)
dim(scRNA)
col.num <- length(levels(scRNA@active.ident))
VlnPlot(scRNA,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        group.by = 'orig.ident', 
        pt.size = 0, 
        ncol = 3) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 

#增加细胞周期评分
s.genes <- cc.genes$s.genes #S期
g2m.genes <- cc.genes$g2m.genes # G2期
scRNA <- CellCycleScoring(scRNA, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


#### 标准化&降维聚类 ####
library(Seurat)
library(tidyverse)
library(patchwork)

# 标准化
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
#高变基因
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNA)
#归一化
scRNA <- ScaleData(scRNA, features = scale.genes)

# PCA降维
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))
DimPlot(scRNA, reduction = "pca", group.by = "orig.ident")
VizDimLoadings(scRNA, dims = 1:2, reduction = "pca")
ElbowPlot(scRNA,ndims = 50)

## harmony去批次
  # 批次效应怎么看：不同样本/分组细胞分布区别很大(没有融合到一起)
library(harmony)
scRNA <- RunHarmony(scRNA, group.by.vars = "orig.ident")
DimPlot(scRNA, reduction = "harmony", group.by = "orig.ident")
ElbowPlot(scRNA, reduction = 'harmony')

# 降维聚类
scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:10) #指定harmony,PCA数选择10
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, reduction = 'harmony', dims = 1:10)
#scRNA <- RunTSNE(scRNA, reduction = 'harmony', dims = 1:10)

# 去批次后绘制聚类图
DimPlot(scRNA,reduction = "umap",group.by = "seurat_clusters",label = T)
DimPlot(scRNA,reduction = "umap",group.by = "group",label = T)
DimPlot(scRNA,reduction = "umap",group.by = "orig.ident",label = T)
DimPlot(scRNA,reduction = "umap",label = T,split.by = 'group')

saveRDS(scRNA,file = 'scRNA_merge_new.rds')


#### 细胞注释 ####
#### 1)自动注释 ####
library(SingleR)
library(celldex)
sc_data<-as.matrix(scRNA@assays$RNA@counts)
hpca.se <- HumanPrimaryCellAtlasData()
hpca.se$label.main
pred.hesc <- SingleR(test = sc_data, ref = hpca.se, assay.type.test=1,
                     labels = hpca.se$label.main)
pred.hesc$labels###预测的结果
#添加单细胞注释标签
scRNA@meta.data$singleR_label <- unlist(lapply(scRNA@meta.data$seurat_clusters, function(x){pred.hesc$labels[x]}))
#绘图
DimPlot(scRNA,reduction = 'umap',group.by = 'singleR_label',label = T)
DimPlot(scRNA,reduction = 'tsne',group.by = 'singleR_label',label = T)

#### 2)手动注释 ####
#根据每种细胞marker确定该细胞对应哪些聚类
Idents(scRNA)=scRNA$seurat_clusters
marker <- c("KRT1","KRT10","KRT14","KRT15", #角质形成细胞 Keratinocyte
            "DCT","PMEL","TYRP1","MLANA", #黑素细胞 Melanocyte
            "COL1A1","DCN","SFRP2","TWIST2", #成纤维细胞 Fibroblast
            "PECAM1","CLEC14A","AQP1","ECSCR.1", #内皮细胞 Endothelial
            "TAGLN","ACTA2","MYL9","NR2F2", #平滑肌细胞 Smooth muscle
            "TRAC","CD3D","TRBC2","CD3E", #T细胞 T cell 
            "NKG7","GNLY",                #NK细胞 NK cell
            "CD207","CD1A","FCGBP","S100B", #朗格汉斯细胞 Langerhans(皮肤组织的树突状细胞)
            "LYZ","CD1C","IL1B","CLEC10A") #单核巨噬细胞 Mono phagocyte
DotPlot(scRNA, 
        features = marker,
        group.by = "seurat_clusters") + coord_flip() +
  scale_color_gradientn(colours = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF')) #颜色

## 要求准确
celltype=c('T & NK',
           'Melanocyte','T & NK','Keratinocyte','Keratinocyte','Keratinocyte',
           'Langerhans','Mono phagocyte','Fibroblast','Keratinocyte','Mono phagocyte',
           'Unknown','Smooth muscle','Unknown','Melanocyte','Unknown',
           'Endothelial','Unknown','Melanocyte','Keratinocyte','Keratinocyte',
           'Melanocyte','Melanocyte','Unknown','Langerhans','Smooth muscle',
           'Melanocyte','Unknown','Fibroblast')

# 重命名Idents
names(celltype) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, celltype)
# 命名celltype
scRNA@meta.data$celltype <- Idents(scRNA)
#自定义celltype顺序
cell = c('Keratinocyte',"Melanocyte","Fibroblast","Endothelial", "Smooth muscle", 
         "T & NK","Langerhans","Mono phagocyte") 
scRNA$celltype <- factor(scRNA$celltype,levels = cell)
Idents(scRNA) <- scRNA$celltype
# 删除unknown
scRNA <- subset(scRNA, celltype != 'Unknown')
scRNA@meta.data$celltype <- Idents(scRNA) #完全去除Unknown
table(scRNA@meta.data$celltype)
scRNA@meta.data[1:5,]


####3)GPT注释####
library(GPTCelltype)
library(openai)
Sys.setenv(OPENAI_API_KEY = '')
suppressWarnings({
  all.markers <- FindAllMarkers(scRNA)
})
res <- gptcelltype(all.markers, 
                   tissuename = 'human PBMC', 
                   model = 'gpt-4')
cat(res)
#将cat(res)放到ChatGPT4上按照提问模板询问

## tsne/umap图
#按照celltype/group/sample映射
DimPlot(scRNA,reduction="umap",label=T,group.by="celltype",cols=cors) #7*5
DimPlot(scRNA,reduction="umap",label=T,group.by="group",cols=cors) #7*5
DimPlot(scRNA,reduction="umap",label=F,group.by="orig.ident") #7*5
#按照group分开呈现
scRNA$group <- factor(scRNA$group,
                      levels = c('Healthy','Melanoma','Vitiligo'))
DimPlot(scRNA,reduction="umap",label=F,group.by="celltype",split.by="group",cols=cors)

#提取白癜风和黑素瘤各自scRNA数据
#scViti <- subset(scRNA, group=='Vitiligo')
#scMela <- subset(scRNA, group=='Melanoma')
#保存scRNA数据
saveRDS(scRNA,file = "scRNA_merge_annotate_new.rds")
#saveRDS(scViti,file = "scRNA_Vitiligo_annotate.rds")
#saveRDS(scMela,file = "scRNA_Melanoma_annotate.rds")
scRNA <- readRDS("scRNA_merge_annotate_new.rds")



####注释后可视化####
 #美化教程：https://mp.weixin.qq.com/s/1C64jsl08oTVoxbNQ49ioA

##绘制nature级别umap图   #7*5
library(scRNAtoolVis)
pdf(file=paste0('1.1.scRNA_celltype_umap.pdf'),width = 7,height = 5)
clusterCornerAxes(object = scRNA, reduction = 'umap',
                  clusterCol = "Cell_type", #分类依据
                  #noSplit = F, groupFacet = 'group', #分组时使用
                  cellLabel = T, cellLabelSize = 3.5,  #标签
                  addCircle = TRUE, cicAlpha = 0.1, nbin = 200) + #画圈
  scale_color_npg() + scale_fill_npg()
dev.off()

colnames(scRNA@meta.data)[12] <- 'Cell_type'
##绘制细胞marker umap图
pdf(file=paste0('2.2.scRNA_celltype_markerumap.pdf'),width = 12,height = 6)
FeatureCornerAxes(object = scRNA, reduction = 'umap',
                  groupFacet = NULL,
                  relLength = 0.5,
                  relDist = 0.2,
                  aspect.ratio = 1,
                  features = c("KRT1","DCT","COL1A1","PECAM1",
                               "TAGLN","TRAC","CD207","LYZ"))
dev.off()

##绘制各细胞marker小提琴图
#确定选用展示的markers（根据经验选择）
features = marker
length(features)
#绘图
pdf(file=paste0('2.1.scRNA_celltype_vlnplot.pdf'),width = 5,height = 8)
VlnPlot(scRNA, features = features, #6*6
        stack = TRUE, 
        sort = TRUE, 
        cols = cors,
        split.by =  "celltype" , #每种cluster 一个颜色
        flip = TRUE,
        split.plot = TRUE) +
  theme(legend.position = "none")
dev.off()

#获得个性化绘图数据
library(reshape2)
vln.df=scRNA[["RNA"]]@data[top10$gene[1:7],] %>% as.matrix %>% t() %>% data.frame()
vln.df$barcode=rownames(vln.df)
vln.df=melt(vln.df,id="barcode")
colnames(vln.df)[c(2,3)]=c("gene","exp")
anno=scRNA@meta.data
anno$barcode=rownames(anno)
vln.df=merge(vln.df,anno,by="barcode")  #metadata与基因表达数据合并
vln.df$gene=factor(vln.df$gene,levels = top10$gene[1:7]) #控制画图的基因顺序
vln.df%>%ggplot(aes(celltype,exp))+geom_violin(aes(fill=gene),scale = "width")+
  facet_grid(vln.df$gene~.,scales = "free_y")+
  scale_fill_brewer(palette = "Set3",direction = 1)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none")

##绘制细胞marker气泡图
pdf(file=paste0('2.2.scRNA_celltype_dotplot.pdf'),width = 6,height = 7)
DotPlot(scRNA, features = marker, assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust = 1))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF')) #颜色
dev.off()

##ggplot个性化绘制
p <- DotPlot(scRNA, features = marker, group.by = "celltype", assay = "RNA")
data <- p$data
data%>%ggplot(aes(x=id,y=features.plot,size=pct.exp,color=avg.exp.scaled))+
  geom_point()+
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))+
  scale_x_discrete("")+scale_y_discrete("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank()
  )


#### 提取单核吞噬细胞 ####
 #https://zhuanlan.zhihu.com/p/375318689 朗格汉斯细胞只是DC的一个亚群！
scRNA_Mono=subset(scRNA,celltype=='Mononulear phagocyte') 
# 提细胞亚群,重新降维聚类
scRNA_Mono <- FindVariableFeatures(scRNA_Mono, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNA_Mono)
scRNA_Mono <- ScaleData(scRNA_Mono, features = scale.genes)
scRNA_Mono<- RunPCA(scRNA_Mono, features = VariableFeatures(scRNA_Mono))
DimPlot(scRNA_Mono, reduction = "pca", group.by = "orig.ident")
ElbowPlot(scRNA_Mono)

## 重新harmony
library(harmony)
set.seed(1000)
scRNA_Mono <- RunHarmony(scRNA_Mono, group.by.vars = "orig.ident")
DimPlot(scRNA_Mono, reduction = "harmony", group.by = "orig.ident")
ElbowPlot(scRNA_Mono,reduction = 'harmony')

scRNA_Mono <- FindNeighbors(scRNA_Mono, reduction = 'harmony',dims = 1:10)
sce_res <- scRNA_Mono
for (i in c(0.01, 0.05, 0.1, 0.15, 0.2, 0.3,0.4, 0.5,0.8,1)){
  sce_res <- FindClusters(sce_res,resolution = i)
}
scRNA_Mono <- FindClusters(scRNA_Mono, resolution = 0.01)
scRNA_Mono <- RunUMAP(scRNA_Mono, reduction = 'harmony',dims = 1:10)
scRNA_Mono <- RunTSNE(scRNA_Mono, reduction = 'harmony',dims = 1:10)
DimPlot(scRNA_Mono,reduction = "umap",label = T) 

## 单核吞噬细胞亚群注释参考
# marker：https://www.jianshu.com/p/16280b2da028,https://zhuanlan.zhihu.com/p/520281423
# DC亚群marker：https://mp.weixin.qq.com/s/-w0e9IhNcYyWUMqajJ9xqQ
# 背景知识：https://zhuanlan.zhihu.com/p/653640358
Idents(scRNA_Mono)=scRNA_Mono$seurat_clusters
# 最新注释
Mono_marker <- c('CXCL5','CCL23','CCL3','CCL20','CCL4', # Mφ_CCL23
                 'IL1B','IL1RN','FCN1','DUSP6', # Mφ_FCN1
                 'SPP1','MMP9','FPR3', # Mφ_SPP1
                 "C1QA","C1QB",'TREM2','APOE','APOC1','GPNMB',  # Mφ_APOE
                 'FCGR3A','CX3CR1','TCF7L2','LRRC25','FCGR3B', # Mono_FCGR3A
                 'CD14','S100A12','MND1') # Mono_CD14
DotPlot(scRNA_Mono, 
        features = Mono_marker,
        group.by = "seurat_clusters") + coord_flip()

##如果文献中提供的marker效果不好，则直接按各cluster的top marker编号
all_markers <- FindAllMarkers(object = scRNA_Mono,
                              only.pos = T, #only.pos改为T则只输出高表达gene
                              min.pct = 0.25,logfc.threshold = 0.5)
# 提取各亚群Top5Marker
top10_markers <- all_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC) #输出每个组的前5个高的gene
# 根据Top5Marker生物学意义注释
Mono_celltype=c('Mφ_CXCL8','Mφ_IRF8','Mφ_FSCN1')

Idents(scRNA_Mono) <- scRNA_Mono@meta.data$seurat_clusters
names(Mono_celltype) <- levels(scRNA_Mono)
scRNA_Mono <- RenameIdents(scRNA_Mono, Mono_celltype)
scRNA_Mono$Mono_celltype <- Idents(scRNA_Mono)
#自定义celltype顺序
cell = c('Mφ_CXCL8','Mφ_IRF8','Mφ_FSCN1') 
scRNA_Mono$Mono_celltype <- factor(scRNA_Mono$Mono_celltype,levels = cell)
Idents(scRNA_Mono) <- scRNA_Mono$Mono_celltype
saveRDS(scRNA_Mono, 'scRNA_Mono.rds')

## umap图  #6*4
library(scRNAtoolVis)
pdf(file=paste0('Mono_celltype_umap.pdf'),width = 7,height = 5)
clusterCornerAxes(object = scRNA_Mono, reduction = 'umap',
                  clusterCol = "Mono_celltype", #分类依据
                  #noSplit = F, groupFacet = 'group', #分组时使用
                  cellLabel = T, cellLabelSize = 3.5)+  #标签
  #addCircle = TRUE, cicAlpha = 0.1, nbin = 200) + #画圈
  scale_color_npg() + scale_fill_npg()
dev.off()

##绘制细胞marker气泡图
pdf(file=paste0('Mono_celltype_dotplot.pdf'),width = 5,height = 7)
DotPlot(scRNA_Mono, features = top10_markers$gene, assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF')) #颜色
dev.off()


## 绘制比例图
prop_df <- table(scRNA_Mono$Mono_celltype,scRNA_Mono$group) %>% melt()
colnames(prop_df) <- c("Cluster","Sample","Number")
prop_df$Cluster <- factor(prop_df$Cluster)
prop_df$Proportion <- ave(prop_df$Number, prop_df$Sample, 
                          FUN = function(x) x/sum(x))

#比例图1  尺寸6*3
p1 = ggplot(data = prop_df, aes(x =Number, y = Sample, fill =  Cluster)) +
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
ggsave(p1, file = 'Mono_proportion_barplot.pdf', width = 6, height = 3) #保存图片



#### 提取树突状细胞 ####
scRNA_DC=subset(scRNA,celltype=='Langerhans cell')
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
scRNA_DC <- FindClusters(scRNA_DC)
scRNA_DC <- RunUMAP(scRNA_DC, reduction = 'harmony',dims = 1:10)
scRNA_DC <- RunTSNE(scRNA_DC, reduction = 'harmony',dims = 1:10)
DimPlot(scRNA_DC,reduction = "umap",label = T) 
DimPlot(scRNA_DC,reduction = "umap",label = T,split.by = 'group') 
dev.off()

Idents(scRNA_DC)=scRNA_DC$seurat_clusters
# 最新注释
DC_marker <- c("GZMB","TSPAN13","LILRA4","ITM2C","IRF4", # pDC
               "CLEC9A","CPNE3", "WDFY4","XCR1", # cDC1
               'CLEC10A',"CD1C",'FCER1A','CD1E', # cDC2
               "LAMP3", "CCR7","FSCN1","IL7R","IDO1") # cDC3

DotPlot(scRNA_DC, 
        features = DC_marker,
        group.by = "seurat_clusters") + coord_flip()

Mono_celltype=c('cDC2',
                'cDC2','cDC2','cDC2','cDC2','cDC2',
                'cDC1','cDC2','pDC','pDC','cDC3',
                'cDC3','pDC')

Idents(scRNA_DC) <- scRNA_DC@meta.data$seurat_clusters
names(Mono_celltype) <- levels(scRNA_DC)
scRNA_DC <- RenameIdents(scRNA_DC, Mono_celltype)
scRNA_DC$Mono_celltype <- Idents(scRNA_DC)
# 删除unknown
scRNA_DC <- subset(scRNA_DC, Mono_celltype != 'Unknown')
scRNA_DC$Mono_celltype <- Idents(scRNA_DC) #完全去除Unknown
table(scRNA_DC$Mono_celltype)
#自定义celltype顺序
cell = c('pDC','cDC1',"cDC2","cDC3") 
scRNA_DC$Mono_celltype <- factor(scRNA_DC$Mono_celltype,levels = cell)
Idents(scRNA_DC) <- scRNA_DC$Mono_celltype
saveRDS(scRNA_DC, 'scRNA_DC_new.rds')


## umap图  #6*4
library(scRNAtoolVis)
pdf(file=paste0('DC_celltype_umap.pdf'),width = 7,height = 5)
clusterCornerAxes(object = scRNA_DC, reduction = 'umap',
                  clusterCol = "Mono_celltype", #分类依据
                  #noSplit = F, groupFacet = 'group', #分组时使用
                  cellLabel = T, cellLabelSize = 3.5)+  #标签
                  #addCircle = TRUE, cicAlpha = 0.1, nbin = 200) + #画圈
  scale_color_npg() + scale_fill_npg()
dev.off()

##绘制细胞marker气泡图
pdf(file=paste0('DC_celltype_dotplot.pdf'),width = 5,height = 4)
DotPlot(scRNA_DC, features = DC_marker, assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF')) #颜色
dev.off()


## 绘制比例图
prop_df <- table(scRNA_DC$Mono_celltype,scRNA_DC$group) %>% melt()
colnames(prop_df) <- c("Cluster","Sample","Number")
prop_df$Cluster <- factor(prop_df$Cluster)
prop_df$Proportion <- ave(prop_df$Number, prop_df$Sample, 
                          FUN = function(x) x/sum(x))

#比例图1  尺寸6*3
p1 = ggplot(data = prop_df, aes(x =Number, y = Sample, fill =  Cluster)) +
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
ggsave(p1, file = 'DC_proportion_barplot.pdf', width = 6, height = 3) #保存图片



#### 提取成纤维细胞亚群 ####
scRNA_Fibro=subset(scRNA,celltype=='Fibroblast')

# 提细胞亚群,重新降维聚类
scRNA_Fibro<- FindVariableFeatures(scRNA_Fibro, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNA_Fibro)
scRNA_Fibro <- ScaleData(scRNA_Fibro, features = scale.genes)
scRNA_Fibro<- RunPCA(scRNA_Fibro, features = VariableFeatures(scRNA_Fibro))
DimPlot(scRNA_Fibro, reduction = "pca", group.by = "orig.ident")
ElbowPlot(scRNA_Fibro)

## 重新harmony
library(harmony)
set.seed(1000)
scRNA_Fibro <- RunHarmony(scRNA_Fibro, group.by.vars = "orig.ident")
DimPlot(scRNA_Fibro, reduction = "harmony", group.by = "orig.ident")
ElbowPlot(scRNA_Fibro,reduction = 'harmony')

scRNA_Fibro <- FindNeighbors(scRNA_Fibro, reduction = 'harmony',dims = 1:10)
scRNA_Fibro <- FindClusters(scRNA_Fibro, resolution = 0.05)
scRNA_Fibro <- RunUMAP(scRNA_Fibro, reduction = 'harmony',dims = 1:10)
scRNA_Fibro <- RunTSNE(scRNA_Fibro, reduction = 'harmony',dims = 1:10)
DimPlot(scRNA_Fibro,reduction = "umap",label = T) 
dev.off()

# 成纤维细胞亚群注释参考
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7181753/
Idents(scRNA_Fibro)=scRNA_Fibro$seurat_clusters
Fibro_marker <- c("WISP2","SLPI", "CTHRC1",'MFAP5','TSPAN8', #Secretory-reticular
                  "APCDD1", "ID1", "WIF1",'COL18A1','PTGDS', #Secretory-papillary
                  "CCL19", "APOE","CXCL2",'CXCL3','EFEMP1', #Pro-inflammatory
                  'ASPN',"POSTN",'GPC3','TNN','SFRP1') #Mesenchymal 
DotPlot(scRNA_Fibro, 
        features = Fibro_marker,
        group.by = "seurat_clusters") + coord_flip()
Fibro_celltype=c('spFib',
                 'piFib','srFib','spFib','mFib','piFib',
                 'piFib','srFib','Unknown','Unknown','mFib',
                 'Unknown','spFib')

##如果文献中提供的marker效果不好，则直接按各cluster的top marker编号
all_markers <- FindAllMarkers(object = scRNA_Fibro,
                              only.pos = T, #only.pos改为T则只输出高表达gene
                              min.pct = 0.25,logfc.threshold = 0.5)
# 提取各亚群Top5Marker
top10_markers <- all_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC) #输出每个组的前5个高的gene
# 根据Top5Marker生物学意义注释
Fibro_celltype=c('irFib','pbFib','mrFib','edFib','irFib')


Idents(scRNA_Fibro) <- scRNA_Fibro@meta.data$seurat_clusters
names(Fibro_celltype) <- levels(scRNA_Fibro)
scRNA_Fibro <- RenameIdents(scRNA_Fibro, Fibro_celltype)
scRNA_Fibro@meta.data$Fibro_celltype <- Idents(scRNA_Fibro)
# 删除edFib（干扰簇）
scRNA_Fibro = subset(scRNA_Fibro, Fibro_celltype != 'edFib')
scRNA_Fibro$Fibro_celltype <- Idents(scRNA_Fibro)  #完全去除
table(scRNA_Fibro$Fibro_celltype)
#自定义celltype顺序
cell = c('piFib','srFib','spFib','mFib') 
cell = c('irFib','pbFib','mrFib') 
scRNA_Fibro$Fibro_celltype <- factor(scRNA_Fibro$Fibro_celltype,levels = cell)
Idents(scRNA_Fibro) <- scRNA_Fibro$Fibro_celltype
table(scRNA_Fibro$group,scRNA_Fibro$Fibro_celltype)

saveRDS(scRNA_Fibro, 'scRNA_Fibro_new2.rds')

## umap图
DimPlot(scRNA_Fibro,reduction = "umap",label = T,cols = cors) 

library(scRNAtoolVis)
pdf(file=paste0('Fibro_celltype_umap.pdf'),width = 7,height = 5)
clusterCornerAxes(object = scRNA_Fibro, reduction = 'umap',
                  clusterCol = "Fibro_celltype", #分类依据
                  #noSplit = F, groupFacet = 'group', #分组时使用
                  cellLabel = T, cellLabelSize = 3.5)+  #标签
  #addCircle = TRUE, cicAlpha = 0.1, nbin = 200) + #画圈
  scale_color_npg() + scale_fill_npg()
dev.off()

##绘制细胞marker气泡图
pdf(file=paste0('Fibro_celltype_dotplot.pdf'),width = 5,height = 7)
DotPlot(scRNA_Fibro, features = genes, assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 60, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF')) #颜色
dev.off()


## 绘制比例图
prop_df <- table(scRNA_Fibro$Fibro_celltype,scRNA_Fibro$group) %>% melt()
colnames(prop_df) <- c("Cluster","Sample","Number")
prop_df$Cluster <- factor(prop_df$Cluster)
prop_df$Proportion <- ave(prop_df$Number, prop_df$Sample, 
                          FUN = function(x) x/sum(x))

#比例图1  尺寸6*3
pdf(file=paste0('Firo_proportion_barplot.pdf'),width = 6,height = 3)
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


