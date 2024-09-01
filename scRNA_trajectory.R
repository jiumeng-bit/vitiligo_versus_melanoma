##拟时序分析目的：进一步细分细胞亚群，鉴定细胞的不同状态
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
library(monocle)
library(ggsci)

cors <- pal_npg()(10) #定义颜色

####一、T细胞####
##流程：找ordering gene，拟时序分析，得到各branch差异基因，分别功能富集分析

####准备文件####
scRNAsub = readRDS('scRNA_T&NK_new.rds')
#根据实际需求选择要分析的细胞亚群
scRNAsub <- scRNAsub[,scRNAsub$T_celltype %in% c('CD8+ Tem','CD8+ Tcm','CD8+ Trm')]

#添加细胞亚型列
scRNAsub$cell_type_val<-Idents(scRNAsub)
scRNAsub@meta.data[1:5,]
scRNAsub <- FindVariableFeatures(scRNAsub)
#提取基因表达矩阵
matrix <- as.matrix(scRNAsub@assays$RNA@counts)
#设置基因注释
gene_ann <- data.frame(
  gene_short_name = row.names(matrix), 
  row.names = row.names(matrix))
#设置细胞注释
sample_ann <- scRNAsub@meta.data

#标准流程
fd <- new("AnnotatedDataFrame",data=gene_ann)
pd<-new("AnnotatedDataFrame",data=sample_ann)
sc_cds_2 <- newCellDataSet(matrix,phenoData = pd,featureData =fd,expressionFamily = negbinomial.size(),lowerDetectionLimit=0.01)
sc_cds_2 <- estimateSizeFactors(sc_cds_2) #size factor：对细胞间mRNA的差异进行归一化，去除由于样本批次带来的影响（去批次效应）
sc_cds_2 <- estimateDispersions(sc_cds_2) #dispersion：帮助我们以后进行差异表达分析
sc_cds_2 <- detectGenes(sc_cds_2, min_expr = 1) #计算表达特定基因的细胞数量
#保存数据集中至少5个细胞中表达的基因
expressed_genes <- row.names(subset(fData(sc_cds_2), num_cells_expressed >= 5))
fData(sc_cds_2)[1:5,]

rm(matrix)

####选择决定细胞进程的基因####
##方法一：以高变基因为ordering
#ordering_genes<-scRNAsub@assays$RNA@var.features

##方法二：以亚群之间差异基因为ordering gene
diff_test_res <- differentialGeneTest(sc_cds_2[expressed_genes,],
                                      fullModelFormulaStr = "~cell_type_val",cores = 4,
                                      verbose = T) #+num_genes_expressed+orig.ident
ordering_genes <- row.names(subset(diff_test_res, qval < 1e-10)) #1e-200

##根据ordering_genes设置细胞顺序
sc_cds2 <- setOrderingFilter(sc_cds_2, ordering_genes)
##绘制ordering_genes散点图
plot_ordering_genes(sc_cds2)


####细胞轨迹构建####
##降维
sc_cds2<- reduceDimension(sc_cds2, max_components = 2, num_dim=6,reduction_method  = "DDRTree")#,residualModelFormulaStr = "~orig.ident")
##对细胞进行排序
sc_cds2 <- orderCells(sc_cds2)
#查看拟时序数据
pData(sc_cds2)[1:5,] #新增Pseudotime，State
#保存sc_cds2
save(sc_cds2, file='trajectory_object_CD8T.rda')
#load('trajectory_object.rda')


##细胞拟时序图（按细胞类型，Pseudotime，State分类绘制）
pdf(file=paste0('拟时序图_cell_type.pdf'),width = 4,height = 4)
plot_cell_trajectory(sc_cds2, color_by = "cell_type_val",cell_size = 0.5,show_branch_points = T) + scale_color_npg() + scale_fill_npg()
dev.off()
pdf(file=paste0('拟时序图_State.pdf'),width = 4,height = 4)
plot_cell_trajectory(sc_cds2, color_by = "State",cell_size = 0.5,show_branch_points = T) + scale_color_npg() + scale_fill_npg()
dev.off()
pdf(file=paste0('拟时序图_Pseudotime.pdf'),width = 4,height = 4)
plot_cell_trajectory(sc_cds2, color_by = "Pseudotime",cell_size = 0.5,show_branch_points = T)
dev.off()
pdf(file=paste0('拟时序图_group.pdf'),width = 4,height = 4)
plot_cell_trajectory(sc_cds2, color_by = "group",cell_size = 0.5,show_branch_points = T) + scale_color_npg() + scale_fill_npg()
dev.off()


#基因表达拟时序图（自己选择基因展示）
  #可以选择时间变化最显著的基因、关键细胞亚群的marker、组间关键DEG
pdf(file=paste0('拟时序图_关键基因.pdf'),width = 4,height = 4)
plot_cell_trajectory(sc_cds2, markers = c("ZNF683"),
                     use_color_gradient = T,show_branch_points = F) + 
  scale_color_gradient2(low="#4DBBD5FF",high="#DC0000FF",midpoint = 0)
dev.off()

#基因表达拟时序散点图（自己选择基因展示）
to_be_tested <- row.names(subset(fData(sc_cds2), 
                                 gene_short_name %in% c("ATPIF1","ATP5F1")))
cds_subset <- sc_cds2[to_be_tested,]

pdf(file=paste0('拟时序散点图_关键基因.pdf'),width = 5,height = 3)
plot_genes_in_pseudotime(cds_subset, color_by = "cell_type_val") + scale_color_npg() + scale_fill_npg()
dev.off()


####按照拟时间值找到差异表达基因####
diff_test_res <- differentialGeneTest(sc_cds2[expressed_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 2) #
# 获取显著差异基因名
# sig_gene_names <- row.names(subset(diff_test_res, qval < 0.05))
# 取前50个显著基因
sig_gene_names <- rownames(diff_test_res[order(diff_test_res$qval, decreasing = FALSE), ])[1:50]

#选取qval前6的基因
phe=pData(sc_cds2)
head(phe)
table(phe$State,phe$group) 

library(dplyr)  
# 保存前6的基因
diff_test_res %>% arrange(qval) %>% head() %>% dplyr::select(gene_short_name) -> my_pseudotime_gene
my_pseudotime_gene=my_pseudotime_gene[,1]
my_pseudotime_gene


##拟时序差异基因热图
dev.new()
pdf(file=paste0('拟时序热图.pdf'),width = 8,height = 10)
pseudoplot <- plot_pseudotime_heatmap(sc_cds2[sig_gene_names,],
                                      num_clusters = 3,###进行调节分群
                                      cores = 3,
                                      show_rownames = T,return_heatmap = T)
dev.off()

##获取热图中每一个cluster的基因名
clusters <- cutree(pseudoplot$tree_row, k = 3)
#clustering储存基因名
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)

write.csv(clustering, file='拟时序差异基因cluster.csv')


####识别具有分支依赖性表达的基因####
#计算建模分支节点(此处以节点1为例，可以做其他节点)
BEAM_branch1 <- BEAM(sc_cds2, branch_point = 1, cores = 4, #branch_point可改
                     progenitor_method = 'duplicate') 
BEAM_branch1 <- BEAM_branch1[order(BEAM_branch1$qval),]
BEAM_branch1 <- BEAM_branch1[,c("gene_short_name", "pval", "qval")]
head(BEAM_branch1) 

BEAM_res = BEAM_branch1
my_branched_heatmap <- plot_genes_branched_heatmap(sc_cds2[row.names(subset(BEAM_res, qval < 1e-4)),], 
                                                   branch_point = 1, #branch_point可改
                                                   num_clusters = 4, #num_clusters可改
                                                   e_short_name = T,show_rownames = T,return_heatmap = T)
pdf('monocle_BEAM_branch1_heatmap.pdf')
print(my_branched_heatmap$ph)
dev.off()

#将所做热图的基因和cluster提取出来
head(my_branched_heatmap$annotation_row)
table(my_branched_heatmap$annotation_row$Cluster) 
my_row <- my_branched_heatmap$annotation_row
my_row <- data.frame(cluster = my_row$Cluster,
                     gene = row.names(my_row),
                     stringsAsFactors = FALSE)
head(my_row[my_row$cluster == 3,'gene']) 

my_gene <- row.names(subset(fData(sc_cds2),
                            gene_short_name %in% head(my_row[my_row$cluster == 2,'gene'])))

# 绘制分支处的基因拟时序轨迹
#分支1
plot_genes_branched_pseudotime(sc_cds2[my_gene,],
                               branch_point = 1,
                               ncol = 1)



####二、黑素细胞####
##流程：找ordering gene，拟时序分析，得到各branch差异基因，分别功能富集分析

####准备文件####
scRNAsub = readRDS('scRNA_Mela.rds')
library(monocle)
Idents(scRNAsub)='Mela_celltype0'
#添加细胞亚型列
scRNAsub$cell_type_val<-Idents(scRNAsub)
scRNAsub@meta.data[1:5,]
scRNAsub <- FindVariableFeatures(scRNAsub)
#提取基因表达矩阵
matrix <- as.matrix(scRNAsub@assays$RNA@counts)
#设置基因注释
gene_ann <- data.frame(
  gene_short_name = row.names(matrix), 
  row.names = row.names(matrix))
#设置细胞注释
sample_ann <- scRNAsub@meta.data

#标准流程
fd <- new("AnnotatedDataFrame",data=gene_ann)
pd<-new("AnnotatedDataFrame",data=sample_ann)
sc_cds_2 <- newCellDataSet(matrix,phenoData = pd,featureData =fd,expressionFamily = negbinomial.size(),lowerDetectionLimit=0.01)
sc_cds_2 <- estimateSizeFactors(sc_cds_2) #size factor：对细胞间mRNA的差异进行归一化，去除由于样本批次带来的影响（去批次效应）
sc_cds_2 <- estimateDispersions(sc_cds_2) #dispersion：帮助我们以后进行差异表达分析
sc_cds_2 <- detectGenes(sc_cds_2, min_expr = 1) #计算表达特定基因的细胞数量
#保存数据集中至少5个细胞中表达的基因
expressed_genes <- row.names(subset(fData(sc_cds_2), num_cells_expressed >= 5))
fData(sc_cds_2)[1:5,]

rm(matrix)

####选择决定细胞进程的基因####
##方法一：以高变基因为ordering
#ordering_genes<-scRNAsub@assays$RNA@var.features

##方法二：以亚群之间差异基因为ordering gene
diff_test_res <- differentialGeneTest(sc_cds_2[expressed_genes,],
                                      fullModelFormulaStr = "~cell_type_val",cores = 4,
                                      verbose = T) #+num_genes_expressed+orig.ident
ordering_genes <- row.names(subset(diff_test_res, qval < 1e-200)) #1e-200

##根据ordering_genes设置细胞顺序
sc_cds2 <- setOrderingFilter(sc_cds_2, ordering_genes)
##绘制ordering_genes散点图
plot_ordering_genes(sc_cds2)


####细胞轨迹构建####
##降维
sc_cds2<- reduceDimension(sc_cds2, max_components = 2, num_dim=6,reduction_method  = "DDRTree")#,residualModelFormulaStr = "~orig.ident")
##对细胞进行排序
sc_cds2 <- orderCells(sc_cds2)
#查看拟时序数据
pData(sc_cds2)[1:5,] #新增Pseudotime，State
#保存sc_cds2
save(sc_cds2, file='trajectory_object_Mela.rda')
load('trajectory_object.rda')


##细胞拟时序图（按细胞类型，Pseudotime，State分类绘制）
pdf(file=paste0('拟时序图_cell_type.pdf'),width = 4,height = 4)
plot_cell_trajectory(sc_cds2, color_by = "cell_type_val",cell_size = 0.5,show_branch_points = T) + scale_color_npg() + scale_fill_npg()
dev.off()
pdf(file=paste0('拟时序图_Phase.pdf'),width = 4,height = 4)
plot_cell_trajectory(sc_cds2, color_by = "Phase",cell_size = 0.5,show_branch_points = T) + scale_color_npg() + scale_fill_npg()
dev.off()
pdf(file=paste0('拟时序图_Pseudotime.pdf'),width = 4,height = 4)
plot_cell_trajectory(sc_cds2, color_by = "Pseudotime",cell_size = 0.5,show_branch_points = T)
dev.off()
pdf(file=paste0('拟时序图_group.pdf'),width = 4,height = 4)
plot_cell_trajectory(sc_cds2, color_by = "group",cell_size = 0.5,show_branch_points = T) + scale_color_npg() + scale_fill_npg()
dev.off()

#基因表达拟时序图（自己选择基因展示）
  #可以选择时间变化最显著的基因、关键细胞亚群的marker、组间关键DEG
   #TYR,MITF代表色素沉着;SERPINE2,TIMP1代表EMT
pdf(file=paste0('trajectory_TIMP1.pdf'),width = 4,height = 4)
plot_cell_trajectory(sc_cds2, markers = c("TIMP1"), 
                     use_color_gradient = T,show_branch_points = F) #+ scale_color_gradient2(low="#4DBBD5FF",high="#DC0000FF",midpoint = 0)
dev.off()

#基因表达拟时序散点图（自己选择基因展示）
to_be_tested <- row.names(subset(fData(sc_cds2), 
                                 gene_short_name %in% c('TYR','MITF','SERPINE2','TIMP1')))
cds_subset <- sc_cds2[to_be_tested,]

pdf(file=paste0('trajectory_scatter_KeyGenes.pdf'),width = 4,height = 4)
plot_genes_in_pseudotime(cds_subset, color_by = "cell_type_val",cell_size = 0.5) + scale_color_npg() + scale_fill_npg()
dev.off()


####按照拟时间值找到差异表达基因####
diff_test_res <- differentialGeneTest(sc_cds2[expressed_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 2) #
# 获取显著差异基因名
# sig_gene_names <- row.names(subset(diff_test_res, qval < 0.05))
# 取前100个显著基因
sig_gene_names <- rownames(diff_test_res[order(diff_test_res$qval, decreasing = FALSE), ])[1:50]

#选取qval前6的基因
phe=pData(sc_cds2)  #获取每个细胞的pseudotime/state
head(phe)
table(phe$State,phe$group) 
library(dplyr)  
diff_test_res %>% arrange(qval) %>% head() 
# 保存前6的基因
diff_test_res %>% arrange(qval) %>% head() %>% dplyr::select(gene_short_name) -> my_pseudotime_gene
my_pseudotime_gene=my_pseudotime_gene[,1]
my_pseudotime_gene


##拟时序差异基因热图
dev.new()
pdf(file=paste0('拟时序热图.pdf'),width = 4,height = 5)
pseudoplot <- plot_pseudotime_heatmap(sc_cds2[sig_gene_names,],
                                      num_clusters = 4,###进行调节分群
                                      cores = 3,
                                      show_rownames = T,return_heatmap = T)
dev.off()

##获取热图中每一个cluster的基因名
clusters <- cutree(pseudoplot$tree_row, k = 4)
#clustering储存基因名
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)

write.csv(clustering, file='拟时序差异基因cluster.csv')




####三、单核吞噬细胞####
##流程：找ordering gene，拟时序分析，得到各branch差异基因，分别功能富集分析

####准备文件####
scRNAsub = readRDS('scRNA_Mf.rds')
library(monocle)
#根据实际需求选择要分析的细胞亚群
scRNAsub <- scRNAsub[,Idents(scRNAsub)%in%c("Mφ_APOE","Mφ_CCL3","Mono_FCN1")]

#添加细胞亚型列
scRNAsub$cell_type_val<-Idents(scRNAsub)
scRNAsub@meta.data[1:5,]
scRNAsub <- FindVariableFeatures(scRNAsub)
#提取基因表达矩阵
matrix <- as.matrix(scRNAsub@assays$RNA@counts)
#设置基因注释
gene_ann <- data.frame(
  gene_short_name = row.names(matrix), 
  row.names = row.names(matrix))
#设置细胞注释
sample_ann <- scRNAsub@meta.data

#标准流程
fd <- new("AnnotatedDataFrame",data=gene_ann)
pd<-new("AnnotatedDataFrame",data=sample_ann)
sc_cds_2 <- newCellDataSet(matrix,phenoData = pd,featureData =fd,expressionFamily = negbinomial.size(),lowerDetectionLimit=0.01)
sc_cds_2 <- estimateSizeFactors(sc_cds_2) #size factor：对细胞间mRNA的差异进行归一化，去除由于样本批次带来的影响（去批次效应）
sc_cds_2 <- estimateDispersions(sc_cds_2) #dispersion：帮助我们以后进行差异表达分析
sc_cds_2 <- detectGenes(sc_cds_2, min_expr = 1) #计算表达特定基因的细胞数量
#保存数据集中至少5个细胞中表达的基因
expressed_genes <- row.names(subset(fData(sc_cds_2), num_cells_expressed >= 5))
fData(sc_cds_2)[1:5,]

rm(matrix)

####选择决定细胞进程的基因####
##方法一：以高变基因为ordering
#ordering_genes<-scRNAsub@assays$RNA@var.features

##方法二：以亚群之间差异基因为ordering gene
diff_test_res <- differentialGeneTest(sc_cds_2[expressed_genes,],
                                      fullModelFormulaStr = "~cell_type_val",cores = 4,
                                      verbose = T) #+num_genes_expressed+orig.ident
ordering_genes <- row.names(subset(diff_test_res, qval < 1e-10)) #1e-200

##根据ordering_genes设置细胞顺序
sc_cds2 <- setOrderingFilter(sc_cds_2, ordering_genes)
##绘制ordering_genes散点图
plot_ordering_genes(sc_cds2)


####细胞轨迹构建####
##降维
sc_cds2<- reduceDimension(sc_cds2, max_components = 2, num_dim=6,reduction_method  = "DDRTree")#,residualModelFormulaStr = "~orig.ident")
##对细胞进行排序
sc_cds2 <- orderCells(sc_cds2)
#查看拟时序数据
pData(sc_cds2)[1:5,] #新增Pseudotime，State
#保存sc_cds2
save(sc_cds2, file='trajectory_object_Mf.rda')
load('trajectory_object.rda')


##细胞拟时序图（按细胞类型，Pseudotime，State分类绘制）
pdf(file=paste0('拟时序图_cell_type.pdf'),width = 4,height =4)
plot_cell_trajectory(sc_cds2, color_by = "cell_type_val",cell_size = 0.5,show_branch_points = T) + scale_color_npg() + scale_fill_npg()
dev.off()
pdf(file=paste0('拟时序图_State.pdf'),width = 4,height =4)
plot_cell_trajectory(sc_cds2, color_by = "State",cell_size = 0.5,show_branch_points = T) + scale_color_npg() + scale_fill_npg()
dev.off()
pdf(file=paste0('拟时序图_Pseudotime.pdf'),width = 4,height =4)
plot_cell_trajectory(sc_cds2, color_by = "Pseudotime",cell_size = 0.5,show_branch_points = T)
dev.off()
pdf(file=paste0('拟时序图_group.pdf'),width = 4,height =4)
plot_cell_trajectory(sc_cds2, color_by = "group",cell_size = 0.5,show_branch_points = T) + scale_color_npg() + scale_fill_npg()
dev.off()

#基因表达拟时序图（自己选择基因展示）
  #可以选择时间变化最显著的基因、关键细胞亚群的marker、组间关键DEG
pdf(file=paste0('拟时序图_关键基因.pdf'),width = 4,height = 4)
plot_cell_trajectory(sc_cds2, markers = c("TXNIP"),
                     cell_size = 0.5,use_color_gradient = T,show_branch_points = F) + 
  scale_color_gradient2(low="#3C5488FF",high="#DC0000FF")
dev.off()

#基因表达拟时序散点图（自己选择基因展示）
to_be_tested <- row.names(subset(fData(sc_cds2), 
                                 gene_short_name %in% c("TXNIP")))
cds_subset <- sc_cds2[to_be_tested,]

pdf(file=paste0('拟时序散点图_关键基因.pdf'),width = 6,height = 4)
plot_genes_in_pseudotime(cds_subset, color_by = "cell_type_val") + scale_color_npg() + scale_fill_npg()
dev.off()


####按照拟时间值找到差异表达基因####
diff_test_res <- differentialGeneTest(sc_cds2[expressed_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 2) #
# 获取显著差异基因名
# sig_gene_names <- row.names(subset(diff_test_res, qval < 0.05))
# 取前50个显著基因
sig_gene_names <- rownames(diff_test_res[order(diff_test_res$qval, decreasing = FALSE), ])[1:50]

#选取qval前6的基因
phe=pData(sc_cds2)
head(phe)
table(phe$State,phe$group) 
diff_test_res %>% dplyr::arrange(qval) %>% head() 
# 保存前6的基因
my_pseudotime_gene = diff_test_res %>% arrange(qval) %>% head() %>% 
  dplyr::select(gene_short_name)
my_pseudotime_gene = my_pseudotime_gene[,1]
my_pseudotime_gene


##拟时序差异基因热图
dev.new()
pdf(file=paste0('拟时序热图.pdf'),width = 4,height = 5)
pseudoplot <- plot_pseudotime_heatmap(sc_cds2[sig_gene_names,],
                                      num_clusters = 4,###进行调节分群
                                      cores = 3,
                                      show_rownames = T,return_heatmap = T)
dev.off()

##获取热图中每一个cluster的基因名
clusters <- cutree(pseudoplot$tree_row, k = 4)
#clustering储存基因名
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)

write.csv(clustering, file='拟时序差异基因cluster.csv')
