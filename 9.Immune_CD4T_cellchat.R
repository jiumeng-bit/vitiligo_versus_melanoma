##健康、白癜风和黑素瘤各做一次
##每次做的内容：注释好的T与Mf、DC与NK、KRT做通讯
#先明确作用对象，再明确高表达通路，再明确配体-受体对
##最后进行细胞通讯差异分析
rm(list=ls())
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
library(stringr)


####准备工作####
#读取文件
setwd("E:/【科研学习】/【皮肤科】/白癜风课题/白癜风+黑色素瘤/思路3原始结果/QC_scRNA_data")
scRNA_T=readRDS('./scRNA_T&NK_new.rds')
scRNA_Myeloid=readRDS('./scRNA_Myeloid.rds')

##合并T细胞与其他细胞
scRNA_T$subcelltype = ifelse(grepl('Treg',scRNA_T$T_celltype),'Treg',
                             ifelse(grepl('CD4',scRNA_T$T_celltype),'CD4+ T','CD8+ T'))

scRNA_Myeloid$subcelltype = scRNA_Myeloid$Mono_celltype

scRNAsub <- merge(scRNA_T, y=c(scRNA_Myeloid))
Idents(scRNAsub) = scRNAsub$subcelltype
table(Idents(scRNAsub))

#将细胞按组织类型拆分
normaldata<-scRNAsub[,scRNAsub$group=="Healthy"]
tumordata<-scRNAsub[,scRNAsub$group=="Melanoma"]
vitidata<-scRNAsub[,scRNAsub$group=="Vitiligo"]
rm(scRNAsub)
gc()




####一、normal组织####
####细胞通讯计算####
library(CellChat)
meta =normaldata@meta.data # a dataframe with rownames containing cell mata data
gc()
data_input <- as.matrix(normaldata@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))

cellchat <- createCellChat(object = data_input, meta = meta, group.by = "subcelltype")

CellChatDB <- CellChatDB.human 
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 

dplyr::glimpse(CellChatDB$interaction)##配体-受体分析
# 提取数据库支持的数据子集
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)
unique(cellchat@idents)
cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

df.net<- subsetCommunication(cellchat)


####细胞通讯数量强度####
#aggregateNet
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
table(cellchat@idents)
#查看细胞通讯的数量/权重矩阵
mat <- cellchat@net$weight
#mat <- cellchat@net$weight
#绘制细胞通讯图（Ligand-receptor） #7*4

##针对其中一个细胞亚型分析其对其他细胞的interaction的Ligand-receptor数量
#创建一个空矩阵，将想要分析的细胞所在行填充进去
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[c(5,6,8),] <- mat[c(5,6,8),] #感兴趣的是第3行细胞（T细胞）
pdf(file=paste0('Healthy_CD8T_cellchat_weight.pdf'),width = 7,height = 7)
netVisual_circle(mat2, vertex.weight = groupSize,arrow.size = 0.2, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[5])
dev.off()

####信号通路通讯分析####
cellchat@netP$pathways   #信号通路查看
pathways.show <- c('MIF')   #以'MIF'信号通路展示为例
levels(cellchat@idents)   #查看细胞亚群及factor顺序
vertex.receiver = c(5:8)  #左侧列展示感兴趣的亚群
#层级图（Hierarchy plot）
pdf(file=paste0('Healthy_CD8T_MIF_层级图.pdf'),width = 10,height = 7)
netVisual_aggregate(cellchat,                #左侧列展示感兴趣的亚群
                    layout = c('hierarchy'), #"circle", "hierarchy", "chord"
                    signaling = pathways.show,
                    vertex.receiver = vertex.receiver)
dev.off()

#计算配受体对在目标信号通路中的贡献条形图
pdf(file=paste0('Healthy_CD8T_MIF_贡献条形图.pdf'),width = 4,height = 3)
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()
#提取细胞对
pairLR.CXCL <- extractEnrichedLR(cellchat,
                                 signaling = pathways.show,
                                 geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] #以贡献度top1的配受体对为例
pairLR.CXCL
LR.show
#Hierarchy plot
netVisual_individual(cellchat,
                     layout = c('hierarchy'),
                     signaling = pathways.show, #目标信号通路
                     pairLR.use = LR.show, #目标配受体对
                     vertex.receiver = vertex.receiver) #感兴趣的细胞亚群


####配体-受体通讯情况####
levels(cellchat@idents) #查看有哪些细胞类型
cellchat@netP$pathways   #信号通路查看
#sources.use为发出信号的细胞，targets.use为接受信号的细胞
#c()里面的数字与levels(cellchat@idents)中细胞的位次对应
##指定信号通路绘制气泡图
pdf(file=paste0('Healthy_CD8T_bubble.pdf'),width = 8,height = 3)
netVisual_bubble(cellchat,
                 sources.use = c(9:11,16,12:14),
                 targets.use = c(5,6,8),
                 #signaling = c("CCL","CXCL","IFN-II"), #指定信号通路
                 remove.isolate = FALSE) +
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=1))
dev.off()
##指定信号通路绘制气泡图
pdf(file=paste0('Healthy_CD8Tem_bubble.pdf'),width = 5,height = 4)
netVisual_bubble(cellchat,
                 sources.use = c(9:11,16,12:14),
                 targets.use = c(5,6,8),
                 signaling = c("CCL","CXCL","IFN-II"), #指定信号通路
                 remove.isolate = FALSE) +
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=1))
dev.off()

#指定配受体对绘制气泡图
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","IFN-II")) #确定在目标信号通路中有重要作用的配受体对
pairLR.use
netVisual_bubble(cellchat,
                 sources.use = c(5,7),
                 targets.use = c(7:12),
                 pairLR.use = pairLR.use,
                 remove.isolate = TRUE) + coord_flip()

#参与目标信号通路的基因在各细胞亚群的表达分布展示
pdf(file=paste0('Healthy_CD8T_CCL_expression.pdf'),width = 6,height = 4.5)
plotGeneExpression(cellchat, signaling = 'CCL', type = 'violin') #小提琴图
dev.off()

saveRDS(cellchat,file="healthy_cellchat_analysis.rds")




####二、tumor组织####
####细胞通讯计算####
library(CellChat)
meta =tumordata@meta.data # a dataframe with rownames containing cell mata data
gc()
data_input <- as.matrix(tumordata@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))

cellchat <- createCellChat(object = data_input, meta = meta, group.by = "subcelltype")

CellChatDB <- CellChatDB.human 
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 

dplyr::glimpse(CellChatDB$interaction)##配体-受体分析
# 提取数据库支持的数据子集
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)
unique(cellchat@idents)
cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

df.net<- subsetCommunication(cellchat)


####细胞通讯数量强度####
#aggregateNet
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
table(cellchat@idents)
#查看细胞通讯的数量/权重矩阵
mat <- cellchat@net$weight
#mat <- cellchat@net$count
#绘制细胞通讯图（Ligand-receptor）
  #针对其中一个细胞亚型分析其对其他细胞的interaction的Ligand-receptor数量
  #创建一个空矩阵，将想要分析的细胞所在行填充进去
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[,c(1)] <- mat[,c(1)]  #行代表信号发出的细胞，列代表接收信号的细胞
pdf(file=paste0('Melanoma_CD8Trm_cellchat_weight.pdf'),width = 6,height = 6)
netVisual_circle(mat2, vertex.weight = groupSize,arrow.size = 0.2, weight.scale = T, edge.weight.max = max(mat))
dev.off()

####信号通路通讯分析####
cellchat@netP$pathways   #信号通路查看
pathways.show <- c('MIF')   #以'MIF'信号通路展示为例
levels(cellchat@idents)   #查看细胞亚群及factor顺序
vertex.receiver = c(1:6)  #左侧列展示感兴趣的亚群
#层级图（Hierarchy plot）
pdf(file=paste0('Melanoma_CD8T_MIF_层级图.pdf'),width = 10,height = 7)
netVisual_aggregate(cellchat,                #左侧列展示感兴趣的亚群
                    layout = c('hierarchy'), #"circle", "hierarchy", "chord"
                    signaling = pathways.show,
                    vertex.receiver = vertex.receiver)
dev.off()

#计算配受体对在目标信号通路中的贡献条形图
pdf(file=paste0('Melanoma_CD8T_MIF_贡献条形图.pdf'),width = 4,height = 3)
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()
#提取细胞对
pairLR.CXCL <- extractEnrichedLR(cellchat,
                                 signaling = pathways.show,
                                 geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] #以贡献度top1的配受体对为例
pairLR.CXCL
LR.show
#Hierarchy plot
netVisual_individual(cellchat,
                     layout = c('hierarchy'),
                     signaling = pathways.show, #目标信号通路
                     pairLR.use = LR.show, #目标配受体对
                     vertex.receiver = vertex.receiver) #感兴趣的细胞亚群


####配体-受体通讯情况####
levels(cellchat@idents) #查看有哪些细胞类型
cellchat@netP$pathways   #信号通路查看
#sources.use为发出信号的细胞，targets.use为接受信号的细胞
#c()里面的数字与levels(cellchat@idents)中细胞的位次对应
##指定信号通路绘制气泡图
pdf(file=paste0('Melanoma_CD4T_bubble.pdf'),width = 3.5,height = 1.5)
netVisual_bubble(cellchat,
                 sources.use = c(2:5),
                 targets.use = c(1),
                 #signaling = c("CCL","CXCL","IFN-II"), #指定信号通路
                 remove.isolate = FALSE) +
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=1))
dev.off()
##指定信号通路绘制气泡图
pdf(file=paste0('Melanoma_CD8Tem_bubble.pdf'),width = 5,height = 4)
netVisual_bubble(cellchat,
                 sources.use = c(1:5,7:14),
                 targets.use = c(6),
                 #signaling = c("CCL","CXCL","IFN-II"), #指定信号通路
                 remove.isolate = FALSE) + coord_flip()
dev.off()

#指定配受体对绘制气泡图
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","IFN-II")) #确定在目标信号通路中有重要作用的配受体对
pairLR.use
netVisual_bubble(cellchat,
                 sources.use = c(6),
                 targets.use = c(1:5,7:15),
                 pairLR.use = pairLR.use,
                 remove.isolate = TRUE) + coord_flip()

#参与目标信号通路的基因在各细胞亚群的表达分布展示
pdf(file=paste0('Melanoma_CD8T_MIF_expression.pdf'),width = 6,height = 4.5)
plotGeneExpression(cellchat, signaling = 'MIF', type = 'violin') #小提琴图
dev.off()

saveRDS(cellchat,file="melanoma_cellchat_analysis.rds")




####三、白癜风组织####
####细胞通讯计算####
library(CellChat)
meta =vitidata@meta.data # a dataframe with rownames containing cell mata data
gc()
data_input <- as.matrix(vitidata@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))

cellchat <- createCellChat(object = data_input, meta = meta, group.by = "subcelltype")

CellChatDB <- CellChatDB.human 
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 

dplyr::glimpse(CellChatDB$interaction)##配体-受体分析
# 提取数据库支持的数据子集
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)
unique(cellchat@idents)
cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

df.net<- subsetCommunication(cellchat)


####细胞通讯数量强度####
#aggregateNet
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
table(cellchat@idents)
#查看细胞通讯的数量/权重矩阵
mat <- cellchat@net$weight
#mat <- cellchat@net$count
#绘制细胞通讯图（Ligand-receptor）
##针对其中一个细胞亚型分析其对其他细胞的interaction的Ligand-receptor数量
#创建一个空矩阵，将想要分析的细胞所在行填充进去
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[,c(2)] <- mat[,c(2)] #感兴趣的是第3行细胞（T细胞）
pdf(file=paste0('Vitiligo_Treg_cellchat_weight.pdf'),width = 6,height = 6)
netVisual_circle(mat2, vertex.weight = groupSize,arrow.size = 0.2, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[5])
dev.off()

####信号通路通讯分析####
cellchat@netP$pathways   #信号通路查看
pathways.show <- c('MIF')   #以'MIF'信号通路展示为例
levels(cellchat@idents)   #查看细胞亚群及factor顺序
vertex.receiver = c(2,3,9)  #左侧列展示感兴趣的亚群
#层级图（Hierarchy plot）
pdf(file=paste0('Vitiligo_T_MIF_hierarchy.pdf'),width = 10,height = 7)
netVisual_aggregate(cellchat,                #左侧列展示感兴趣的亚群
                    layout = c('hierarchy'), #"circle", "hierarchy", "chord"
                    signaling = pathways.show,
                    vertex.receiver = vertex.receiver)
dev.off()
#和弦图
netVisual_aggregate(cellchat,                #左侧列展示感兴趣的亚群
                    layout = c('chord'), #"circle", "hierarchy", "chord"
                    signaling = pathways.show,
                    #sources.use = c(1,3:9), 
                    targets.use = c(2,3,9), 
                    vertex.receiver = vertex.receiver)


#计算配受体对在目标信号通路中的贡献条形图
pdf(file=paste0('Vitiligo_CD8T_MIF_贡献条形图.pdf'),width = 4,height = 3)
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()
#提取细胞对
pairLR.CXCL <- extractEnrichedLR(cellchat,
                                 signaling = pathways.show,
                                 geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] #以贡献度top1的配受体对为例
pairLR.CXCL
LR.show
#Hierarchy plot
netVisual_individual(cellchat,
                     layout = c('hierarchy'),
                     signaling = pathways.show, #目标信号通路
                     pairLR.use = LR.show, #目标配受体对
                     vertex.receiver = vertex.receiver) #感兴趣的细胞亚群


####配体-受体通讯情况####
levels(cellchat@idents) #查看有哪些细胞类型
cellchat@netP$pathways   #信号通路查看
#sources.use为发出信号的细胞，targets.use为接受信号的细胞
#c()里面的数字与levels(cellchat@idents)中细胞的位次对应
##指定信号通路绘制气泡图
pdf(file=paste0('Vitiligo_CD4T_bubble.pdf'),width = 4,height = 2)
netVisual_bubble(cellchat,
                 sources.use = c(6:8,4:5,1),
                 targets.use = c(2),
                 #signaling = c("CCL","CXCL","IFN-II"), #指定信号通路
                 remove.isolate = FALSE) +
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=1))
dev.off()
##指定信号通路绘制气泡图
pdf(file=paste0('Vitiligo_CD8Tem_bubble.pdf'),width = 5,height = 4)
netVisual_bubble(cellchat,
                 sources.use = c(8:10,6:7,1),
                 targets.use = c(3:5),
                 signaling = c("CCL","CXCL","IFN-II"), #指定信号通路
                 remove.isolate = FALSE) + coord_flip()
dev.off()

#指定配受体对绘制气泡图
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","IFN-II")) #确定在目标信号通路中有重要作用的配受体对
pairLR.use
netVisual_bubble(cellchat,
                 sources.use = c(6),
                 targets.use = c(1:5,7:15),
                 pairLR.use = pairLR.use,
                 remove.isolate = TRUE) + coord_flip()

#参与目标信号通路的基因在各细胞亚群的表达分布展示
pdf(file=paste0('Vitiligo_CD8T_CCL_expression.pdf'),width = 6,height = 4.5)
plotGeneExpression(cellchat, signaling = 'CCL', type = 'violin') #小提琴图
dev.off()

saveRDS(cellchat,file="vitiligo_cellchat_analysis.rds")
