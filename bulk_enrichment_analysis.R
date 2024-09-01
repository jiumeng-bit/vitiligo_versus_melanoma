library(stringr)

setwd("E:/【科研学习】/【皮肤科】/白癜风课题/白癜风+黑色素瘤/原始结果")
#define the color
library(ggsci)
cors <- pal_npg()(10) #定义颜色

####获取交集基因####
###共有DEG
DEG_Viti = read.table("./1.差异分析/白癜风/DEG.txt")
DEG_Mel = read.table("./1.差异分析/黑素瘤/DEG.txt")

DEG_Viti_UP = subset(DEG_Viti, change == 'UP')
DEG_Viti_DOWN = subset(DEG_Viti, change == 'DOWN')
DEG_Mel_UP = subset(DEG_Mel, change == 'UP')
DEG_Mel_DOWN = subset(DEG_Mel, change == 'DOWN')

#变化趋势相同的DEG
same_DEG = union(intersect(rownames(DEG_Viti_UP), rownames(DEG_Mel_UP)), 
                 intersect(rownames(DEG_Viti_DOWN), rownames(DEG_Mel_DOWN)))
#变化趋势相反的DEG
inverse_DEG = union(intersect(rownames(DEG_Viti_UP), rownames(DEG_Mel_DOWN)), 
                    intersect(rownames(DEG_Viti_DOWN), rownames(DEG_Mel_UP)))

###共有Module genes
#白癜风
module_Viti_turquoise = read.table("./2.WGCNA/白癜风/turquoise.txt")
module_Viti_blue = read.table("./2.WGCNA/白癜风/blue.txt")
module_Viti = rbind(module_Viti_turquoise, module_Viti_blue)
#黑素瘤
module_Mel_blue = read.table("./2.WGCNA/黑素瘤/blue.txt")
module_Mel_grey = read.table("./2.WGCNA/黑素瘤/grey.txt")

#变化趋势相同的module
same_module = intersect(module_Viti$V1, module_Mel_grey$V1)
#变化趋势相反的module
inverse_module = intersect(module_Viti$V1, module_Mel_blue$V1)

##common DEG与common module取并集
same_gene = union(same_DEG, same_module)
inverse_gene = union(inverse_DEG, inverse_module)

setwd("E:/【科研学习】/【皮肤科】/白癜风课题/白癜风+黑色素瘤/原始结果/3.富集分析")
save(same_gene, inverse_gene, file = 'common_genes.rda')
write.csv(same_gene, file = 'same_gene.csv')
write.csv(inverse_gene, file = 'inverse_gene.csv')


##韦恩图
library(ggvenn)
venn_dat <- list(
  "DEG_Viti_UP" = rownames(DEG_Viti_UP),
  "DEG_Mel_DOWN" = rownames(DEG_Mel_DOWN),
  "DEG_Viti_DOWN" = rownames(DEG_Viti_DOWN),
  "DEG_Mel_UP" = rownames(DEG_Mel_UP))

ggvenn(
  data = venn_dat,         # 数据列表
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
  text_size = 4             # 交集个数文字大小
)



##富集分析
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)


#产生ENTREZID（富集分析用）
genelist <- bitr(inverse_gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')

####GO分析####
ego <- enrichGO(gene = genelist$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.99, 
                qvalueCutoff =0.99,
                readable = TRUE)
ego_res <- ego@result
ego_res$LogP <- -log(ego_res$p.adjust) #计算-logP
#只提取与xxx有关的基因集   JAK|IFN|cytokine|chemokine|oxidative|pathway|T cell|immune
selected_ego_res = ego_res  %>%
  dplyr::filter(grepl('pathway|signal|T cell|immune|oxidative|antigen|mitochondrion', Description)) %>%
  dplyr::arrange(dplyr::desc(LogP),dplyr::desc(Description)) %>%
  mutate(Description =forcats:: fct_inorder(Description))

# 气泡图
dotplot(ego, showCategory = 20)
# 分类展示
dotplot(ego,showCategory = 5, split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')
barplot(ego, drop = TRUE, showCategory =5,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')
# 个性化柱状图
pdf('inverseDEG_GO_BP_barplot.pdf', width = 7, height = 3.5)
ggplot(selected_ego_res[1:15,],aes(Count,Description))+ #只展示前十条通路
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=LogP),stat='identity')+
  scale_fill_gradient(low="#F39B7FFF",high="#DC0000FF")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=0.5))
dev.off()

####KEGG分析####
kk <- enrichKEGG(gene = genelist$ENTREZID,
                 keyType = 'kegg',
                 organism = 'hsa',
                 pvalueCutoff = 0.99,
                 qvalueCutoff =0.99)
kk_res <- kk@result
kk_res$LogP <- -log(kk_res$p.adjust) #计算-logP
#只提取与xxx有关的基因集  JAK|IFN|cytokine|chemokine|oxidative|pathway|T cell|immune
selected_kk_res = kk_res  %>%
  dplyr::filter(grepl('pathway|signal|T cell|immune|oxidative|antigen|mitochondrion', Description)) %>%
  dplyr::arrange(dplyr::desc(LogP),dplyr::desc(Description)) %>%
  mutate(Description =forcats:: fct_inorder(Description))

#气泡图
dotplot(kk, showCategory = 20)
# 个性化柱状图
pdf('inverseDEG_KEGG_barplot.pdf', width = 5, height = 2.5)
ggplot(selected_kk_res[1:10,],aes(Count,Description))+ #只展示前十条通路
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=LogP),stat='identity')+
  scale_fill_gradient(low="#F39B7FFF",high="#DC0000FF")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=0.5))
dev.off()

####关键通路GSEA####
#教程：https://zhuanlan.zhihu.com/p/667669372
#msigdbr包教程：https://www.jianshu.com/p/f2febb3123d8
library(clusterProfiler)
# 白癜风和黑素瘤各做一次，后续使用logFC进行基因排序
DEG_limma <- DEG_Mel
head(DEG_limma)
# 加载基因集，基因集介绍往下滑
geneSet = msigdbr(species = "Homo sapiens", category = "C2")  #category = "C2"
geneSet %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
geneSet = geneSet %>% dplyr::select(gs_name, gene_symbol)
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
                 'REACTOME_INTERFERON_GAMMA_SIGNALING', #BIOCARTA_IFNG_PATHWAY
                 'GOBP_RESPONSE_TO_OXIDATIVE_STRESS',
                 #'GOBP_T_CELL_MEDIATED_CYTOTOXICITY',
                 'GOBP_AUTOPHAGY_OF_MITOCHONDRION',
                 'GOBP_RESPONSE_TO_TYPE_I_INTERFERON',
                 'WP_FERROPTOSIS',
                 'BIOCARTA_CYTOKINE_PATHWAY')
geneSet = subset(geneSet, gs_name %in% geneSet_name)
head(geneSet)

# 接下来我们进行基因排序
geneList <- DEG_limma$log2FoldChange      # 获取GeneList
names(geneList) <- rownames(DEG_limma)      # 对GeneList命名
geneList <- sort(geneList, decreasing = T)  # 从高到低排序

# 开始GSEA富集分析
GSEA_enrichment <- GSEA(geneList,                 # 排序后的gene
                        TERM2GENE = geneSet,      # 基因集
                        pvalueCutoff = 0.99,      # P值阈值
                        minGSSize = 10,           # 最小基因数量
                        maxGSSize = 1000,         # 最大基因数量
                        pAdjustMethod = "BH")     # 校正P值的计算方法

result <- data.frame(GSEA_enrichment)

#只提取感兴趣的基因集
selected_result=result  %>%
  dplyr::filter(stringr::str_detect(pattern = "CD8_TCELL",Description)) %>%
  group_by(Description) %>% add_count() %>%
  dplyr::arrange(dplyr::desc(n),dplyr::desc(Description)) %>%
  mutate(Description =forcats:: fct_inorder(Description))

# 特定通路绘图
library(enrichplot) 
#gseaplot2(GSEA_enrichment, c(2), color = "red3", pvalue_table = T)
pdf('3.mela_Chemokine_GSEA.pdf', width = 5, height = 4)
gseaplot2(GSEA_enrichment,
          geneSetID = 'KEGG_CHEMOKINE_SIGNALING_PATHWAY',
          color = cors,
          title = 'Chemokine signaling',
          rel_heights = c(1.3, 0.3, 0.6),
          pvalue_table = F)
dev.off()

# 展示富集到的通路，我们这里选择展示前15个
dotplot(GSEA_enrichment, showCategory = 15, color = "p.adjust")

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



####GSVA####
#教程：https://www.jianshu.com/p/bc1715b1ff4e
library(GSVA)

##选择msigdb基因集
#读取表达矩阵
load('GSE_merge.rda')
exp <- exp1_merge
#读取基因集
geneSet = msigdbr(species = "Homo sapiens")  #category = "C2"
geneSet %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
geneSet = subset(geneSet, gs_name %in% geneSet_name)
geneSet = geneSet %>% split(x = .$gene_symbol, f = .$gs_name)#基因集是list
str(head(geneSet))

##GSVA分析
exp_gsva <- gsva(as.matrix(exp), geneSet) #表达矩阵需转换为matrix格式
#method=c("gsva", "ssgsea", "zscore", "plage")

##热图
library(pheatmap) 
#设置参考水平
group_list <- c(rep('Normal',25),rep('Vitiligo',25))
group_list = factor(group_list,levels = c("Normal","Vitiligo"))

#sim_clinical=rbind(subset(sim_clinical,group=='Healthy'),subset(sim_clinical,group=='Melanoma'))
#counts=counts[,rownames(sim_clinical)]
#identical(colnames(counts),rownames(sim_clinical))
#group_list <- c(rep('Healthy',356),rep('Melanoma',98))
#group_list = factor(group_list,levels = c("Healthy","Melanoma"))

annotation_col=data.frame(group=group_list)
rownames(annotation_col) <- colnames(exp_gsva)

pdf('melanoma_GSVA_interest.pdf', width = 6, height = 4)
pheatmap::pheatmap(exp_gsva,
                   show_colnames = F, 
                   cluster_rows = F, cluster_cols = F, 
                   color = colorRampPalette(c("#3C5488FF", "white", "#DC0000FF"))(50),
                   annotation_col = annotation_col)
dev.off()

##自制基因集
gs = read.csv("geneset.csv", stringsAsFactors = FALSE, check.names = FALSE) #geneset.csv每列首行是基因集名称，基因集包含的基因列在名称下面
a = read.table("RNA.csv", stringsAsFactors = FALSE, header = TRUE, row.names = 1, sep = ",")
a = as.matrix(a)
gs = as.list(gs)
gs = lapply(gs, function(x) x[!is.na(x)])
ssgsea_score = gsva(a, gs, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'


####关键通路基因气泡图####
exp1_df <- as.data.frame(t(exp1_merge))
exp1_df$type <-rt2_merge$group
data_new0 <- aggregate(. ~ type, exp1_df, mean) #aggregate函数计算每种细胞类型的基因表达平均值
data_new <- as.data.frame(t(data_new0))
colnames(data_new) <- data_new[1,]
data_new <- data_new[-1,]
data_new <- data_new[genes,]
data_new <- data_new %>%
  mutate(across(where(is.character), as.numeric))

annotation_col = data.frame(group=factor(c('Healthy','Vitiligo')))
rownames(annotation_col) <- colnames(data_new)
pheatmap::pheatmap(as.matrix(data_new),
                   show_colnames = F, 
                   cluster_rows = F, cluster_cols = F, 
                   color = colorRampPalette(c("#3C5488FF", "white", "#DC0000FF"))(50),
                   annotation_col = annotation_col)




####细胞特征生存分析####

#读取表达矩阵
load(file="TCGA_SKCM.rda")
exp <- log(counts+1)
##准备细胞特征基因集
geneSet <- list('Keratinocyte'=c("KRT1","KRT10","KRT14","KRT15"),
               'Melanocyte'=c("DCT","PMEL","TYRP1","MLANA"),
               'Fibroblast'=c("COL1A1","DCN","SFRP2","TWIST2"),
               'Endothelial'=c("PECAM1","CLEC14A","AQP1","ECSCR.1"),
               'T cell'=c("TRAC","CD3D","TRBC2","CD3E"),
               'NK cell'=c("NKG7","GNLY"),
               'Langerhans'=c("CD207","CD1A","FCGBP","S100B"),
               'Mono phagocyte'=c("LYZ","CD1C","IL1B","CLEC10A"))

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

Myeloid_marker <- list('Mφ_CXCL8'=c('CXCL8','CCL23','CCL3','CCL20','CCL4'),
                       'Mono_FCN1'=c('IL1B','IL1RN','FCN1','DUSP6'),
                       'Mφ_APOE'=c("C1QA","C1QB",'TREM2','APOE','APOC1','GPNMB'),
                       'pDC'=c("GZMB","TSPAN13","LILRA4","ITM2C"),
                       'cDC1'=c("CLEC9A","CPNE3", "HLA-DRB1"),
                       'cDC2'=c('FCER1A',"CD1C",'CD1A','CD1E'),
                       'cDC3'=c("LAMP3","FSCN1","CD83","CSF2RA"))

##GSVA评分
exp_gsva <- gsva(as.matrix(exp), Myeloid_marker) #表达矩阵需转换为matrix格式
cellinfo <- t(exp_gsva) %>% data.frame() %>%
  mutate(time=sim_clinical$time,status=sim_clinical$status)
   #colnames(cellinfo)=str_replace(colnames(cellinfo),'CD4..','CD4_')
   #colnames(cellinfo)=str_replace(colnames(cellinfo),'CD8..','CD8_')
##细胞比例做批量cox回归分析
library(survival)
library(survminer)
library(dplyr)

lev<-rownames(exp_gsva)
#lev<-colnames(cellinfo)[1:12]
res<-list()
length(lev)
for (i in 1:length(lev)){
  #i=1
  gene_exp<-cellinfo[,lev[i]]
  group<-ifelse(gene_exp>median(gene_exp),"Group_1","Group_0")
  survival_dat <- data.frame(group=group,time=cellinfo$time,status=cellinfo$status,
                             stringsAsFactors = F)
  m=coxph(Surv(time, status) ~  group, data =  survival_dat)#age + stage+
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  #每次该细胞比例对应的cox回归结果
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  
  res[[i]]=tmp["groupGroup_1",]
}
res_data<-as.data.frame(res)
colnames(res_data)<-lev
res_data<-t(res_data)

##HR值汇总森林图
library(forestplot)
res_data<-res_data %>% data.frame() %>% arrange(coef)
res_data<-subset(res_data,!HRCIUL=='Inf') 

pdf('T_cell_forestplot.pdf', width = 4, height = 3)
forestplot::forestplot(
  labeltext = rownames(res_data),
  res_data[, c("HR", "HRCILL", "HRCIUL")],
  boxsize = 0.2,
  zero = 1,
  xlog = F,txt_gp=fpTxtGp(cex=1),xlab="HR (95%CI)",
  col = fpColors(lines = "#3C5488FF", box = "#DC0000FF")
)
dev.off()

##感兴趣细胞特征的生存曲线
gene_exp <- cellinfo$Mφ_CXCL8
group <- ifelse(gene_exp > median(gene_exp),"Mφ_CXCL8_high","Mφ_CXCL8_low")
survival_dat <- data.frame(group=group,time=cellinfo$time,status=cellinfo$status,
                           stringsAsFactors = F)
diff=survdiff(Surv(time, status)~group, data=survival_dat)
fit <- survfit(Surv(time, status)~group, data=survival_dat) #age + stage+

pdf('Mφ_CXCL8_KM_plot.pdf', width = 5, height = 4)
ggsurvplot(fit,
           data = survival_dat,
           conf.int = TRUE, 
           risk.table = F,
           pval=T,
           palette =c("#DC0000FF","#3C5488FF"))
dev.off()

