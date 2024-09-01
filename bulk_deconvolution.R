rm(list = ls())

#BiocManager::install("TCGAbiolinks")
#BiocManager::install("DESeq2")
#BiocManager::install("DEFormats")
#BiocManager::install("tximport")

library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
library(tximport)
library(readr)
library(dplyr)
library(stringr)
library(data.table)
library(biomaRt)
library(tibble)


####一、TCGA数据准备####
####1.TCGA数据读取####
load(file="TCGA_SKCM.rda")

####2.DESeq2标准化数据####
library(DESeq2)
dds <- DESeqDataSetFromMatrix(counts,
                              colData = sim_clinical,
                              design = ~ gender)
dds <- DESeq(dds) #数据标准化处理
#获取标准化的counts矩阵
counts <- data.frame(counts(dds, normalized=TRUE))
counts[1:5,1:5]


####2.DEFormats准备反卷积分析格式####
library(DEFormats)
dge <- as.DGEList(dds)
rownames(dge$samples)
colnames(counts) <- rownames(dge$samples)
pheno <- as(dge$samples, "AnnotatedDataFrame") # coerces a data.frame to an AnnotatedDataFrame.
eset <- ExpressionSet(assayData = as.matrix(counts), phenoData = pheno)
#保存eset文件
saveRDS(eset, file = "eset_tcga_SKCM.rds")
eset <- readRDS('eset_tcga_SKCM.rds')


####二、反卷积分析####
library(BisqueRNA) # for Seurat to ExpressionSet
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape)
library(Biobase)
library(survival)

####1.读取并处理单细胞数据####
#读取文件
setwd("~/科研/黑素瘤与白癜风/data")
scRNA<-readRDS("scRNA_merge_annotate_new.rds")
scRNA_T<-readRDS("scRNA_T&NK_new2.rds")
scRNA_Mono<-readRDS("scRNA_Mono.rds")
scRNA_DC<-readRDS("scRNA_DC_new.rds")
scRNA_Mela<-readRDS("scRNA_Mela.rds")
scRNA_Fibro<-readRDS("scRNA_Fibro_new2.rds")

sccell<-readRDS("scRNA_Fibro_new2.rds")
scRNA<-scRNA[,! scRNA$celltype %in% c('Melanocyte',"T & NK",'Langerhans cell','Mononulear phagocyte',"Fibroblast")]

scRNA$subcelltype<-scRNA$celltype
scRNA_T$subcelltype<-scRNA_T$T_celltype
scRNA_Mono$subcelltype<-scRNA_Mono$Mono_celltype
scRNA_DC$subcelltype<-scRNA_DC$Mono_celltype
scRNA_Fibro$subcelltype<-scRNA_Fibro$Fibro_celltype
scRNA_Mela$subcelltype<-scRNA_Mela$Mela_celltype


sccell<-merge(scRNA,y=c(scRNA_T,scRNA_Mono,scRNA_DC,scRNA_Fibro,scRNA_Mela))
sccell@meta.data[1:5,]
#随机抽取5000个细胞（仅用于演示）
#sccell<-sccell[,sample(nrow(sccell),5000)] 
rm(scRNA,scRNA_T,scRNA_Mono,scRNA_DC,scRNA_Fibro,scRNA_Mela)
gc()
#将sccell行名进行更改简化
scrna_eset <- SeuratToExpressionSet(sccell, delimiter="_", position=1,version="v3")
scrna_eset@phenoData$celltype <- Idents(sccell)
levels(Idents(sccell))


#建立单个细胞类型表达谱的参考基矩阵
####2.估计细胞比例####
res <- BisqueRNA::ReferenceBasedDecomposition(eset, scrna_eset, markers=NULL, use.overlap=F)
bulk<-res$bulk.props
table(rownames(sim_clinical)==colnames(bulk))
cellinfo<-cbind(sim_clinical,t(bulk))  #细胞比例
write.csv(cellinfo,file="TCGA_SKCM_deconvolution.csv")


####3.成纤维亚群与T细胞相关性####
library(vegan)
library(dplyr)
library(ggcor)
library(ggplot2)
library(corrplot) 

#提取成纤维细胞
df <- cellinfo[,16:18]
df = exp(df) #取指数，使得数据可以分析
#提取T细胞
env <- cellinfo[,19:27]
env = exp(env)

df_mantel = fortify_mantel(df, env, 
                           spec.select = list(irFib = 1, #依次定义四种物种作为Mantel的分析对象
                                              mrFib = 2,
                                              emFib = 3)) %>% 
  mutate(r = cut(r, breaks = c(-Inf, 0.25, 0.5, Inf), 
                 labels = c("<0.25", "0.25-0.5", ">=0.5"), #定义Mantel的R值范围标签，便于出图
                 right = FALSE),
         p.value = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),  
                       labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">=0.05"), #定义Mantel检验的p值范围标签，便于出图
                       right = FALSE))

quickcor(env,method = "spearman", type = "upper", cor.test = T, cluster.type = "all") +
  geom_square() +#相关性显示形式
  geom_mark(r = NA,sig.thres = 0.05, size = 3.5, colour = "black")+#显著性标签
  scale_fill_gradient2( high = '#DC0000FF', mid = 'white',low = '#4DBBD5FF') + #颜色设置
  add_link(df_mantel, mapping = aes(colour = p.value, size = r),
           diag.label = TRUE)+
  #scale_color_manual(values = c("#00A087FF","#3C5488FF","#E64B35FF"))+#线条颜色设置
  scale_size_manual(values = c(0.5, 1, 1.5, 2))+#连线粗细设置
  guides(fill = guide_colorbar(title = "correlation", order = 1),#图例相关设置
         size = guide_legend(title = "Mantel's r",order = 2),
         color = guide_legend(title = "Mantel's p", order = 3),
         linetype = "none")


####三、生存分析####
library(survival)
library(survminer)
library(dplyr)


##细胞比例做批量cox回归分析
lev<-levels(Idents(sccell))
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
res_data<-res_data%>%data.frame()%>%arrange(coef)
res_data

forestplot::forestplot(
  labeltext = rownames(res_data),
  res_data[, c("HR", "HRCILL", "HRCIUL")],
  boxsize = 0.2,
  zero = 1,
  xlog = F,txt_gp=fpTxtGp(cex=1),xlab="HR (95%CI)",
  col = fpColors(lines = "#3C5488FF", box = "#DC0000FF")
)
dev.off()

##感兴趣细胞比例的生存曲线
gene_exp<-cellinfo$mrFib
group<-ifelse(gene_exp>median(gene_exp),"mrFib_high","mrFib_low")
survival_dat <- data.frame(group=group,time=cellinfo$time,status=cellinfo$status,
                           stringsAsFactors = F)
diff=survdiff(Surv(time, status)~group, data=survival_dat)
fit <- survfit(Surv(time, status)~group, data=survival_dat) #age + stage+
ggsurvplot(fit,
           data = survival_dat,
           conf.int = TRUE, 
           palette =c("#DC0000FF","#3C5488FF"),
           risk.table = T,
           pval=T)


##分析成纤维亚群跟其他细胞亚群比例的相关性
cell_prop = cellinfo[,15:33]
cordata = cor(cell_prop)
corrplot(cordata, method = 'square', order = 'FPC', 
         type = 'lower', diag = FALSE)

colnames(cellinfo)
Fibro_prop = cellinfo[,15:17]
othercell_prop = cellinfo[,18:33]

ciber <- cibe1
exp <- as.data.frame(t(exp1_machine[hub_genes,]))
identical(rownames(Fibro_prop),rownames(othercell_prop))

cor<-sapply(othercell_prop,function(x,y) cor(x,y,method="spearman"), Fibro_prop) 
rownames(cor) <- colnames(Fibro_prop)
cor_res <- cor.mtest(cor,#计算p值
                     conf.level = 0.95)#置信区间

pdf("Fib_corplot.pdf",height=4,width=7)
corrplot(cor,  #尺寸10*8
         method = "color",#相关性矩阵展示的图形
         col=colorRampPalette(c("#01468b","white","#ee0000"))(100),
         addCoef.col = "black",#为相关系数添加颜色
         tl.col="black",#设置文本标签的颜色
         number.cex = 0.5,
         tl.cex = 0.7,
         cl.align = "l")
dev.off()



####四、细胞亚群的药物敏感性####
##https://cloud.tencent.com/developer/article/2353766?areaId=106001
library(pRRophetic)
library(ggplot2)
load("TCGA_SKCM.rda")
exp=as.matrix(counts)
data(cgp2016ExprRma)
data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)
allDrugs=unique(drugData2016$Drug.name)
allDrugs
# 264个样品的表达量矩阵
table(studyResponse)

###预测一个药物
predictedPtype <- pRRopheticPredict(testMatrix=as.matrix(counts), 
                                    drug="Tamoxifen", #Cisplatin
                                    tissueType = "all", 
                                    batchCorrect = "eb",
                                    selection=1,
                                    dataset = "cgp2016")
#boxplot(predictedPtype)
#添加分组信息
df <- data.frame(values = predictedPtype, 
                 group = ifelse(cellinfo$emFib>median(cellinfo$emFib), "High", "Low"))
ggplot(data = df,
       aes(y = values,
           x = group))+
  geom_boxplot(alpha = 0.3,
               fill = c('#e94753','#47a4e9'))+
  theme_bw()+
  ylab('Predicted Bortezomib Sensitivity') +
  xlab('Clinical Response') +
  stat_compare_means( method = "t.test")


###预测全部药物
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsci)

group=ifelse(cellinfo$emFib>median(cellinfo$emFib), "High", "Low")
group=group %>% data.frame()
colnames(group)='group'

##网站在线分析：http://www.sxdyc.com/singleCollectionTool

res.df <- read.table('predicted.drugs.score.txt')

allres_ttest = NULL
for (i in 1:ncol(res.df)){
  res_ttest = t.test(subset(plot_df,group=='High')[,i],
                     subset(plot_df,group=='Low')[,i],
                     alternative="two.sided")
  allres_ttest = c(allres_ttest, res_ttest$p.value)
}

plot_df <- res.df[, which(allres_ttest<5e-10)] %>% 
  bind_cols(group=df$group) #分组信息

plot_df = melt(plot_df,id.vars=c("group"))
colnames(plot_df) = c("Group","Drug", "IC50")
plot_df = na.omit(plot_df)

pdf(file=paste0('emFib_IC50_boxplot.pdf'),width = 10,height = 6)
ggboxplot(plot_df, x="Group", y="IC50", width = 0.6, #按group分组
          color = "black",#轮廓颜色
          fill="Group",#填充
          palette = cors,
          xlab = F, #不显示x轴的标签
          bxp.errorbar=T,#显示误差条
          bxp.errorbar.width=0.5, #误差条大小
          size=0.5, #箱型图边线的粗细
          outlier.shape=NA, #不显示outlier
          legend = "right") + 
  ylim(0,20) +
  facet_wrap(~Drug,scales = "free_y",nrow = 2) + # 按照细胞类型分面
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  stat_compare_means(size = 3) # 添加t检验，修改p值大小
dev.off()
