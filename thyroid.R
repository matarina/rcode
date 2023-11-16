rm(list= ls())
library(survminer)
library(dplyr)
library(Seurat)
library(patchwork)
library(monocle3)
library(scWGCNA)
# read data in
setwd('/data/dk/proj/thyroid/figures/')
para1 = Read10X('~/dk/proj/thyroid/raw/para1/')
para2 = Read10X('~/dk/proj/thyroid/raw/para2/')
para3 = Read10X('~/dk/proj/thyroid/raw/para3/')
para5 = Read10X('~/dk/proj/thyroid/raw/para5/')
para8 = Read10X('~/dk/proj/thyroid/raw/para8/')
para9 = Read10X('~/dk/proj/thyroid/raw/para9/')
tumor1 = Read10X('~/dk/proj/thyroid/raw/tumor1/')
tumor2 = Read10X('~/dk/proj/thyroid/raw/tumor2/')
tumor3 = Read10X('~/dk/proj/thyroid/raw/tumor3/')
tumor5 = Read10X('~/dk/proj/thyroid/raw/tumor5/')
tumor8 = Read10X('~/dk/proj/thyroid/raw/tumor8/')
tumor9 = Read10X('~/dk/proj/thyroid/raw/tumor9/')

p1 <- CreateSeuratObject(counts = para1, project = "para1" , min.cells = 3,min.features = 200)
p2 <- CreateSeuratObject(counts = para2, project = "para2" , min.cells = 3,min.features = 200)
p3 <- CreateSeuratObject(counts = para3, project = "para3" , min.cells = 3,min.features = 200)
p5 <- CreateSeuratObject(counts = para5, project = "para5" , min.cells = 3,min.features = 200)
p8 <- CreateSeuratObject(counts = para8, project = "para8" , min.cells = 3,min.features = 200)
p9 <- CreateSeuratObject(counts = para9, project = "para9" , min.cells = 3,min.features = 200)
t1 <- CreateSeuratObject(counts = tumor1, project = "tumor1" , min.cells = 3,min.features = 200)
t2 <- CreateSeuratObject(counts = tumor2, project = "tumor2" , min.cells = 3,min.features = 200)
t3 <- CreateSeuratObject(counts = tumor3, project = "tumor3" , min.cells = 3,min.features = 200)
t5 <- CreateSeuratObject(counts = tumor5, project = "tumor5" , min.cells = 3,min.features = 200)
t8 <- CreateSeuratObject(counts = tumor8, project = "tumor8" , min.cells = 3,min.features = 200)
t9 <- CreateSeuratObject(counts = tumor9, project = "tumor9" , min.cells = 3,min.features = 200)

obj <- merge(p1,y = c(p2,p3,p5,p8,p9,t1,t2,t3,t5,t8,t9), add.cell.ids=c("p1","p2","p3","p5","p8","p9","t1","t2","t3","t5","t8","t9"))
rm(p1,p2,p3,p5,p8,p9,t1,t2,t3,t5,t8,t9)


#preprocess , tidy data 
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(obj)
obj <- ScaleData(obj, features = all.genes)
obj <- RunPCA(obj, features = VariableFeatures(object = obj))
ElbowPlot(obj)
obj <- FindNeighbors(obj, dims = 1:15)
obj <- FindClusters(obj, resolution = 0.3)
obj <- RunUMAP(obj, dims = 1:15)
#plot qc of features statistics
png("~/dk/proj/thyroid/figures/qc.png",width = 600,height = 400, pointsize = 1)
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
dev.off()

#add a new gorup, plot umap between group
groupcol = ifelse(grepl("para", obj$orig.ident), "para","tumor" )
obj <- AddMetaData(obj, metadata = groupcol, col.name = "group")
DimPlot(object = obj, reduction = 'umap', group.by = 'group')
obj.markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)




#find cell markers
table(obj.markers$cluster)
allmarker = obj.markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)
for (i in 0:22){
  cat("\n\n")
  print(i)
  cat(allmarker %>% filter(cluster==i) %>% pull(gene) %>% paste(collapse = ','))
}

obj = RenameIdents(object = obj,'2'='B cell')
obj = RenameIdents(object = obj,'4'='Myeloid cell')
obj = RenameIdents(object = obj,'6'='Thyrocytes cell')
obj = RenameIdents(object = obj,'7'='Thyrocytes cell')
obj = RenameIdents(object = obj,'8'='Thyrocytes cell')
obj = RenameIdents(object = obj,'9'='Filbroblast cell')
obj = RenameIdents(object = obj,'11'='Endothelial cell')
obj = RenameIdents(object = obj,'12'='Thyrocytes cell')
obj = RenameIdents(object = obj,'13'='B cell')
obj = RenameIdents(object = obj,'14'='Thyrocytes cell')
obj = RenameIdents(object = obj,'15'='Endothelial cell')
obj = RenameIdents(object = obj,'17'='Myeloid cell')
obj = RenameIdents(object = obj,'20'='Myeloid cell')
obj = RenameIdents(object = obj,'21'='Thyrocytes cell')
obj = RenameIdents(object = obj,'22'='Filbroblast cell')
obj = RenameIdents(object = obj,'0'='Others')
obj = RenameIdents(object = obj,'1'='Others')
obj = RenameIdents(object = obj,'3'='Others')
obj = RenameIdents(object = obj,'5'='Others')
obj = RenameIdents(object = obj,'10'='Others')
obj = RenameIdents(object = obj,'16'='Others')
obj = RenameIdents(object = obj,'18'='Others')
obj = RenameIdents(object = obj,'19'='Others')
#plot umap by group
png("~/dk/proj/thyroid/figures/umap_group.png",width = 500,height = 400, pointsize = 1)
DimPlot(object = obj, reduction = 'umap',group.by = 'group',label = TRUE)
dev.off()

DimPlot(object = obj, reduction = 'umap',label = TRUE)
levels(obj)
saveRDS(obj,"~/dk/proj/thyroid/obj.RDS")
thyroid = subset(obj, idents = 'Thyrocytes cell')
saveRDS(thyroid, "~/dk/proj/thyroid/thyobj.RDS")
thyroidmarker = FindMarkers(obj ,logfc.threshold = 1, ident.1 = "tumor", ident.2 = 'para', group.by = 'group', subset.ident = "Thyrocytes cell")
head(thyroidmarker)

thyroid = subset(obj, idents = 'Thyrocytes cell')
express = data.frame(AverageExpression(thyroid, group.by = 'group'))
express[1:3,1:2]
diffexpress = merge(express, thyroidmarker, by = "row.names")
head(diffexpress)

library('ggpubr')
colnames(diffexpress)[4]  = 'padj'
colnames(diffexpress)[5]  = 'log2FoldChange'

diffexpress$baseMean = (diffexpress$RNA.para+diffexpress$RNA.tumor)/2
png("~/dk/proj/thyroid/figures/maplot.png",width = 500,height = 400, pointsize = 1)
ggmaplot(diffexpress, size= 1.5, fdr = 0.05, fc =1,top = 20)
dev.off()
saveRDS(obj,'obj.RDS')
obj = readRDS('obj.RDS')




#WGCNA
astro <- merge(t1,y = c(t2,t3,t5,t8,t9),add.cell.ids=c("tumor1","tumor2","tumor3","tumor5","tumor8","tumor9"))
astro[["percent.mt"]] <- PercentageFeatureSet(astro, pattern = "^MT-")
astro <- NormalizeData(astro, normalization.method = "LogNormalize", scale.factor = 10000)
astro <- FindVariableFeatures(astro, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(astro)
astro <- ScaleData(astro, features = all.genes)
astro <- RunPCA(astro, features = VariableFeatures(object = astro))
ElbowPlot(astro)
astro <- FindNeighbors(astro, dims = 1:15)
astro <- FindClusters(astro, resolution = 0.3)
astro <- RunUMAP(astro, dims = 1:15)
astro = RunTSNE(object = astro,dims.use = 1:15)
DimPlot(astro, reduction = "tsne")
MmLimbE155.pcells = calculate.pseudocells(s.cells = astro, # Single cells in Seurat object
                                          seeds=0.2, # Fraction of cells to use as seeds to aggregate pseudocells
                                          nn = 10, # Number of neighbors to aggregate
                                          reduction = "pca", # Reduction to use
                                          dims = 1:15) # The dimensions to use
MmLimbE155.scWGCNA = run.scWGCNA(p.cells = astro, # Pseudocells (recommended), or Seurat single cells
                                 s.cells = astro, # single cells in Seurat format
                                 is.pseudocell = F, # We are using single cells twice this time
                                 features = VariableFeatures(astro))
setwd('~/dk/proj/thyroid/figures/')
saveRDS(MmLimbE155.scWGCNA,'scwgcna.RDS')
MmLimbE155.scWGCNA = readRDS('scWGCNA.RDS')

png("~/dk/proj/thyroid/figures/wgcna.png",width = 500,height = 400)
scW.p.dendro(scWGCNA.data = MmLimbE155.scWGCNA)
dev.off()
scW.p.expression(s.cells = my.small_MmLimbE155, # Single cells in Seurat format
                 scWGCNA.data = MmLimbE155.scWGCNA, # scWGCNA list dataset
                 modules = 11, # Which modules to plot?
                 reduction = "tsne", # Which reduction to plot?
                 ncol=3) 
MmLimbE155.scWGCNA = scWGCNA.networks(scWGCNA.data = MmLimbE155.scWGCNA)
scW.p.network(MmLimbE155.scWGCNA, module=11)
#get intersect genes for pathway analysis
length(diffexpress$Row.names)
length((rownames(MmLimbE155.scWGCNA$modules[[11]])))
inter <- c()
for (i in names(MmLimbE155.scWGCNA$modules)) {interset  = length(intersect(rownames(MmLimbE155.scWGCNA$modules[[i]]), diffexpress$Row.names));inter = append(inter,interset);print(interset)}
names(inter) = paste('module',seq(1:12))
intersect_bar
intersect_bar = data.frame(module = paste('module',seq(1:12)),internumber = inter)
intersect_bar$module = reorder(intersect_bar$module,intersect_bar$inter)
png("~/dk/proj/thyroid/figures/intersect_bar.png",width = 600,height = 400)
ggplot(intersect_bar,aes(x=module,y=internumber,fill = 'module'))+geom_bar(stat = 'identity')+theme_minimal()+theme(legend.position = "none")
dev.off()
#plot intersect venn plot
library(venneuler)
MyVenn <- venneuler(c(diffgene=122,Module11=397,
                       "diffgene&Module11"=87))
MyVenn$labels <- c("diffgene\n122","Module11\n397")
png('venn.png',width =400, height =400)
plot(MyVenn)
text(0.53,0.52,"87")
dev.off()


inge = intersect(rownames(MmLimbE155.scWGCNA$modules[[11]]),diffexpress$Row.names)
cat(inge)
intersect_matrix = diffexpress[diffexpress$Row.names %in% inge,]
write.table(intersect_matrix,file = '~/dk/proj/thyroid/network_matrx.txt',quote = FALSE,sep = '\t',row.names = FALSE)
cat(diffexpress$Row.names)


#enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)
write.csv(diffexpress,'/data/dk/proj/thyroid/diffexpress.csv',quote = FALSE,row.names = FALSE)
diffexpress = read.csv('/data/dk/proj/thyroid/diffexpress.csv')
head(diffexpress)
intersect_matrix = read.csv('/data/dk/proj/thyroid/analyst_network_input.txt',sep = ' ')
inge = intersect_matrix$Row.names

inge_diff = diffexpress[diffexpress$Row.names %in% inge,]
upgene = inge_diff[inge_diff$avg_log2FC > 0,]
downgene = inge_diff[inge_diff$avg_log2FC < 0,]
inge =inge_diff$Row.names 
gene <- bitr(inge,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Hs.eg.db')
go = enrichGO(gene = inge,
	      keyType = 'SYMBOL',
	      OrgDb = org.Hs.eg.db,
	      ont = 'BP',
	      pAdjustMethod = 'BH',
	      pvalueCutoff = 0.05,
	      qvalueCutoff = 0.05)

kegg = enrichKEGG(gene =  gene$ENTREZID,
		  organism = 'hsa',
		  pvalueCutoff = 0.05,
	          qvalueCutoff = 0.05)
png("~/dk/proj/thyroid/figures/go.png",width = 400,height = 600)
dotplot(go)
dev.off()
png("~/dk/proj/thyroid/figures/kegg.png",width = 400,height = 600)
dotplot(kegg)
dev.off()

head(gene)
hormone = c(2878,389434,5172,7038,7173)
hormonegene = gene[gene$ENTREZID %in% hormone,]
hormonegene
diffexpress[diffexpress$Row.names %in% hormonegene$SYMBOL,]
#trajacet analysis
library(SeuratWrappers)

monobj <- as.cell_data_set(obj)
mono_tumor <- monobj[, monobj$orig.ident %in% c("tumor1","tumor2","tumor3","tumor5","tumor8","tumor9")]

mono_tumor <- cluster_cells(cds = mono_tumor, reduction_method = "UMAP")
mono_tumor <- learn_graph(mono_tumor, use_partition = TRUE)
root_group = colnames(subset(mono_tumor, monobj$indent == "Thyrocytes cell"))
mono_tumor <- order_cells(mono_tumor, reduction_method = "UMAP",root_cells = root_group)
mono_tumor <- order_cells(mono_tumor, reduction_method = "UMAP")

VlnPlot(obj, features = c('FN1','LRRK2','TSC22D1','CCND1','PRDX1','CLU','APOE','DSP','LGALS3','EPS8','MET','ANXA1','CRYAB'))
                          
FeaturePlot(obj, features = lasso_gene,ncol = 3)
ciliated_genes <- c('FN1','LRRK2','TSC22D1','CCND1','PRDX1','CLU','APOE','DSP','LGALS3','EPS8','MET','ANXA1','CRYAB')
rowData(mono_tumor)$gene_short_name <-row.names(rowData(mono_tumor))
plot_cells(mono_tumor)
plot_cells(mono_tumor,genes=ciliated_genes,label_cell_groups=FALSE,show_trajectory_graph=TRUE,trajectory_graph_segment_size =0.5,graph_label_size = 0.1,label_branch_points = FALSE)
png("~/dk/proj/thyroid/figures/trajacet.png",width = 400,height = 400)
plot_cells(mono_tumor,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,alpha = 0.1,
           graph_label_size=0.1)
dev.off()





library(tidyverse)

#integrate TCGA data and perform lasso regression
matrix = read.csv('~/dk/proj/thyroid/TCGA.THCA.sampleMap%2FHiSeqV2',header=TRUE,sep = '\t',row.names = 1)
matrix[1:4,1:4]
head(colnames(matrix))
matrix = t(matrix)
rownames(matrix) = gsub('\\.','-',rownames(matrix))
gene = read.csv('~/dk/proj/thyroid/network_matrx.txt',sep='\t')
head(gene)
matrix = matrix[,colnames(matrix) %in% gene$Row.names]
meta = read.csv('~/dk/proj/thyroid/meta.txt',sep = '\t')
type = c()
for (i in 1:nrow(meta)){
	string = meta[i,1]
	print(string)
	last_two <- as.numeric(substr(string, nchar(string) - 1, nchar(string)))
	if (last_two > 9){
		print('Normal')
		type = append(type,'Normal')
	}else
		{print('Tumor')
		type = append(type,'Tumor')
		}
}
meta$type = type
class(meta)
meta = meta %>% select_('sample','type','OS','OS.time') %>% column_to_rownames(var = 'sample')
me = merge(meta,matrix,by = 'row.names')
me[1:4,1:9]
class(meta$OS)
me$OS.time = round(me$OS.time/12,digits = 0)
me = filter(me, !OS.time == 0)
rownames(me) = me$Row.names
me = me[,-1]
library(pROC)

me[1:4,1:4]
res = roc(type~GPX3+IYD+SLC26A4+TG+TPO,data=me,aur=TRUE,
         smooth=TRUE,# 是否平滑曲线
         levels=c('Normal','Tumor'),direction=">" #设置分组方向
         )

library(ggplot2)
p<- ggroc(res, legacy.axes = TRUE)+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype=4)+
  theme_bw() + # 设置背景
  ggtitle("Genes ROC Curve")+
  theme(plot.title = element_text(hjust = 0.5,size = 16),
        axis.text=element_text(size=12,colour = "black"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
p
png("~/dk/proj/thyroid/figures/roc.png",width = 500,height = 400)
p+annotate("text",x=0.85,y=0.055,label=paste("GPX3-AUC = ", round(res$GPX3$auc,3)))+
	annotate("text",x=0.85,y=0.085,label=paste("IYD-AUC = ", round(res$IYD$auc,3)))+
	annotate("text",x=0.85,y=0.115,label=paste("SLC26A4-AUC = ", round(res$SLC26A4$auc,3)))+
	annotate("text",x=0.85,y=0.145,label=paste("TG-AUC = ", round(res$TG$auc,3)))+
	annotate("text",x=0.85,y=0.175,label=paste("TPO-AUC = ", round(res$TPO$auc,3)))
dev.off()



data("ToothGrowth")
df <- ToothGrowth

g1 = ggboxplot(me, x = "type", y = "GPX3",color = 'type', add = 'jitter',width = 0.8)
g2 = ggboxplot(me, x = "type", y = "IYD",color = 'type', add = 'jitter',width = 0.8)
g3 = ggboxplot(me, x = "type", y = "SLC26A4",color = 'type', add = 'jitter',width = 0.8)
g4 = ggboxplot(me, x = "type", y = "TG",color = 'type', add = 'jitter',width = 0.8)
g5 = ggboxplot(me, x = "type", y = "TPO",color = 'type', add = 'jitter',width = 0.8)

png("~/dk/proj/thyroid/figures/box.png",width = 900,height = 400)
ggarrange(g1,g2,g3,g4,g5,ncol = 5,common.legend = TRUE)
dev.off()




library('survival')
library(glmnet)
x = as.matrix(me[,c(4:ncol(me))])
y = data.matrix(Surv(me$OS.time,me$OS)) 
lasso = glmnet(x,y,family = "cox",alpha = 1 ,nlanmbda = 1000)
png("~/dk/proj/thyroid/figures/lasso_lambda.png",width = 400,height = 400)
plot(lasso)
dev.off()
cvfit = cv.glmnet(x,y,family = 'cox',type.measure = 'deviance',nfolds = 10)
library(ROCR)
library(pROC)
preds <- predict(fit, newx = x, type = 'response')
perf <- performance(prediction(preds, y$status), 'tpr', 'fpr')
head(y)
plot(perf)
data("aSAH")
head(aSAH)
png("~/dk/proj/thyroid/figures/lasso_cv.png",width = 400,height = 400)
plot(cvfit, xvar = 'lambda', label = TRUE)
dev.off()

cvfit$lambda.min
lambda.1se = cvfit$lambda.1se 
lambda.1se
library(survminer)
fit = survfit(Surv(OS.time,OS) ~ OS, data = me)
png("~/dk/proj/thyroid/figures/survival.png",width = 400,height = 400)
ggsurvplot(fit)
dev.off()
coefficient = coef(cvfit,s = cvfit$lambda.min)
coefficient
index = which(as.numeric(coefficient) != 0 )
index
indexcoff  = as.numeric(coefficient)[index]
indexcoff
rownames(coefficient)[index]
lasso_gene = c("TPD52L1","APOE")
png("~/dk/proj/thyroid/figures/umap_feature.png",width = 400,height = 800)
FeaturePlot(obj, features =hormonegene$SYMBOL ,ncol = 2, split.by = "group")
dev.off()
png("~/dk/proj/thyroid/figures/violin_feature.png",width = 800,height = 400)
VlnPlot(obj, features = lasso_gene)
dev.off()


