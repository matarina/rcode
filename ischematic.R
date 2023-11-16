rm(list= ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(monocle3)
library(scWGCNA)
setwd('/data/dk/proj/ischematic_stroke/')
filelist = list.files("/data/dk/ischematic_stroke/rawdata")
filelist


mcao1 = Read10X(file.path("/data/dk/ischematic_stroke/rawdata",filelist[1]))
mcao2 = Read10X(file.path("/data/dk/ischematic_stroke/rawdata",filelist[2]))
mcao3 = Read10X(file.path("/data/dk/ischematic_stroke/rawdata",filelist[3]))
sham1 = Read10X(file.path("/data/dk/ischematic_stroke/rawdata",filelist[4]))
sham2 = Read10X(file.path("/data/dk/ischematic_stroke/rawdata",filelist[5]))
sham3 = Read10X(file.path("/data/dk/ischematic_stroke/rawdata",filelist[6]))

m1 <- CreateSeuratObject(counts = mcao1, project = "mcao1" , min.cells = 3,min.features = 200)
m2 <- CreateSeuratObject(counts = mcao2, project = "mcao2" , min.cells = 3,min.features = 200)
m3 <- CreateSeuratObject(counts = mcao3, project = "mcao3" , min.cells = 3,min.features = 200)
s1 <- CreateSeuratObject(counts = sham1, project = "sham1" , min.cells = 3,min.features = 200)
s2 <- CreateSeuratObject(counts = sham2, project = "sham2" , min.cells = 3,min.features = 200)
s3 <- CreateSeuratObject(counts = sham3, project = "sham3" , min.cells = 3,min.features = 200)

obj <- merge(s1,y = c(s2,s3,m1,m2,m3),add.cell.ids=c("sham1","sham2","sham3","mcao1","mcao2","mcao3"))


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

groupcol = ifelse(grepl("sham", obj$orig.ident), "sham","mcao" )
obj <- AddMetaData(obj, metadata = groupcol, col.name = "group")

DimPlot(object = obj, reduction = 'umap', group.by = 'group')


obj.markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
table(obj.markers$cluster)
allmarker = obj.markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)
for (i in 0:17){
  cat("\n\n")
  cat(allmarker %>% filter(cluster==i) %>% pull(gene) %>% paste(collapse = ','))
}
all = data.frame(obj.markers)
schwann =  all %>% filter(cluster == 'Schwann cell') %>% pull(gene)












new.cluster.ids <- c('Endothelial cell','Choroid cell','CD14+ monocyte','mricroglia cell',
                     'Vascular cell','Choroid plexus cell','Quiescent neural stem cell','Vascular cell',
                     'Antigen-presenting cell','macrophage','Schwann cell','Pericyte-like cell','cDC2b',
                     'Immature astrocyte','microglia cell','Lymphocyte','Pericyte-like cell',
                     'Quiescent neural stem cell')

names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)
DimPlot(obj, reduction = "umap", label = TRUE,repel = TRUE, pt.size = 0.5) + NoLegend()


saveRDS(obj,'obj.RDS')
obj = readRDS('obj.RDS')

for (i in seq_along(new.cluster.ids)) {
  diff <- FindMarkers(obj, logfc.threshold = 1, ident.1 = "mcao", ident.2 = "sham", group.by = "group", subset.ident = new.cluster.ids[i])
  print(paste(paste('cell:',new.cluster.ids[i]),nrow(diff)))
}

apc  <- FindMarkers(obj, logfc.threshold = 1, ident.1 = "mcao", ident.2 = "sham", group.by = "group", subset.ident = 'Antigen-presenting cell')
rownames(apc)


#WGCNA
astro <- merge(s1,y = c(s2,s3),add.cell.ids=c("mcao1","mcao2","mcao3"))
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
MmLimbE155.scWGCNA = run.scWGCNA(p.cells = MmLimbE155.pcells, # Pseudocells (recommended), or Seurat single cells
                                 s.cells = astro, # single cells in Seurat format
                                 is.pseudocell = T, # We are using single cells twice this time
                                 features = VariableFeatures(astro))
scW.p.dendro(scWGCNA.data = MmLimbE155.scWGCNA)
scW.p.expression(s.cells = my.small_MmLimbE155, # Single cells in Seurat format
                 scWGCNA.data = MmLimbE155.scWGCNA, # scWGCNA list dataset
                 modules = 3, # Which modules to plot?
                 reduction = "tsne", # Which reduction to plot?
                 ncol=3) 
MmLimbE155.scWGCNA = scWGCNA.networks(scWGCNA.data = MmLimbE155.scWGCNA)
scW.p.network(MmLimbE155.scWGCNA, module=2)
#get intersect genes for pathway analysis


inter <- c()
for (i in names(MmLimbE155.scWGCNA$modules)) {
interset  = length(intersect(rownames(MmLimbE155.scWGCNA$modules[[i]]), rownames(apc)))
inter = append(inter,interset)
}
inter

names(inter) = paste('module',seq(1:21))

png('bar.png',width =400, height =400)
barplot(inter, main = " ", xlab = "Antigen-representing cell", ylab = "Intersected gene numbers")
dev.off()
intersectname = intersect(rownames(MmLimbE155.scWGCNA$modules[[13]]),rownames(apc))
cat(intersectname,sep='\n')
png('venn.png',width =400, height =400)
library(ggVennDiagram)
library(ggplot2)
x = list( diffgene = rownames(apc), Module13 =rownames(MmLimbE155.scWGCNA$modules[[13]])) 
ggVennDiagram(x) + scale_fill_gradient(low="#8EFFDA",high = "#F6AE33")
dev.off()



#enrichment analysis 

library(clusterProfiler)
library(org.Mm.eg.db)
gene = rownames(MmLimbE155.scWGCNA$modules[[13]])
geneid <- bitr(intersectname,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Mm.eg.db')
go = enrichGO(gene = gene,
	      keyType = 'SYMBOL',
	      OrgDb = org.Mm.eg.db,
	      ont = 'BP',
	      pAdjustMethod = 'BH',
	      pvalueCutoff = 0.05,
	      qvalueCutoff = 0.05)

kegg = enrichKEGG(gene = geneid$ENTREZID,
		  organism = 'mmu',
		  pvalueCutoff = 0.05,
	          qvalueCutoff = 0.05)
cat(intersectname,sep = '\n')
png("~/dk/proj/ischematic_stroke/figure/mod13go.png",width = 400,height = 600)
dotplot(go)
dev.off()
png("~/dk/proj/ischematic_stroke/figure/mod13kegg.png",width = 400,height = 600)
dotplot(kegg)
dev.off()
head(kegg)
apc[rownames(apc) %in% fea,]
fea = geneid[geneid$ENTREZID %in% c('20308','12768','13051'),][,1]
png("~/dk/proj/ischematic_stroke/figure/fea.png",width = 800,height = 1200,res = 150)
FeaturePlot(obj, features =fea,ncol = 3, split.by = "group")
dev.off()

#external validation
ev = read.csv('GSE199066_All.counts.exp.txt',sep = '\t')
ev = ev[,c(1,2,3,4,5,6,7)]
rownames(ev) = ev[,1]
ev = ev[,-1]
ev = as.data.frame(edgeR::cpm(ev))
ev = as.data.frame(t(ev))
ev[1:3,1:3]
ev$group = rep(c('sham','MCAO'),each =3,)
library(ggpubr)

g1 = ggboxplot(ev, x = "group", y = "Ccl9",color = 'group', add = 'jitter',width = 0.8)
g2 = ggboxplot(ev, x = "group", y = "Ccr1",color = 'group', add = 'jitter',width = 0.8)
g3 = ggboxplot(ev, x = "group", y = "Cx3cr1",color = 'group', add = 'jitter',width = 0.8)

png("~/dk/proj/ischematic_stroke/figure/box.png",width = 900,height = 400)
ggarrange(g1,g2,g3,ncol = 3,common.legend = TRUE)
dev.off()








#trajacet analysis
library(SeuratWrappers)

monobj <- as.cell_data_set(obj)
mono_mcao <- monobj[, monobj$orig.ident %in% c("mcao1","mcao2","mcao3")]

mono_mcao <- cluster_cells(cds = mono_mcao, reduction_method = "UMAP")
mono_mcao <- learn_graph(mono_mcao, use_partition = TRUE)
root_group = colnames(subset(mono_mcao, monobj$indent == "Antigen representing cell"))
mono_mcao <- order_cells(mono_mcao, reduction_method = "UMAP",root_cells = root_group)
mono_mcao <- order_cells(mono_mcao, reduction_method = "UMAP")


VlnPlot(obj, features = c("Hsp90aa1", "Pak1",'Prkcz','Gnai1','Mapt','Ap2a1'))
                          
FeaturePlot(obj, features = c("Hsp90aa1", "Pak1",'Prkcz','Gnai1','Mapt','Ap2a1'),ncol = 3)
ciliated_genes <- c("Hsp90aa1", "Pak1",'Prkcz','Gnai1','Mapt','Ap2a1')
rowData(mono_mcao)$gene_short_name <-row.names(rowData(mono_mcao))
plot_cells(mono_mcao)
plot_cells(mono_mcao,
           genes=ciliated_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=TRUE,
           trajectory_graph_segment_size =0.5,
           graph_label_size = 0.1,
           label_branch_points = FALSE)

plot_cells(mono_mcao,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,alpha = 0.1,
           graph_label_size=0.1)





nrow(apc)






###intersect gene with



