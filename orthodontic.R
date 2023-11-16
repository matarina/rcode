library(Seurat)
library(EnhancedVolcano)
library(ggplot2)
library(dplyr)
setwd('/data/dk/proj/orthodontic/figures')
ctrl <- read.csv("/data/dk/proj/orthodontic/GSM5640072_sample2_raw_Expression.csv", header = T)
rownames(ctrl) <- ctrl[, 2]
ctrl <- ctrl[, -c(1, 2)]
ctrl <- CreateSeuratObject(counts = ctrl, project = "ctrl", min.cells = 3, min.features = 200)
case <- read.csv("/data/dk/proj/orthodontic/GSM5640073_sample3_raw_Expression.csv", header = T)
rownames(case) <- case[, 2]
case <- case[, -c(1, 2)]
case <- CreateSeuratObject(counts = case, project = "case", min.cells = 3, min.features = 200)
obj <- merge(ctrl, case, add.cell.ids = c("ctrl", "case"))


obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(obj)
obj <- ScaleData(obj, features = all.genes)
obj <- RunPCA(obj, features = VariableFeatures(object = obj))
ElbowPlot(obj)
obj <- FindNeighbors(obj, dims = 1:10)
obj <- FindClusters(obj, resolution = 0.3)
obj <- RunUMAP(obj, dims = 1:10)
head(obj@meta.data)
png('group_umap.png',width =400, height =400)
DimPlot(object = obj, reduction = "umap", group.by = "orig.ident")
dev.off()
diffgene <- FindMarkers(obj, logfc.threshold = 0.01, ident.1 = "case",
                         ident.2 = "ctrl", group.by = "orig.ident")
head(diffgene)
png('volcano.png',width =400, height =400)
EnhancedVolcano(diffgene,
    lab = rownames(diffgene),
    pCutoff = 0.05,
    labSize = 0.0,
    FCcutoff = 0.7,
    x = 'avg_log2FC',
    y = 'p_val')
dev.off()
upgene <- as.data.frame(diffgene) %>% filter(avg_log2FC > 0.7 & p_val < 0.05)
downgene <- as.data.frame(diffgene) %>% filter(avg_log2FC < -0.7 & p_val < 0.05)
diff <- rbind(upgene, downgene)
nrow(upgene)
diff$gene <- rownames(diff)
cat(rownames(diff), sep = "\n")
head(diff)
write.table(diffgene, "diffgene.csv", quote = FALSE, row.names = TRUE, sep = "\t")


expressmatrix <- AverageExpression(obj, features = rownames(diffgene), group.by = "orig.ident")
expressmatrix
write.table(expressmatrix, "expressmatrix.csv", quote = FALSE, row.names = FALSE, sep = "\t")

#perform enrich analysis

library(clusterProfiler)
library(org.Mm.eg.db)
transid  <- bitr(rownames(diff),fromType = 'ENSEMBL',toType = 'ENTREZID',OrgDb = 'org.Mm.eg.db') 
diff = merge(transid,diff,by.x = 'ENSEMBL', by.y = 'row.names')
diff = diff %>% arrange(desc(avg_log2FC))
diff
geneList = diff$avg_log2FC
names(geneList) = as.character(diff$ENTREZID)
geneList
go <- gseGO(
  geneList     = geneList,
  OrgDb        = org.Mm.eg.db,
  ont          = "ALL",
  minGSSize    = 100,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)
kegg <- enrichKEGG(
  gene = diff$ENTREZID,
  organism = "mmu",
  pvalueCutoff = 0.05
)
dotplot(kegg)
kegg <- gseKEGG(
  geneList = geneList,
  organism     = 'mmu',
  minGSSize    = 120,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)
gseaplot(kegg2, 1)








#scWGCNA
library(scWGCNA)
MmLimbE155.pcells <- calculate.pseudocells(
  s.cells = obj, # Single cells in Seurat object
  seeds = 0.2, # Fraction of cells to use as seeds to aggregate pseudocells
  nn = 10, # Number of neighbors to aggregate
  reduction = "pca", # Reduction to use
  dims = 1:15
) # The dimensions to use
MmLimbE155.scWGCNA <- run.scWGCNA(
  p.cells = obj, # Pseudocells (recommended), or Seurat single cells
  s.cells = obj, # single cells in Seurat format
  is.pseudocell = F, # We are using single cells twice this time
  features = VariableFeatures(obj)
)
png("~/dk/proj/orthodontic/figures/wgcna.png", width = 500, height = 400)
scW.p.dendro(scWGCNA.data = MmLimbE155.scWGCNA)
dev.off()
scW.p.expression(
  s.cells = my.small_MmLimbE155, # Single cells in Seurat format
  scWGCNA.data = MmLimbE155.scWGCNA, # scWGCNA list dataset
  modules = 1, # Which modules to plot?
  reduction = "tsne", # Which reduction to plot?
  ncol = 3
)
# plot module network
inter <- c()
MmLimbE155.scWGCNA <- scWGCNA.networks(scWGCNA.data = MmLimbE155.scWGCNA)
for (i in names(MmLimbE155.scWGCNA$modules)) {
  interset <- length(intersect(rownames(MmLimbE155.scWGCNA$modules[[i]]), diff$ENSEMBL))
  inter <- append(inter, interset)
  print(paste('total number:',length(rownames(MmLimbE155.scWGCNA$modules[[i]]))))
  print(paste("intersect:",interset))
}

MmLimbE155.scWGCNA  <-  scWGCNA.networks(scWGCNA.data = MmLimbE155.scWGCNA)
scW.p.network(MmLimbE155.scWGCNA, module = 1)

tfs = read.csv('/data/dk/proj/orthodontic/Mus_musculus_TF.txt',sep= '\t')
head(tfs)
inflammation = c('ENSMUSG00000027398','ENSMUSG00000049130','ENSMUSG00000045551','ENSMUSG00000021025','ENSMUSG00000062825','ENSMUSG00000018930','ENSMUSG00000025746')

library(venneuler)
MyVenn <- venneuler(c(diffgene=51,Module1=63,
                       "diffgene&Module1"=16))
MyVenn$labels <- c("diffgene\n35","Module1\n63")
png('venn.png',width =400, height =400)
plot(MyVenn)
text(0.51,0.52,"16")
dev.off()

cat(intersect(rownames(MmLimbE155.scWGCNA$modules[[1]]), diff$ENSEMBL),sep = '\n')

a = intersect(rownames(MmLimbE155.scWGCNA$modules[[1]]), diff$ENSEMBL)
length(a)
diff[diff$ENSEMBL %in% a,]
ENSMUSG00000071561
ENSMUSG00000054905
ENSMUSG00000022902
ENSMUSG00000079597
ENSMUSG00000033508
ENSMUSG00000059657
ENSMUSG00000079594
ENSMUSG00000045551
ENSMUSG00000063234
ENSMUSG00000044103
ENSMUSG00000032369
ENSMUSG00000093973
ENSMUSG00000096719
ENSMUSG00000037095
ENSMUSG00000022126
ENSMUSG00000026073

nrow(diff %>% filter(avg_log2FC >0))
40
nrow(diff %>% filter(avg_log2FC <0))
11

pdf('fea.pdf')
FeaturePlot(obj,features =c('ENSMUSG00000044103','ENSMUSG00000026073'),split.by = 'orig.ident')
dev.off()
