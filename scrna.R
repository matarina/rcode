rm(list= ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(monocle3)


count = readRDS('GSE182256_Export_counts.rds')
meta = read.csv('GSE182256_Export_Metadata.txt',sep = '\t',row.names = 1)
meta$condition = gsub('\\d','',meta$orig.ident)
obj <- CreateSeuratObject(counts = count,meta.data = meta, project = "renal_fibrosis", min.cells = 3, min.features = 200)

obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(obj)
obj <- ScaleData(obj, features = all.genes)

obj <- RunPCA(obj, features = VariableFeatures(object = obj))
obj <- FindNeighbors(obj, dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5)
head(Idents(obj), 5)
obj <- RunUMAP(obj, dims = 1:15)
DimPlot(obj, reduction = "umap",label=TRUE,split.by = 'condition')
DimPlot(obj, group.by = 'condition',reduction = "umap")



obj.markers <- FindAllMarkers(obj, only.pos = FALSE, min.pct = 0.25)
# a = obj.markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 10, order_by = avg_log2FC) %>% filter(cluster == '19') %>% select(gene)
# cat(a$gene,sep = ",")


table(obj.markers$cluster)

DimPlot(obj, reduction = "umap",label=TRUE,group.by = 'Cluster',  repel = TRUE,split.by = 'condition')+ NoLegend()

DimPlot(obj, reduction = "umap", label = TRUE,group.by = 'Cluster',  repel = TRUE,pt.size = 0.5) + NoLegend()
diff = obj.markers %>% filter(cluster == '16') 
saveRDS(obj, file = "obj.RDS")









library(SeuratWrappers)
cds  = as.cell_data_set(obj)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- cluster_cells(cds, resolution=1e-5)

cds <- reduce_dimension(cds)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "Cluster")

cds_subset <- choose_cells(cds)
cds_subset <- cluster_cells(cds_subset, resolution=1e-5)
cds_subset <- learn_graph(cds_subset)
plot_cells(cds_subset,
           color_cells_by = "Cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           cell_size = 0.55,
           group_label_size = 3,
           label_branch_points=FALSE)


##wgcna co expression analysis 
diff2 = obj.markers %>% filter(cluster == '16' ) %>% filter(avg_log2FC < -1 | avg_log2FC > 1) %>% filter(p_val < 0.05)
subobj = subset(x = obj, subset = (condition == "UUO"|Cluster == "CD PC"))
subobj <- subset(subobj, features = diff2$gene)


library(scWGCNA)
# Calculate the pseudocells
subobj <- JackStraw(subobj, num.replicate = 100)
subobj <- ScoreJackStraw(subobj, dims = 1:20)
subobj <- RunUMAP(subobj, dims = 1:15)
subobj <- RunTSNE(subobj, dims = 1:15)


# pseudocell = calculate.pseudocells(s.cells = subobj, # Single cells in Seurat object
#                                           seeds=0.2, # Fraction of cells to use as seeds to aggregate pseudocells
#                                           nn = 10, # Number of neighbors to aggregate
#                                           reduction = "umap" ,# Reduction to use
#                                           dims = 1:15) # The dimensions to use
scWGCNA = run.scWGCNA(p.cells = subobj, # Pseudocells (recommended), or Seurat single cells
                                 s.cells = subobj, # single cells in Seurat format
                                 is.pseudocell = F, # We are using single cells twice this time
                                 features = rownames(subobj)) # Recommended: variable genes
scW.p.dendro(scWGCNA.data = scWGCNA)

scW.p.expression(s.cells = subobj, # Single cells in Seurat format
                 scWGCNA.data = scWGCNA, # scWGCNA list dataset
                 modules = "all", # Which modules to plot?
                 reduction = "umap", # Which reduction to plot?
                 ncol=2) # How m

m3 = scWGCNA$modules$'3_turquoise'



library(clusterProfiler)
library(org.Mm.eg.db)
# Perform GO enrichment analysis using the enrichGO function
ego <- enrichGO(gene = diff2$gene,
                OrgDb = org.Mm.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2)
dotplot(ego)

# Perform GO enrichment analysis using the enrichGO function
ego <- enrichGO(gene = rownames(m3),
                OrgDb = org.Mm.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2)
dotplot(ego)



##Volcano plot
subobj2 = subset(x = obj, subset = Cluster == "CD PC")
subcount2 =as.data.frame(subobj2@assays[["RNA"]]@counts)
for (colname in colnames(subcount2)){
 if (meta[colname,6] == "UUO"){
   colnames(subcount2)[colnames(subcount2) == colname] = "UUO"
 } else{
   colnames(subcount2)[colnames(subcount2) == colname] = "Control"
 }
}
head(colname(subcount2))
target = intersect(hsaGenes,rownames(m4))
DoHeatmap(
  obj,
  features = target,
  cells = NULL,
  group.by = "Cluster",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = TRUE,
  size = 4.5,
  hjust = 0,
  angle = 90,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)
library(EnhancedVolcano)

EnhancedVolcano(diff,
                lab = rownames(diff),
                x = 'avg_log2FC',
                y = 'p_val',
                title = 'N061011 versus N61311',
                pCutoff = 10e-32,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0)


df <- as.data.frame(as.matrix(GetAssayData(subobj, assay="RNA", slot="scale.data")))
df = df[target,]
tdf = t(df)
write.csv(tdf,'tdf.csv',quote = FALSE,row.names = FALSE)

cat(target)


library(nichenetr)
frg = read.csv('FRG_Articles.csv')
hsaGenes <- convert_human_to_mouse_symbols(unique(frg$GeneName))
intersect(hsaGenes,rownames(m4))
