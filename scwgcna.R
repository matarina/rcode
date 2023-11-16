library(Seurat)
library(rliger)
library(scWGCNA)
filelist = list.files("/data/dk/ischematic_stroke/rawdata")
filelist

mcao1 = Read10X(file.path("/data/dk/ischematic_stroke/rawdata",filelist[1]))
mcao2 = Read10X(file.path("/data/dk/ischematic_stroke/rawdata",filelist[2]))
mcao3 = Read10X(file.path("/data/dk/ischematic_stroke/rawdata",filelist[3]))
sham1 = Read10X(file.path("/data/dk/ischematic_stroke/rawdata",filelist[4]))
sham2 = Read10X(file.path("/data/dk/ischematic_stroke/rawdata",filelist[5]))
sham3 = Read10X(file.path("/data/dk/ischematic_stroke/rawdata",filelist[6]))
s1 <- CreateSeuratObject(counts = sham1, project = "sham1" , min.cells = 3,min.features = 200)
s2 <- CreateSeuratObject(counts = sham2, project = "sham2" , min.cells = 3,min.features = 200)
s3 <- CreateSeuratObject(counts = sham3, project = "sham3" , min.cells = 3,min.features = 200)
m1 <- CreateSeuratObject(counts = mcao1, project = "mcao1" , min.cells = 3,min.features = 200)
m2 <- CreateSeuratObject(counts = mcao2, project = "mcao2" , min.cells = 3,min.features = 200)
m3 <- CreateSeuratObject(counts = mcao3, project = "mcao3" , min.cells = 3,min.features = 200)
class(s1)
#seurat_obj <- merge(s1,y = c(s2,s3,m1,m2,m3),add.cell.ids=c("sham1","sham2","sham3","mcao1","mcao2","mcao3"))
seurat_obj <- merge(m1,y = c(m2,m3),add.cell.ids=c("mcao1","mcao2","mcao3"))
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
expression_list <- list(
  'sham1' = GetAssayData(s1, slot='counts'),
  'sham2' = GetAssayData(s2, slot='counts'),
  'sham3' = GetAssayData(s3, slot='counts'),
  'mcao1' = GetAssayData(m1, slot='counts'),
  'mcao2' = GetAssayData(m2, slot='counts'),
  'mcao3' = GetAssayData(m3, slot='counts')
)
seurat_obj = NormalizeData(seurat_obj)
seurat_obj = ScaleData(seurat_obj)



seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))


MmLimbE155.pcells = calculate.pseudocells(s.cells = seurat_obj, # Single cells in Seurat object
                                          seeds=0.2, # Fraction of cells to use as seeds to aggregate pseudocells
                                          nn = 10, # Number of neighbors to aggregate
                                          reduction = "pca", # Reduction to use
                                          dims = 1:10) # The dimensions to use
MmLimbE155.scWGCNA = run.scWGCNA(p.cells = MmLimbE155.pcells, # Pseudocells (recommended), or Seurat single cells
                                 s.cells = seurat_obj, # single cells in Seurat format
                                 is.pseudocell = T, # We are using single cells twice this time
                                 features = rownames(seurat_obj))
scW.p.dendro(scWGCNA.data = MmLimbE155.scWGCNA)

# Plot the expression of all modules at once
scW.p.expression(s.cells = seurat_obj, # Single cells in Seurat format
                       scWGCNA.data = seurat_obj, # scWGCNA list dataset
                       modules = "all", # Which modules to plot?
                       reduction = "tsne", # Which reduction to plot?
                       ncol=3) # How many columns to use, in case we're plotting several?
