library(Seurat)
a = Read10X_h5('/data/dk/proj/diabetes_dy/raw/GSM4453619_Control1_HPAP022_molecule_info.h5', use.names = TRUE, unique.features = TRUE)
library(hdf5r)
Read10X_h5('/data/dk/proj/diabetes_dy/raw/GSM4453620_Control2_HPAP034_molecule_info.h5')
a = h5file('/data/dk/proj/diabetes_dy/raw/GSM4453620_Control2_HPAP034_molecule_info.h5')
a = Read10X_h5('/data/dk/proj/diabetes_dy/raw/GSM4453621_Control3_HPAP035_molecule_info.h5')

library(scCustomize)
Read10X_h5_GEO('/data/dk/proj/diabetes_dy/raw/')

tmp = DropletUtils::read10xMolInfo('/data/dk/proj/diabetes_dy/raw/GSM4453621_Control3_HPAP035_molecule_info.h5')
mtx = DropletUtils::makeCountMatrix(tmp$data$gene, tmp$data$cell, value = tmp$data$reads)
dimnames(mtx) = list(tmp$genes[tmp$data$gene], tmp$cell)
mtx[1:3, 1:10]
