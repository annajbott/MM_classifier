## Build classifiers

library(Seurat)
library(tidyverse)
library(scClassify)

naive <- readRDS("/stopgap/cribbslab/proj012/RBC_removed/RDS_objects.dir/Harmony_2_covariates_scclassify_annotated_SeuratObject.rds")

relapse <- readRDS("/stopgap/cribbslab/proj012/relapsed_RBC_removed/RDS_objects.dir/Harmony_clustifyr_annotated_SeuratObject.rds")

# Add manual cell types naive
cell_type_ano2 <- tibble(seurat_clusters = 0:17)

cell_type_ano2$manual_annotation2 <- c("CD4_T_cell", "CD8_T_cell", "MM_cell", "Monocyte", "NK_cell",
                                    "Lymphoid_origin","B_cell", "MM_cell", "Lymphoid_origin",
                                     "Myeloid_origin","Plasma_cell","Myeloid_origin","Myeloid_origin",
                                     "MM_cell","Myeloid_origin",
                                    "Dendritic_cell","CD4_T_cell", "Dendritic_cell")

meta2 <- plyr::join(naive@meta.data, cell_type_ano2,
                                 by = "seurat_clusters", type = "left")

naive@meta.data$manual_annotation_general <- meta2$manual_annotation2

# Add manual cell types relapse
cell_type_ano <- tibble(seurat_clusters = 0:14)

cell_type_ano$manual_annotation <- c("CD4+_T_cell","CD4+_T_cell", "Myeloid_origin", "Myeloid_origin",
                                    "MM_cell", "NK_cell", "Myeloid_origin", "CD8_T_cell", "Monocyte",
                                    "Relapsed_B_cell", "T_cell", "CD4+_T_cell", "B_cell", "Myeloid_origin",
                                    "T_reg")


meta <- plyr::join(relapse@meta.data, cell_type_ano,
                                 by = "seurat_clusters", type = "left")

relapse@meta.data$manual_annotation <- meta$manual_annotatio


exprs_naive <- GetAssayData(object = naive, slot = "data") 
exprs_relapse <- GetAssayData(object = relapse, slot = "data") 
#
cellTypes_naive <- naive@meta.data$manual_annotation_general
cellTypes_relapse <- relapse@meta.data$manual_annotation

trainClass_naive <- train_scClassify(exprsMat_train = exprs_naive,
                               cellTypes_train = cellTypes_naive,
                               selectFeatures = c("limma", "BI"),
                               returnList = FALSE
)

trainClass_relapse <- train_scClassify(exprsMat_train = exprs_relapse,
                               cellTypes_train = cellTypes_relapse,
                               selectFeatures = c("limma", "BI"),
                               returnList = FALSE
)

saveRDS(object = trainClass_naive, file = "MM_naive_scclassify_model.rds" )

saveRDS(object = trainClass_relapse, file = "MM_relapsed_scclassify_model.rds")