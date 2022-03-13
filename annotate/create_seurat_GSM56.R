library(Seurat)
library(tidyverse)

files_GSM56 <- list.files(pattern = "^GSM56")

barcodes_file <- files_GSM56[1]
features_file <- files_GSM56[2]
mat_file <- files_GSM56[3]

barcode <- read_tsv(barcodes_file, col_names = "barcodes")
feature <- read_tsv(features_file, col_names = FALSE)
colnames(feature) <- c("ENSEMBL", "SYMBOL", "type")

matrix2 <- Matrix::readMM(mat_file)
rownames(matrix2) <- feature$ENSEMBL
colnames(matrix2) <- barcode$barcodes

so_gsm56 <- CreateSeuratObject(counts = matrix2, min.cells = 20, min.features = 200) %>% 
    NormalizeData(verbose=FALSE) %>% 
    ScaleData(verbose=FALSE) %>% 
    FindVariableFeatures(verbose=FALSE)

saveRDS(object = so_gsm56, file = "GSM5687372_filtered_SeuratObject.rds") 