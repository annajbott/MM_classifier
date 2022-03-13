## Build classifiers relapse MM dataset

library(Seurat)
library(optparse)
library(tidyverse)
library(scClassify)

option_list <- list(
    make_option(c("-i", "--input"), default=NULL),
    make_option(c("-o", "--output"), default=NULL)
)

# Read in options
opt <- parse_args(OptionParser(option_list=option_list))

input_file <- opt$input
output_file <- opt$output

seurat_object <- readRDS(input_file)
    
dgc_mat <- GetAssayData(object = seurat_object, slot = "data") 
    
relapse_pretrain <- readRDS("MM_relapsed_scclassify_model.rds")

# list_train <- list(naive_pretrain, relapse_pretrain)


similarity  <- "pearson" # Error with scClassify where multiple metrics seem to make predict_scClassifyJoint fail
pred_res <- predict_scClassify(exprsMat_test = dgc_mat,
                             trainRes = relapse_pretrain, # or pre-made reference using
                             cellTypes_test = NULL,
                             algorithm = "WKNN",
                             features = c("limma"),
                             similarity = similarity,
                             prob_threshold = 0.7,
                             verbose = FALSE)

column_name <- paste0(similarity, "_WKNN_limma")
cells_assigned <- as.vector(pred_res[[column_name]]$predRes)

# pred_results  <- pred_res_naive$jointRes
# cells_assigned <- pred_results$cellTypes
    

seurat_object@meta.data[['MM_relapse_labels']] <- cells_assigned
saveRDS(seurat_object, output_file) # Save seurat object with labels