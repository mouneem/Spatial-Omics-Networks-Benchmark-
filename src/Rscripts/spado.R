remotes::install_github("bm2-lab/SpaDo")
library(SpaDo)
### make sure the version of Seurat and SeuratObject is 4.0.4
library(Seurat)
packageVersion("Seurat")
packageVersion("SeuratObject")

### load data
library(SpaDo)
load("~/Downloads/./osmFISH_demo.RData")


### If you lack cell type labels, you can utilize the Seurat clustering results to obtain preliminary cell type labels. Set the parameter "user_offered" to False by using "user_offered=F".
# test_expression_normalized<-SpatialNormalize(expression_profile = test_expression,ST_method = "osmFISH")
# initial_result<-InitialClustering(expression_profile = test_expression_normalized,user_offered = F)

### If you have your own cell type labels, you can directly execute the following code using the "user_offered=T" parameter.
initial_result<-InitialClustering(expression_profile = test_expression,user_offered = T,sample_information_user_offered = sample_information_cellType)
