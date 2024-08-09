print('R script Loaded')

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(SpatialExperiment)
  library(Giotto)
})

args <- commandArgs(trailingOnly = TRUE)

coord_file <- args[1]
matrix_file <- args[2]
feature_file <- args[3]
observation_file <- args[4]
out_dir <- args[5]
config_file <- args[6]
n_clusters <- args[7]
technology <- args[8]
seed <- args[9]

# Output files
label_file <- file.path(out_dir, "domains.tsv")
embedding_file <- file.path(out_dir, "embedding.tsv")
# if additional output files are required write it also to out_dir

config <- fromJSON(config_file)

# You can get SpatialExperiment directly
get_SpatialExperiment <- function(
    feature_file,
    observation_file,
    coord_file,
    matrix_file = NA,
    assay_name = "counts",
    reducedDim_name = "reducedDim") {
  rowData <- read.csv(feature_file, stringsAsFactors = FALSE, row.names = 1)
  colData <- read.csv(observation_file, stringsAsFactors = FALSE, row.names = 1)
  coordinates <- read.csv(coord_file, row.names = 1)
  coordinates <- as.matrix(coordinates[rownames(colData), ])
  coordinates[,c(1:2)] <- as.numeric(coordinates[,c(1:2)])
  
  colData <- cbind(coordinates, colData)
  spe <- SpatialExperiment::SpatialExperiment(
    rowData = rowData, colData = colData, spatialCoords = coordinates
  )
  
  exp <-  t(read.csv(matrix_file, row.names = 1))
  assay(spe, assay_name, withDimnames = FALSE) <- as(as.matrix(exp), "CsparseMatrix")
  assay(spe, "logcounts", withDimnames = FALSE) <- log1p(as(exp, "CsparseMatrix"))
  
  # Filter features and samples
  if ("selected" %in% colnames(rowData(spe))) {
    spe <- spe[as.logical(rowData(spe)$selected), ]
  }
  if ("selected" %in% colnames(colData(spe))) {
    spe <- spe[, as.logical(colData(spe)$selected)]
  }
  
  return(spe)
}

# Seed
set.seed(seed)

spe <- get_SpatialExperiment(
  feature_file = feature_file,
  observation_file = observation_file,
  coord_file = coord_file,
  matrix_file = matrix_file
)
print('get_SpatialExpr done')

## Configuration
method <- "HMRF"
k <- config$k
dims_used <- config$dims_used

## Giotto instructions
python_path <- Sys.which(c("python"))
instrs <- createGiottoInstructions(save_dir = out_dir,
                                   save_plot = FALSE,
                                   show_plot = FALSE,
                                   python_path = python_path)

## raw expression counts expected
createGiotto_fn = function(spe, annotation = FALSE, selected_clustering = NULL, instructions = NULL){
  raw_expr <- SummarizedExperiment::assay(spe, "counts")
  #colnames(raw_expr) <- colData(sce)[,"Barcode"]
  norm_expression <- SummarizedExperiment::assay(spe, "logcounts")
  
  cell_metadata <- SingleCellExperiment::colData(spe)
  cell_metadata$cell_ID <- rownames(SingleCellExperiment::colData(spe))
  print(cell_metadata)

  colnames(cell_metadata)[c(1,2)] <- c("sdimx", "sdimy")
  cell_metadata <- as.data.frame(cell_metadata[,c(5,1,2,3, 4)])
  feat_metadata <- tibble::as_tibble(SingleCellExperiment::rowData(spe),rownames = "feat_ID")
  if (annotation) {
    rownames(raw_expr) <- c(SingleCellExperiment::rowData(spe)[, "SYMBOL"])
    #rownames(norm_expression) <- c(SingleCellExperiment::rowData(sce)[,"SYMBOL"])
  }
  gobj = Giotto::createGiottoObject(
    expression = list("raw" = raw_expr#,
                      #"normalized" = norm_expression
    ),
    cell_metadata = cell_metadata,
    spatial_locs = as.data.frame(SpatialExperiment::spatialCoords(spe)),
    feat_metadata = feat_metadata,
    instructions = instrs#,
    # Add dimred (doesn't quite work)
    #dimension_reduction = GiottoClass::createDimObj(
    #    coordinates = SingleCellExperiment::reducedDim(spe, "reducedDim"),
    #    name = "PCA",
    #    method = "pca")
  )
  return(gobj)
}

# Convert to Giotto object
gobj <- createGiotto_fn(spe, instructions = instrs)
print('Converted to Giotto object')

# Normalize
gobj <- Giotto::normalizeGiotto(gobj)
# Alternatively, use the Giotto normalization
print('Normalized giotto')

# PCA
gobj <- runPCA(gobject = gobj, center = TRUE, scale_unit = TRUE, name = "PCA", feats_to_use = NULL)
print('PCA done')

message("Running ", method, " clustering")

gobj <- Giotto::createSpatialNetwork(
  gobject = gobj,
  minimum_k = 10
)
print('Spatial Network created')

# identify genes with a spatial coherent expression profile (Not necessary - default uses all 'selected' features)
# km_spatialgenes <- Giotto::binSpect(gobj, bin_method = 'rank')
# print('BinSpect created')
HMRF_init_obj <- initHMRF_V2(gobject = g, cl.method = "km")
install.packages('tidygraph')
print(HMRF_init_obj)
#my_spatial_genes <- km_spatialgenes[1:100]$feats
HMRF_spatial_genes <- Giotto::doHMRF_V2 (
  HMRF_init_obj = HMRF_init_obj,
  # spat_unit = "cell",
  # feat_type = "rna",
  betas = c(0, 2, config$beta),
  # spatial_dimensions = c("sdimx", "sdimy"),
  # expression_values = "scaled", # not used when dim_reduction_to_use is given
  # spatial_genes = km_spatialgenes,
  # dim_reduction_to_use = "pca",
  # dim_reduction_name = "PCA",
  # dimensions_to_use = 1,
  # k = n_clusters,
  # name = method,
  # seed = seed
)
print('HMRF ')


## Add HMRF
gobj <- addHMRF(
  gobject = gobj,
  HMRFoutput = HMRF_spatial_genes,
  k = k, betas_to_add = c(config$beta),
  hmrf_name = method
)

## Write output
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
label_df <- as.data.frame(Giotto::getCellMetadata(gobj, output = "data.table"))
label_df <- data.frame(label = label_df[length(colnames(label_df))], row.names = label_df[[1]])
colnames(label_df) <- "label"

#print(table(label_df$label))

write.table(label_df, file = label_file, sep = "\t", col.names = NA, quote = FALSE)

# Clean HMRF_output folder
system("rm -rf HMRF_output")
