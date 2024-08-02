if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# BiocManager::install("imcRtools")
# BiocManager::install("spatstat")
# BiocManager::install("imcRtools")
remotes::install_github("BodenmillerGroup/imcRtools")

example3 <- read.csv('./data/simulation_coordiantes/mosna/isolated_niches/net_dir/nodes_patient-patient.csv')

# Load libraries
library(imcRtools)
library(dplyr)
library(ggplot2)
library(viridis)
library(spatstat)

df <- data.frame(
  x = example3$x,
  y = example3$y,
  CellType = example3$phenotype,
  ImageNb = rep('1',nrow(example3)),
  marker1 = example3$type1,
  marker2 = example3$type2
)

# Convert marker data to a matrix
marker_data <- as.matrix(df[, grep("marker", names(df))])
marker_data <- t(marker_data) # Transposing so columns are cells

# Create SCE object
sce <- SingleCellExperiment(assays = list(counts = marker_data))

# Add metadata
colData(sce) <- DataFrame(x = df$x, y = df$y, coords = df[c('x','y')] , CellType = df$CellType, ImageNb = df$ImageNb)
rowData(sce) <- DataFrame(markerNames = rownames(marker_data))
colData(sce)
# Build spatial graphs
sce <- buildSpatialGraph(sce, coords = c('x','y'), img_id = "ImageNb", type = "expansion", threshold = 20)
sce <- buildSpatialGraph(sce, coords = c('x','y'), img_id = "ImageNb", type = "knn", k = 5)
sce <- buildSpatialGraph(sce, coords = c('x','y'), img_id = "ImageNb", type = "delaunay")

# Aggregate neighbors
sce <- aggregateNeighbors(sce, colPairName = "knn_interaction_graph", aggregate_by = "metadata", count_by = "CellType")

# View results
head(sce$aggregatedNeighbors)

# Clustering
set.seed(22)
cur_cluster <- kmeans(sce$aggregatedNeighbors, centers = 2)
sce$clustered_neighbors <- factor(cur_cluster$cluster)

# Visualization
library(ggraph)
plotSpatial(sce,
            img_id = "ImageNb",
            node_color_by = "CellType",
            node_size_fix = 4,
            coords = c("x", "y"),
            
            edge_width_fix = 2,
            edge_color_by = "clustered_neighbors",
            draw_edges = TRUE,
            colPairName = "knn_interaction_graph",
            directed = FALSE,
            nodes_first = FALSE,
            scales = "free") +
  scale_color_brewer(palette = "Set2") +
  scale_edge_color_brewer(palette = "Set1")


write.csv(sce@colData , './data/simulation_coordiantes/imcRtools/isolated_niches.csv')

TMP <- sce@colData


####### PATCHES
sce2 <- patchDetection(sce, 
                       coords = c('x','y'),
                              patch_cells = sce$CellType %in% c("X"),
                              colPairName = "expansion_interaction_graph",
                              expand_by = .005, 
                              img_id = "ImageNb")

plotSpatial(sce2, 
            coords = c('x','y'),
            img_id = "ImageNb", 
            node_color_by = "patch_id",
            scales = "free")










example3 <- read.csv('./data/simulation_coordiantes/mosna/gradual/net_dir/nodes_patient-patient.csv')

df <- data.frame(
  x = example3$x,
  y = example3$y,
  CellType = example3$phenotype,
  ImageNb = rep('1',nrow(example3)),
  marker1 = example3$type1,
  marker2 = example3$type2
)

# Convert marker data to a matrix
marker_data <- as.matrix(df[, grep("marker", names(df))])
marker_data <- t(marker_data) # Transposing so columns are cells

# Create SCE object
sce <- SingleCellExperiment(assays = list(counts = marker_data))

# Add metadata
colData(sce) <- DataFrame(x = df$x, y = df$y, coords = df[c('x','y')] , CellType = df$CellType, ImageNb = df$ImageNb)
rowData(sce) <- DataFrame(markerNames = rownames(marker_data))
colData(sce)
# Build spatial graphs
sce <- buildSpatialGraph(sce, coords = c('x','y'), img_id = "ImageNb", type = "expansion", threshold = 20)
sce <- buildSpatialGraph(sce, coords = c('x','y'), img_id = "ImageNb", type = "knn", k = 5)
sce <- buildSpatialGraph(sce, coords = c('x','y'), img_id = "ImageNb", type = "delaunay")

# Aggregate neighbors
sce <- aggregateNeighbors(sce, colPairName = "knn_interaction_graph", aggregate_by = "metadata", count_by = "CellType")

# View results
head(sce$aggregatedNeighbors)

# Clustering
set.seed(22)
cur_cluster <- kmeans(sce$aggregatedNeighbors, centers = 5)
sce$clustered_neighbors <- factor(cur_cluster$cluster)

# Visualization
library(ggraph)
plotSpatial(sce,
            img_id = "ImageNb",
            node_color_by = "CellType",
            node_size_fix = 4,
            coords = c("x", "y"),
            
            edge_width_fix = 2,
            edge_color_by = "clustered_neighbors",
            draw_edges = TRUE,
            colPairName = "knn_interaction_graph",
            directed = FALSE,
            nodes_first = FALSE,
            scales = "free") +
  scale_color_brewer(palette = "Set2") +
  scale_edge_color_brewer(palette = "Set1")


write.csv(sce@colData , './data/simulation_coordiantes/imcRtools/gradual.csv')













example3 <- read.csv('./data/simulation_coordiantes/mosna/doublets/net_dir/nodes_patient-patient.csv')

df <- data.frame(
  x = example3$x,
  y = example3$y,
  CellType = example3$phenotype,
  ImageNb = rep('1',nrow(example3)),
  marker1 = example3$A,
  marker3 = example3$C,
  marker4 = example3$D,
  marker2 = example3$B
)

# Convert marker data to a matrix
marker_data <- as.matrix(df[, grep("marker", names(df))])
marker_data <- t(marker_data) # Transposing so columns are cells

# Create SCE object
sce <- SingleCellExperiment(assays = list(counts = marker_data))

# Add metadata
colData(sce) <- DataFrame(x = df$x, y = df$y, coords = df[c('x','y')] , CellType = df$CellType, ImageNb = df$ImageNb)
rowData(sce) <- DataFrame(markerNames = rownames(marker_data))
colData(sce)
# Build spatial graphs
sce <- buildSpatialGraph(sce, coords = c('x','y'), img_id = "ImageNb", type = "expansion", threshold = 20)
sce <- buildSpatialGraph(sce, coords = c('x','y'), img_id = "ImageNb", type = "knn", k = 5)
sce <- buildSpatialGraph(sce, coords = c('x','y'), img_id = "ImageNb", type = "delaunay")

# Aggregate neighbors
sce <- aggregateNeighbors(sce, colPairName = "knn_interaction_graph", aggregate_by = "metadata", count_by = "CellType")

# View results
head(sce$aggregatedNeighbors)

# Clustering
set.seed(22)
cur_cluster <- kmeans(sce$aggregatedNeighbors, centers = 2)
sce$clustered_neighbors <- factor(cur_cluster$cluster)

# Visualization
library(ggraph)
plotSpatial(sce,
            img_id = "ImageNb",
            node_color_by = "CellType",
            node_size_fix = 4,
            coords = c("x", "y"),
            
            edge_width_fix = 2,
            edge_color_by = "clustered_neighbors",
            draw_edges = TRUE,
            colPairName = "knn_interaction_graph",
            directed = FALSE,
            nodes_first = FALSE,
            scales = "free") +
  scale_color_brewer(palette = "Set2") +
  scale_edge_color_brewer(palette = "Set1")


write.csv(sce@colData , './data/simulation_coordiantes/imcRtools/doublets.csv')

















example3 <- read.csv('./data/simulation_coordiantes/mosna/boundries/net_dir/nodes_patient-patient.csv')

df <- data.frame(
  x = example3$x,
  y = example3$y,
  CellType = example3$phenotype,
  ImageNb = rep('1',nrow(example3)),
  marker1 = example3$A,
  marker3 = example3$C,
  marker4 = example3$D,
  marker2 = example3$B
)

# Convert marker data to a matrix
marker_data <- as.matrix(df[, grep("marker", names(df))])
marker_data <- t(marker_data) # Transposing so columns are cells

# Create SCE object
sce <- SingleCellExperiment(assays = list(counts = marker_data))

# Add metadata
colData(sce) <- DataFrame(x = df$x, y = df$y, coords = df[c('x','y')] , CellType = df$CellType, ImageNb = df$ImageNb)
rowData(sce) <- DataFrame(markerNames = rownames(marker_data))
colData(sce)
# Build spatial graphs
sce <- buildSpatialGraph(sce, coords = c('x','y'), img_id = "ImageNb", type = "expansion", threshold = 20)
sce <- buildSpatialGraph(sce, coords = c('x','y'), img_id = "ImageNb", type = "knn", k = 5)
sce <- buildSpatialGraph(sce, coords = c('x','y'), img_id = "ImageNb", type = "delaunay")

# Aggregate neighbors
sce <- aggregateNeighbors(sce, colPairName = "knn_interaction_graph", aggregate_by = "metadata", count_by = "CellType")

# View results
head(sce$aggregatedNeighbors)

# Clustering
set.seed(22)
cur_cluster <- kmeans(sce$aggregatedNeighbors, centers = 2)
sce$clustered_neighbors <- factor(cur_cluster$cluster)

# Visualization
library(ggraph)
plotSpatial(sce,
            img_id = "ImageNb",
            node_color_by = "CellType",
            node_size_fix = 4,
            coords = c("x", "y"),
            
            edge_width_fix = 2,
            edge_color_by = "clustered_neighbors",
            draw_edges = TRUE,
            colPairName = "knn_interaction_graph",
            directed = FALSE,
            nodes_first = FALSE,
            scales = "free") +
  scale_color_brewer(palette = "Set2") +
  scale_edge_color_brewer(palette = "Set1")


write.csv(sce@colData , './data/simulation_coordiantes/imcRtools/boundries.csv')





















example3 <- read.csv('./data/simulation_coordiantes/mosna/rare/net_dir/nodes_patient-patient.csv')

df <- data.frame(
  x = example3$x,
  y = example3$y,
  CellType = example3$phenotype,
  ImageNb = rep('1',nrow(example3)),
  marker1 = example3$A,
  marker3 = example3$C,
  # marker4 = example3$D,
  marker2 = example3$B
)

# Convert marker data to a matrix
marker_data <- as.matrix(df[, grep("marker", names(df))])
marker_data <- t(marker_data) # Transposing so columns are cells

# Create SCE object
sce <- SingleCellExperiment(assays = list(counts = marker_data))

# Add metadata
colData(sce) <- DataFrame(x = df$x, y = df$y, coords = df[c('x','y')] , CellType = df$CellType, ImageNb = df$ImageNb)
rowData(sce) <- DataFrame(markerNames = rownames(marker_data))
colData(sce)
# Build spatial graphs
sce <- buildSpatialGraph(sce, coords = c('x','y'), img_id = "ImageNb", type = "expansion", threshold = 20)
sce <- buildSpatialGraph(sce, coords = c('x','y'), img_id = "ImageNb", type = "knn", k = 5)
sce <- buildSpatialGraph(sce, coords = c('x','y'), img_id = "ImageNb", type = "delaunay")

# Aggregate neighbors
sce <- aggregateNeighbors(sce, colPairName = "knn_interaction_graph", aggregate_by = "metadata", count_by = "CellType")

# View results
head(sce$aggregatedNeighbors)

# Clustering
set.seed(22)
cur_cluster <- kmeans(sce$aggregatedNeighbors, centers = 2)
sce$clustered_neighbors <- factor(cur_cluster$cluster)

# Visualization
library(ggraph)
plotSpatial(sce,
            img_id = "ImageNb",
            node_color_by = "CellType",
            node_size_fix = 4,
            coords = c("x", "y"),
            
            edge_width_fix = 2,
            edge_color_by = "clustered_neighbors",
            draw_edges = TRUE,
            colPairName = "knn_interaction_graph",
            directed = FALSE,
            nodes_first = FALSE,
            scales = "free") +
  scale_color_brewer(palette = "Set2") +
  scale_edge_color_brewer(palette = "Set1")


write.csv(sce@colData , './data/simulation_coordiantes/imcRtools/rare.csv')













###### SPAT

example3 <- read.csv('./data/simulation_coordiantes/mosna/example_3_niches/net_dir/nodes_patient-patient.csv', row.names = 1)

df <- data.frame(
  x = example3$x,
  y = example3$y,
  CellType = example3$cell_type,
  ImageNb = rep('1',nrow(example3)),
  marker1 = example3$A,
  marker2 = example3$B,
  marker3 = example3$C
)

# Convert marker data to a matrix
marker_data <- as.matrix(df[, grep("marker", names(df))])
marker_data <- t(marker_data) # Transposing so columns are cells

# Create SCE object
sce <- SingleCellExperiment(assays = list(counts = marker_data))

# Add metadata
colData(sce) <- DataFrame(x = df$x, y = df$y, coords = df[c('x','y')] , CellType = df$CellType, ImageNb = df$ImageNb)
rowData(sce) <- DataFrame(markerNames = rownames(marker_data))
colData(sce)
# Build spatial graphs
sce <- buildSpatialGraph(sce, coords = c('x','y'), img_id = "ImageNb", type = "expansion", threshold = 20)
sce <- buildSpatialGraph(sce, coords = c('x','y'), img_id = "ImageNb", type = "knn", k = 5)
sce <- buildSpatialGraph(sce, coords = c('x','y'), img_id = "ImageNb", type = "delaunay")

# Aggregate neighbors
sce <- aggregateNeighbors(sce, colPairName = "knn_interaction_graph", aggregate_by = "metadata", count_by = "CellType")

# View results
head(sce$aggregatedNeighbors)

# Clustering
set.seed(22)
cur_cluster <- kmeans(sce$aggregatedNeighbors, centers = 5)
sce$clustered_neighbors <- factor(cur_cluster$cluster)

# Visualization
library(ggraph)
plotSpatial(sce,
            img_id = "ImageNb",
            node_color_by = "CellType",
            node_size_fix = 4,
            coords = c("x", "y"),
            
            edge_width_fix = 2,
            edge_color_by = "clustered_neighbors",
            draw_edges = TRUE,
            colPairName = "knn_interaction_graph",
            directed = FALSE,
            nodes_first = FALSE,
            scales = "free") +
  scale_color_brewer(palette = "Set2") +
  scale_edge_color_brewer(palette = "Set1")


write.csv(sce@colData , './data/simulation_coordiantes/imcRtools/aggregateNeighbors-Dep.csv')








