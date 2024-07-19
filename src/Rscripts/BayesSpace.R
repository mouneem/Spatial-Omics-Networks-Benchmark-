library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
# BiocManager::install('BayesSpace', force)

library(SpatialExperiment)

# Example data loading
# Assume 'data' is a SpatialExperiment object, you would normally load your data like this:
# data <- readRDS("path_to_your_data.rds")

# For demonstration, let's create a mock SpatialExperiment with random data
set.seed(123)

example3 <- read.csv('./data/simulation_coordiantes/mosna/example_3_niches/net_dir/nodes_patient-patient.csv', row.names = 1)
coordinates = example3[c('x','y')]
colnames(coordinates) <- c('row','col')

counts = example3[c('X','Y','Z')]
N = ncol(counts)
counts['tmp1'] <- runif(N)
counts['tmp2'] <- runif(N)
counts['tmp3'] <- runif(N)
counts['tmp4'] <- runif(N)
rownames(counts) <- paste0("Gene_", seq_len(nrow(counts)))
counts = t(counts)


# Create a SpatialExperiment object
se <- SpatialExperiment(
  assays = list(counts = counts),
  colData = DataFrame(col = coordinates[,1], row = coordinates[,2])
)

rownames(se) <- rownames(counts)
rownames(se)

se <- spatialPreprocess(se, platform="Visium",
                        n.PCs = 2,
                        n.HVGs = 5,
                        skip.PCA = FALSE,
                        # log.normalize = FALSE,
                        )
# 'res' controls the granularity of the clustering
results <- spatialCluster(se, q = 4)

head(colData(results))


clusterPlot(results)

write.csv(results@colData, './data/simulation_coordiantes/BayesSpace/example_3.csv')








example3 <- read.csv('./data/simulation_coordiantes/mosna/example_3_niches/net_dir/nodes_patient-patient.csv', row.names = 1)
coordinates = example3[c('x','y')]
colnames(coordinates) <- c('row','col')

counts = example3[c('A','B','C')]
N = ncol(counts)
counts['tmp1'] <- runif(N)
counts['tmp2'] <- runif(N)
counts['tmp3'] <- runif(N)
counts['tmp4'] <- runif(N)
rownames(counts) <- paste0("Gene_", seq_len(nrow(counts)))
counts = t(counts)


# Create a SpatialExperiment object
se <- SpatialExperiment(
  assays = list(counts = counts),
  colData = DataFrame(col = coordinates[,1], row = coordinates[,2])
)

rownames(se) <- rownames(counts)
rownames(se)

se <- spatialPreprocess(se, platform="Visium",
                        n.PCs = 2,
                        n.HVGs = 5,
                        skip.PCA = FALSE,
                        # log.normalize = FALSE,
)
# 'res' controls the granularity of the clustering
results <- spatialCluster(se, q = 6)

head(colData(results))


clusterPlot(results)

write.csv(results@colData, './data/simulation_coordiantes/BayesSpace/example_3_dep.csv')


