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

example3 <- read.csv('./data/simulation_coordiantes/mosna/isolated_niches/net_dir/nodes_patient-patient.csv')
coordinates = example3[c('x','y')]
colnames(coordinates) <- c('row','col')

counts = example3[c('type1','type2')]
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
results <- spatialCluster(se, q = 2)

head(colData(results))


clusterPlot(results)

write.csv(results@colData, './data/simulation_coordiantes/BayesSpace/isolated_niches.csv')








example3 <- read.csv('./data/simulation_coordiantes/mosna/overlap/net_dir/nodes_patient-patient.csv')
coordinates = example3[c('x','y')]
colnames(coordinates) <- c('row','col')

counts = example3[c('type1','type2')]
N = ncol(counts)
counts['tmp1'] <- runif(N)
counts['tmp2'] <- runif(N)
counts['tmp3'] <- runif(N)
counts['tmp4'] <- runif(N)
counts['tmp5'] <- runif(N)
counts['tmp6'] <- runif(N)
counts['tmp7'] <- runif(N)
counts['tmp8'] <- runif(N)
counts['tmp9'] <- runif(N)
counts['tmp10'] <- runif(N)
counts = t(counts)
rownames(counts) <- paste0("Gene_", seq_len(nrow(counts)))

# Create a SpatialExperiment object
se <- SpatialExperiment(
  assays = list(counts = counts),
  colData = DataFrame(col = coordinates[,1], row = coordinates[,2])
)

rownames(se) <- rownames(counts)
rownames(se)

se <- spatialPreprocess(se, platform="Visium",
                        n.PCs = 2,
                        n.HVGs = 6,
                        skip.PCA = FALSE,
                        # log.normalize = FALSE,
)
# 'res' controls the granularity of the clustering
results <- spatialCluster(se, q = 2)

head(colData(results))


clusterPlot(results)

write.csv(results@colData, './data/simulation_coordiantes/BayesSpace/overlap')














example3 <- read.csv('./data/simulation_coordiantes/mosna/gradual/net_dir/nodes_patient-patient.csv')
coordinates = example3[c('x','y')]
colnames(coordinates) <- c('row','col')

counts = example3[c('type1','type2')]
N = ncol(counts)
counts['tmp1'] <- runif(N)
counts['tmp2'] <- runif(N)
counts['tmp3'] <- runif(N)
counts['tmp4'] <- runif(N)
counts['tmp5'] <- runif(N)
counts['tmp6'] <- runif(N)
counts['tmp7'] <- runif(N)
counts['tmp8'] <- runif(N)
counts['tmp9'] <- runif(N)
counts['tmp10'] <- runif(N)
counts = t(counts)
rownames(counts) <- paste0("Gene_", seq_len(nrow(counts)))

# Create a SpatialExperiment object
se <- SpatialExperiment(
  assays = list(counts = counts),
  colData = DataFrame(col = coordinates[,1], row = coordinates[,2])
)

rownames(se) <- rownames(counts)
rownames(se)

se <- spatialPreprocess(se, platform="Visium",
                        n.PCs = 2,
                        n.HVGs = 6,
                        skip.PCA = FALSE,
                        # log.normalize = FALSE,
)
# 'res' controls the granularity of the clustering
results <- spatialCluster(se, q = 2)

head(colData(results))

clusterPlot(results)

write.csv(results@colData, './data/simulation_coordiantes/BayesSpace/gradual.csv')










example3 <- read.csv('./data/simulation_coordiantes/mosna/doublets/net_dir/nodes_patient-patient.csv')
coordinates = example3[c('x','y')]
colnames(coordinates) <- c('row','col')

counts = example3[c('A','B','C','D')]
N = ncol(counts)
counts['tmp1'] <- runif(N)
counts['tmp2'] <- runif(N)
counts['tmp3'] <- runif(N)
counts['tmp4'] <- runif(N)
counts['tmp5'] <- runif(N)
counts['tmp6'] <- runif(N)
counts['tmp7'] <- runif(N)
counts['tmp8'] <- runif(N)
counts['tmp9'] <- runif(N)
counts['tmp10'] <- runif(N)
counts = t(counts)
rownames(counts) <- paste0("Gene_", seq_len(nrow(counts)))

# Create a SpatialExperiment object
se <- SpatialExperiment(
  assays = list(counts = counts),
  colData = DataFrame(col = coordinates[,1], row = coordinates[,2])
)

rownames(se) <- rownames(counts)
rownames(se)

se <- spatialPreprocess(se, platform="Visium",
                        n.PCs = 2,
                        n.HVGs = 6,
                        skip.PCA = FALSE,
                        # log.normalize = FALSE,
)
# 'res' controls the granularity of the clustering
results <- spatialCluster(se, q = 2)

head(colData(results))

clusterPlot(results)

write.csv(results@colData, './data/simulation_coordiantes/BayesSpace/doublets.csv')


X <- results@colData





















example3 <- read.csv('./data/simulation_coordiantes/mosna/boundries/net_dir/nodes_patient-patient.csv')
coordinates = example3[c('x','y')]
colnames(coordinates) <- c('row','col')

counts = example3[c('A','B','C','D')]
N = nrow(counts)
counts['tmp1'] <- runif(N)
counts['tmp2'] <- runif(N)
counts['tmp3'] <- runif(N)
counts['tmp4'] <- runif(N)
counts = t(counts)
rownames(counts) <- paste0("Gene_", seq_len(nrow(counts)))


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
results <- spatialCluster(se, q = 2)

head(colData(results))


clusterPlot(results)

write.csv(results@colData, './data/simulation_coordiantes/BayesSpace/boundries.csv')








example3 <- read.csv('./data/simulation_coordiantes/mosna/rare/net_dir/nodes_patient-patient.csv')
coordinates = example3[c('x','y')]
colnames(coordinates) <- c('row','col')

counts = example3[c('A','B','C')]
N = nrow(counts)
counts['tmp1'] <- runif(N)
counts['tmp2'] <- runif(N)
counts['tmp3'] <- runif(N)
counts['tmp4'] <- runif(N)
counts = t(counts)
rownames(counts) <- paste0("Gene_", seq_len(nrow(counts)))


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
results <- spatialCluster(se, q = 3)

head(colData(results))


clusterPlot(results)

write.csv(results@colData, './data/simulation_coordiantes/BayesSpace/rare.csv')


