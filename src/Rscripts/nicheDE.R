devtools::install_github("kaishumason/NicheDE") # install

example3 <- read.csv('./data/simulation_coordiantes/mosna/example_3_niches/net_dir/nodes_patient-patient.csv', row.names = 1)


library(nicheDE)
# Assuming df is your initial dataframe with coordinates (x, y) and markers (gene1, gene2, ...)
coordinates <- example3[, c("x", "y")]
markers <- example3[, c("X", "Y", "Z")]  # adjust column names as per your dataframe

library(SingleCellExperiment)

# Example dataframe structure
df <- data.frame(
  x = c(1, 2, 3, 4, 5),
  y = c(1, 2, 3, 4, 5),
  marker1 = rnorm(5),
  marker2 = rnorm(5),
  marker3 = rnorm(5)
)

# Create the SCE object
sce <- SingleCellExperiment(assays = list(counts = as.matrix(df[, -c(1, 2)])))

# Add spatial coordinates
spatialCoords(sce) <- df[, 1:2, drop = FALSE]
