# remotes::install_github('JEFworks-Lab/MERINGUE')
suppressMessages(library(MERINGUE))
data(mOB)

example3 <- read.csv('./data/simulation_coordiantes/mosna/example_3_niches/net_dir/nodes_patient-patient.csv', row.names = 1)
counts <- example3[c('X','Y','Z')]

# CPM normalize
mat <- normalizeCounts(counts = t(counts), 
                       log=FALSE,
                       verbose=TRUE)

# Dimensionality reduction by PCA on log10 CPM expression values
pcs.info <- prcomp(t(log10(as.matrix(mat)+1)), center=TRUE)
nPcs <- 3
pcs <- pcs.info$x[,1:nPcs]

# 2D embedding by tSNE
emb <- Rtsne::Rtsne(pcs,
                    is_distance=FALSE,
                    check_duplicates = FALSE,
                    perplexity=3,
                    num_threads=1,
                    verbose=FALSE)$Y

rownames(emb) <- rownames(pcs)

# Graph-based cluster detection
k <- 30
com <- getClusters(pcs, k, weight=TRUE)

# Manually annotate identified clusters with cell-types
annot <- as.character(com); names(annot) <- names(com)


# Plot
par(mfrow=c(1,2), mar=rep(1,4))
plotEmbedding(emb, groups=annot, 
              show.legend=TRUE, xlab=NA, ylab=NA,
              verbose=FALSE)
plotEmbedding(pos, groups=annot, 
              cex=1, xlab=NA, ylab=NA,
              verbose=FALSE)



# Sample 2000 genes for demo purposes only to minimize runtime for demo only
set.seed(0)
test <- mat

# Identify significantly differentially upregulated genes
# in each identified cluster by Wilcox test
dg <- getDifferentialGenes(as.matrix(mat[test,]), annot)

dg.sig <- lapply(dg, function(x) {
  x <- x[x$p.adj < 0.05,]
  x <- na.omit(x)
  x <- x[x$highest,]
  rownames(x)
})
print(lapply(dg.sig, length))



dg.genes <- unlist(dg.sig)
ggroup <- unlist(lapply(1:length(dg.sig), function(i) { 
  rep(names(dg.sig)[i], length(dg.sig[[i]]))
}))
names(ggroup) <- dg.genes
ggroup <- factor(ggroup)

# Plot
ccol <- rainbow(length(levels(annot)))[annot]
names(ccol) <- names(annot) # column colors

lgcol <- rainbow(length(levels(ggroup)), v=0.5)[ggroup]
names(gcol) <- names(ggroup) # row colors

m <- as.matrix(mat[dg.genes, names(sort(annot))])
m <- winsorize(t(scale(t(m))))
heatmap(m, scale="none", 
        Colv=NA, Rowv=NA, labRow=NA, labCol=NA,
        ColSideColors=ccol[colnames(m)],
        RowSideColors=gcol[rownames(m)],
        col=colorRampPalette(c('blue', 'white', 'red'))(100)
)



# Get neighbor-relationships
w <- getSpatialNeighbors(pos, filterDist = 2.5)
plotNetwork(pos, w)





# Identify sigificantly spatially auto-correlated genes
start_time <- Sys.time()
I <- getSpatialPatterns(mat[test,], w)
end_time <- Sys.time()
print(end_time - start_time)

## Time difference of 1.59313 secs



results.filter <- filterSpatialPatterns(mat = mat[test,],
                                        I = I,
                                        w = w,
                                        adjustPv = TRUE,
                                        alpha = 0.05,
                                        minPercentCells = 0.05,
                                        verbose = TRUE)

## Number of significantly autocorrelated genes: 256

## ...driven by > 13 cells: 230




# Compute spatial cross correlation matrix
scc <- spatialCrossCorMatrix(mat = as.matrix(mat[results.filter,]), 
                             weight = w)
# Identify primary patterns
par(mfrow=c(2,2), mar=rep(2,4))
ggroup <- groupSigSpatialPatterns(pos = pos, 
                                  mat = as.matrix(mat[results.filter,]), 
                                  scc = scc, 
                                  power = 1, 
                                  hclustMethod = 'ward.D', 
                                  deepSplit = 2,
                                  zlim=c(-1.5,1.5))
