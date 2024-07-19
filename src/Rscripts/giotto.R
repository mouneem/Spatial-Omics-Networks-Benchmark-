print('running giotto from ./Rscripts/')

library(Giotto)
# Ensure the Python environment for Giotto has been installed.
genv_exists = checkGiottoEnvironment()
if(!genv_exists){
  # The following command need only be run once to install the Giotto environment.
  installGiottoEnvironment()
}

library(GiottoData)

# Specify path from which data may be retrieved/stored
data_directory = paste0(getwd(),'/gobject_visual_data/')
# alternatively, "/path/to/where/the/data/lives/"

# Specify path to which results may be saved
results_directory = paste0(getwd(),'/gobject_visual_results/')
# alternatively, "/path/to/store/the/results/"

# Optional: Specify a path to a Python executable within a conda or miniconda
# environment. If set to NULL (default), the Python executable within the previously
# installed Giotto environment will be used.
my_python_path = NULL # alternatively, "/local/python/path/python" if desired.

# Get the dataset
getSpatialDataset(dataset = 'merfish_preoptic', directory = data_directory, method = 'curl')

### Giotto instructions and data preparation
# Optional: Set Giotto instructions
instrs = createGiottoInstructions(save_plot = TRUE,
                                  show_plot = TRUE,
                                  save_dir = results_directory,
                                  python_path = my_python_path)

# Create file paths to feed data into Giotto object
expr_path = paste0(data_directory, "merFISH_3D_data_expression.txt.gz")
loc_path = paste0(data_directory, "merFISH_3D_data_cell_locations.txt")
meta_path = paste0(data_directory, "merFISH_3D_metadata.txt")

### Create Giotto object
testobj <- createGiottoObject(expression = expr_path,
                              spatial_locs = loc_path,
                              instructions = instrs)


# Add additional metadata
metadata = data.table::fread(meta_path)

testobj = addCellMetadata(testobj,
                          new_metadata = metadata$layer_ID,
                          vector_name = 'layer_ID')

testobj = addCellMetadata(testobj,
                          new_metadata = metadata$orig_cell_types,
                          vector_name = 'orig_cell_types')

### Process the Giotto Object
# Note that for the purposes of this tutorial, the entire dataset will be visualized.
# Thus, filter parameters are set to 0, so as to not remove any cells.
# Note that since adjustment is not required, adjust_params is set to NULL.

testobj <- processGiotto(testobj,
                         filter_params = list(expression_threshold = 0,
                                              feat_det_in_min_cells = 0,
                                              min_det_feats_per_cell = 0),
                         norm_params = list(norm_methods = 'standard',
                                            scale_feats = TRUE,
                                            scalefactor = 1000),
                         stat_params = list(expression_values = 'normalized'),
                         adjust_params = NULL)




# Create a spatial network (optional, depending on your analysis)
giotto_object <- testobj

# Extract statistics
# Example: Calculate cell type proportions
cell_type_proportions <- table(giotto_object@cell_metadata$cell) / nrow(giotto_object@cell_metadata)

giotto_object <- createSpatialNetwork(gobject = giotto_object, minimum_k = 2)
# Example: Calculate spatial enrichment for cell types
spatial_enrichment <- binSpect(giotto_object)

# Example: Calculate nearest neighbor network statistics
spatial_network_stats <- cellProximityEnrichment(gobject = giotto_object, spatial_network_name = 'Delaunay_network', cluster_column = 'cell_type')

