import pandas as pd
import seaborn as sns
import os
import random
import math
import matplotlib.pyplot as plt

print('Script loaded!')

def distance(p1, p2):
    return math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)




def select_annotation(annotation):
    total_ratio = sum(annotation.values())
    rand_num = random.uniform(0, total_ratio)
    cumulative_ratio = 0
    for label, ratio in annotation.items():
        cumulative_ratio += ratio
        if rand_num < cumulative_ratio:
            return label

    # In case of rounding issues, return the last label
    return label


def generate_points(max_x, max_y, num_points, min_x=0, min_y=0, min_distance=None, equidistant=False, max_attempts_per_point=100, annotation=None):
    points = []
    annotation_list = []

    if equidistant:
        points = []
        grid_size = int(math.ceil(math.sqrt(num_points)))
        if min_distance is None:
            raise ValueError("Minimum distance must be specified for equidistant grid generation.")
        
        step_x = (max_x - min_x) / grid_size
        step_y = (max_y - min_y) / grid_size
        
        for i in range(grid_size):
            for j in range(grid_size):
                if len(points) < num_points:
                    points.append((min_x + i * step_x, min_y + j * step_y))
                    if annotation:
                        annotation_list.append(select_annotation(annotation))
        return points, annotation_list

    attempts = 0
    max_attempts = num_points * max_attempts_per_point

    while len(points) < num_points and attempts < max_attempts:
        new_point = (random.uniform(min_x, max_x), random.uniform(min_y, max_y))
        if min_distance is not None:
            if all(distance(new_point, point) >= min_distance for point in points):
                points.append(new_point)
                if annotation:
                    annotation_list.append(select_annotation(annotation))
        else:
            points.append(new_point)
            if annotation:
                annotation_list.append(select_annotation(annotation))
        
        attempts += 1

    if attempts >= max_attempts:
        raise ValueError(f"Could not place {num_points} points with the given constraints.")

    return points, annotation_list

def add_neighboor(point_coordinates, distance_range = [1,10], annotation=['C'], number_neighbors = [1,5]):
    x, y = point_coordinates
    if type(distance_range) is int :
        distance_range = [distance_range, distance_range]
    min_distance, max_distance = distance_range
    
    if type(number_neighbors) is list or  type(number_neighbors) is tuple  :
        number_neighbors = random.randint(number_neighbors[0], number_neighbors[0])  # Choose between 1 to 4 neighbors randomly

    neighbors = []
    types = []
    for _ in range(number_neighbors):
        distance = random.uniform(min_distance, max_distance)
        angle = random.uniform(0, 2 * math.pi)
        neighbor_x = x + distance * math.cos(angle)
        neighbor_y = y + distance * math.sin(angle)
        neighbors.append((neighbor_x, neighbor_y))
        types.append(random.choice(annotation))

    return neighbors, types

def neighboring_points(coordiantes, types_list, target_type, new_type, dist_range, ratio):
    # coordiantes: x and y coordiantes of input cells
    # types_list: all type of input cells
    # target type: which type of cells (from annotation) to target
    # ratio: ratio of neighboors to add for each target cell
    coordiantes = pd.DataFrame(coordiantes)
    # select type
    selected_coords = pd.Series(types_list).isin(target_type)
    coordiantes_select = coordiantes.loc[selected_coords,:]

    coords = []
    types = []

    for key, row in coordiantes_select.iterrows():
        x = row[0]
        y = row[1]
        nei, nei_types = add_neighboor([x,y], distance_range = dist_range, annotation=new_type, number_neighbors = ratio)
        coords =  coords + list(nei)
        types = types + list(nei_types)
    return coords, types



def export_df_to_csv(df, file_path):
    # Create the directory path if it doesn't exist
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    
    # Export the DataFrame to a CSV file
    df.to_csv(file_path, index=False)



def plot_network(points, edges, types = None, linewidth = 1):

    if type(points) is pd.DataFrame:
        # Plot scatter points
        sns.scatterplot(x=points.iloc[:, 0], y=points.iloc[:, 1], hue=types)
    else:
        sns.scatterplot(x=points[:, 0], y=points[:, 1], hue=types)

    # Plot the edges
    for edge in edges:
        try:
            point1, point2 = edge
            x_values = [points.iloc[point1,0], points.iloc[point2,0]]
            y_values = [points.iloc[point1,1], points.iloc[point2,1]]
            plt.plot(x_values, y_values, color = 'black', linewidth=linewidth)
        except:
            print(f"The edge: {edge} does not correspond to any points")
            print(points.iloc[point1,0], points.iloc[point1,1])
    # Remove axis and grid
    plt.axis('off')
    # Set labels and title
    plt.show()



import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay

import numpy as np
from scipy.spatial import Delaunay

def delaunay_edges(points):
    """
    Perform Delaunay triangulation on a set of points and return the edges.

    Parameters:
    - points (DataFrame or ndarray): Input points in 2D space. If DataFrame, 
      should have two columns representing x and y coordinates. If ndarray, 
      should have shape (n, 2) where n is the number of points.

    Returns:
    - edges (list of tuples): List of edges (as tuples of point indices) 
      representing the Delaunay triangulation of the input points.

    """
    # Convert DataFrame to ndarray if points is a DataFrame
    if isinstance(points, pd.DataFrame):
        points = points.to_numpy()
    
    tri = Delaunay(points)
    
    # Collect unique edges from the Delaunay triangulation
    edges = {tuple(sorted([simplex[i], simplex[(i + 1) % 3]])) 
             for simplex in tri.simplices for i in range(3)}
    
    return list(edges)



from scipy.spatial import distance_matrix
from scipy.sparse.csgraph import minimum_spanning_tree

def mst_edges(points):
    """
    Compute Minimum Spanning Tree (MST) edges for a set of points.

    Parameters:
    - points (DataFrame or ndarray): Input points in 2D space. If DataFrame, 
      should have two columns representing x and y coordinates. If ndarray, 
      should have shape (n, 2) where n is the number of points.

    Returns:
    - mst_edges (list of tuples): List of edges (as tuples of point indices) 
      representing the edges in the Minimum Spanning Tree of the input points.

    Example:
    >>> import pandas as pd
    >>> points_df = pd.DataFrame([[0, 0], [1, 0], [0, 1], [1, 1]], columns=['x', 'y'])
    >>> edges = mst_edges(points_df)
    >>> print(edges)
    [(0, 1), (0, 2), (1, 3)]
    """
    # Convert DataFrame to ndarray if points is a DataFrame
    if isinstance(points, pd.DataFrame):
        points = points.to_numpy()
    
    # Compute distance matrix
    dist_matrix = distance_matrix(points, points)
    
    # Compute Minimum Spanning Tree
    mst = minimum_spanning_tree(dist_matrix).toarray().astype(bool)
    
    # Extract edges from the MST
    def extract_edges(matrix):
        edges = []
        num_vertices = matrix.shape[0]
        for i in range(num_vertices):
            for j in range(i + 1, num_vertices):
                if matrix[i, j]:
                    edges.append((i, j))
        return edges
    
    mst_edges = extract_edges(mst)
    
    return mst_edges



def gabriel_edges(points):
    """
    Calculate Gabriel edges for a set of points in 2D space.

    Parameters:
    - points (DataFrame or ndarray): Input points in 2D space. If DataFrame, 
      should have two columns representing x and y coordinates. If ndarray, 
      should have shape (n, 2) where n is the number of points.

    Returns:
    - edges_gabriel (list of tuples): List of edges (as tuples of point indices) 
      representing the Gabriel graph edges of the input points.

    Example:
    >>> import pandas as pd
    >>> points_df = pd.DataFrame([[0, 0], [1, 0], [0, 1], [1, 1]], columns=['x', 'y'])
    >>> gabriel_edges(points_df)
    [(0, 2), (1, 3)]
    """
    # Convert DataFrame to ndarray if points is a DataFrame
    if isinstance(points, pd.DataFrame):
        points = points.to_numpy()
    
    edges_gabriel = []
    dist_matrix = distance_matrix(points, points)
    
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            d = dist_matrix[i, j] / 2
            mid_point = (points[i] + points[j]) / 2
            if np.all(np.linalg.norm(points - mid_point, axis=1) >= d):
                edges_gabriel.append((i, j))
    
    return edges_gabriel


def rng_edges(points):
    """
    Compute the edges of the Relative Neighborhood Graph (RNG) for a set of points.

    Parameters:
    - points (DataFrame or ndarray): Input points in 2D space. If DataFrame, 
      should have two columns representing x and y coordinates. If ndarray, 
      should have shape (n, 2) where n is the number of points.

    Returns:
    - edges_rng (list of tuples): List of edges (as tuples of point indices) 
      representing the Relative Neighborhood Graph (RNG) of the input points.
    """
    if isinstance(points, pd.DataFrame):
        points = points.to_numpy()


    n = len(points)
    edges_rng = []

    # Compute distance matrix
    dist_matrix = np.linalg.norm(points[:, np.newaxis] - points, axis=2)

    # Compute RNG edges
    for i in range(n):
        for j in range(i + 1, n):
            is_rng_edge = True
            for k in range(n):
                if k != i and k != j:
                    if dist_matrix[i, k] < dist_matrix[i, j] and dist_matrix[j, k] < dist_matrix[i, j]:
                        is_rng_edge = False
                        break
            if is_rng_edge:
                edges_rng.append((i, j))

    return edges_rng


from sklearn.neighbors import kneighbors_graph
import numpy as np


def multi_sample_edges(points, sample_col = 'sample', coords_col = ['X','Y'], method = 'knn', param = False):
    all_edges = []
    last_edge = 0
    if method == 'knn':
        if not param:
            param = 1
        for sample in points[sample_col].unique():
            points_sample = points.loc[points[sample_col]==sample,:][coords_col]
            edges = knn_edges(points_sample, k=param)
            edges = pd.DataFrame(edges) + last_edge
            all_edges.append(edges)
            last_edge = edges.max()
        return pd.concat(all_edges)


def knn_edges(points, k=3):
    """
    Create a k-nearest neighbors graph from points and extract the edges.

    Parameters:
    - points (DataFrame or ndarray): Input points in 2D space. If DataFrame, 
      should have two columns representing x and y coordinates. If ndarray, 
      should have shape (n, 2) where n is the number of points.
    - k (int, optional): Number of nearest neighbors for the k-nearest neighbors graph. Default is 3.

    Returns:
    - edges_knn (list of tuples): List of edges (as tuples of point indices) 
      representing the edges of the k-nearest neighbors graph.

    Example:
    >>> import pandas as pd
    >>> points_df = pd.DataFrame([[0, 0], [1, 0], [0, 1], [1, 1]], columns=['x', 'y'])
    >>> edges = knn_edges(points_df, k=2)
    >>> print(edges)
    [(0, 1), (0, 2), (1, 0), (1, 3), (2, 0), (2, 3), (3, 1), (3, 2)]
    """
    # Convert DataFrame to ndarray if points is a DataFrame
    if isinstance(points, pd.DataFrame):
        points = points.to_numpy()
    
    # Create k-nearest neighbors graph
    knn_graph = kneighbors_graph(points, k, mode='connectivity', include_self=False)
    
    # Extract edges from the graph
    def extract_edges(adj_matrix):
        edges = np.column_stack(adj_matrix.nonzero())
        return [(edges[i, 0], edges[i, 1]) for i in range(edges.shape[0])]
    
    edges_knn = extract_edges(knn_graph)
    
    return edges_knn

from sklearn.metrics import silhouette_score

def find_optimal_k(points, max_k=10):
    """
    Find the optimal number of nearest neighbors (k) using silhouette score.

    Parameters:
    - points (ndarray): Input points in 2D space with shape (n, 2), where n is the number of points.
    - max_k (int): Maximum number of nearest neighbors to consider. Default is 10.

    Returns:
    - best_k (int): Optimal number of nearest neighbors (k) based on silhouette score.
    """
    silhouette_scores = []
    for k in range(2, max_k + 1):
        knn_graph = kneighbors_graph(points, k, mode='connectivity', include_self=False)
        adj_matrix = knn_graph.toarray()
        # Calculate silhouette score
        silhouette = silhouette_score(points, np.argmax(adj_matrix, axis=1))
        silhouette_scores.append(silhouette)
    
    # Find the index of the maximum silhouette score
    best_k = np.argmax(silhouette_scores) + 2  # add 2 because k starts from 2
    
    return best_k


def optimal_knn_edges(points, max_k=10):
    """
    Create a k-nearest neighbors graph from points and extract the edges.

    Parameters:
    - points (DataFrame or ndarray): Input points in 2D space. If DataFrame, 
      should have two columns representing x and y coordinates. If ndarray, 
      should have shape (n, 2) where n is the number of points.
    - max_k (int): Maximum number of nearest neighbors to consider. Default is 10.

    Returns:
    - edges_knn (list of tuples): List of edges (as tuples of point indices) 
      representing the edges of the k-nearest neighbors graph.
    """
    # Convert DataFrame to ndarray if points is a DataFrame
    if isinstance(points, pd.DataFrame):
        points = points.to_numpy()
    
    # Find optimal k using silhouette score
    best_k = find_optimal_k(points, max_k)
    
    # Create k-nearest neighbors graph with optimal k
    knn_graph = kneighbors_graph(points, best_k, mode='connectivity', include_self=False)
    
    # Extract edges from the graph
    def extract_edges(adj_matrix):
        edges = np.column_stack(adj_matrix.nonzero())
        return [(edges[i, 0], edges[i, 1]) for i in range(edges.shape[0])]
    
    edges_knn = extract_edges(knn_graph)
    
    return edges_knn


def radius_edges(points, radius):
    """
    Construct edges based on a radius threshold in a radius-based graph.

    Parameters:
    - points (DataFrame or ndarray): Input points in 2D space. If DataFrame, 
      should have two columns representing x and y coordinates. If ndarray, 
      should have shape (n, 2) where n is the number of points.
    - radius (float): Radius threshold for edge construction.

    Returns:
    - edges_radius (list of tuples): List of edges (as tuples of point indices) 
      where the distance between points is less than the specified radius.

    """
    # Convert DataFrame to ndarray if points is a DataFrame
    if isinstance(points, pd.DataFrame):
        points = points.to_numpy()

    # Calculate distance matrix
    dist_matrix = distance_matrix(points, points)

    # Extract edges where distance is less than radius
    edges_radius = np.column_stack(np.where(dist_matrix < radius)).tolist()

    # Ensure each edge is represented with sorted indices
    edges_radius = [tuple(sorted(edge)) for edge in edges_radius if edge[0] != edge[1]]
    
    # Remove duplicates by converting to set and back to list
    edges_radius = list(set(edges_radius))
    
    return edges_radius



import numpy as np
from sklearn.metrics.pairwise import pairwise_distances

def rips_complex_edges(points, threshold):
    """
    Construct edges based on the Vietoris-Rips complex method with a specified threshold.

    Parameters:
    - points (DataFrame or ndarray): Input points in 2D space. If DataFrame, 
      should have two columns representing x and y coordinates. If ndarray, 
      should have shape (n, 2) where n is the number of points.
    - threshold (float): Distance threshold for edge construction.

    Returns:
    - edges_rips (list of tuples): List of edges (as tuples of point indices) 
      where the pairwise distance between points is less than the specified threshold.

    """
    # Convert DataFrame to ndarray if points is a DataFrame
    if isinstance(points, pd.DataFrame):
        points = points.to_numpy()

    # Calculate pairwise distances
    dist_matrix = pairwise_distances(points)

    # Extract edges where distance is less than threshold
    edges_rips = np.column_stack(np.where(dist_matrix < threshold)).tolist()

    # Ensure each edge is represented with sorted indices
    edges_rips = [tuple(sorted(edge)) for edge in edges_rips if edge[0] != edge[1]]
    
    # Remove duplicates by converting to set and back to list
    edges_rips = list(set(edges_rips))
    
    return edges_rips



def plot_facetgrid_networks(edges_df, col_wrap = 4, height=4):
    """
    Plot scatter plots and connecting lines for edges data using Seaborn's FacetGrid.

    Parameters:
    - edges_df (DataFrame): DataFrame containing edge data with columns 'x1', 'y1', 'x2', 'y2', and 'method'.

    """
    # Initialize FacetGrid
    g = sns.FacetGrid(edges_df, col='method', height=height, col_wrap=col_wrap)

    # Define plotting function
    def plot_scatter_and_line(data, **kwargs):
        ax = plt.gca()
        for index, row in data.iterrows():
            ax.scatter(row['x1'], row['y1'])  # Scatter plot for (x1, y1)
            ax.scatter(row['x2'], row['y2'])  # Scatter plot for (x2, y2)
            ax.plot([row['x1'], row['x2']], [row['y1'], row['y2']], color='gray')  # Line connecting (x1, y1) to (x2, y2)

    # Apply plotting function to FacetGrid
    g.map_dataframe(plot_scatter_and_line)
    g.set_titles(col_template='{col_name}')

    # Adjust layout and display plot
    plt.tight_layout()
    plt.show()



import numpy as np
from scipy.spatial import cKDTree

def EpsNet_edges(points, epsilon):
    """
    Compute ε-Nets edges for a set of points.

    Parameters:
    - points (DataFrame or ndarray): Input points in 2D space. If DataFrame, 
      should have two columns representing x and y coordinates. If ndarray, 
      should have shape (n, 2) where n is the number of points.
    - epsilon (float): Maximum distance ε defining the ε-Nets.

    Returns:
    - eps_nets (list of tuples): List of edges (as tuples of point indices) 
      representing the ε-Nets edges.

    """
    # Convert DataFrame to ndarray if points is a DataFrame
    if isinstance(points, pd.DataFrame):
        points = points.to_numpy()
    
    # Build a KDTree for efficient nearest neighbor search
    tree = cKDTree(points)
    
    eps_nets = set()
    
    # Iterate through each point and find its neighbors within ε distance
    for i, point in enumerate(points):
        # Find neighbors within ε distance
        neighbors = tree.query_ball_point(point, epsilon)
        
        # Add edges between the point and its neighbors
        for neighbor_index in neighbors:
            if neighbor_index > i:  # Only add edges once per pair
                eps_nets.add((i, neighbor_index))
    
    return list(eps_nets)



def euclidean_distance(p1, p2):
    """Calculate the Euclidean distance between two points."""
    return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)

def get_distances(points, edges):
    """Calculate the distances of the edges."""
    distances = []
    if isinstance(edges, pd.DataFrame):
        edges = edges.values
    for edge in edges:
        point1, point2 = points.iloc[edge[0],:].values, points.iloc[edge[1],:].values
        distance = euclidean_distance(point1, point2)
        distances.append(distance)
    return distances


def filter_top(integers, percentile=95):
    """
    Filters out the top percentage of values from a list of integers based on the given percentile.
        - integers (list): List of integers to filter.
        - percentile (int): The percentile cutoff to filter the top values. Default is 95.
    """
    # Sort the list of integers
    sorted_integers = sorted(integers)
    
    # Calculate the cutoff based on the given percentile
    cutoff = np.percentile(sorted_integers, percentile)
    
    # Filter the list to include only values below the calculated cutoff
    filtered_list = [x for x in sorted_integers if x <= cutoff]
    
    return filtered_list



"""
# Compute features
"""
from scipy.stats import entropy
import itertools


def compute_abundance(df, phenotype_col = 'celltypes'):
    return df[phenotype_col].value_counts()

def compute_proportions(df, phenotype_col = 'celltypes'):
    return df[phenotype_col].value_counts(normalize=True)

def compute_ratio(df, phenotype_col = 'celltypes'):
    counts = df[phenotype_col].value_counts()
    ratios = {}
    for type1, type2 in itertools.combinations(counts.index, 2):
        # Prevent division by zero and ensure that the ratio is defined by checking count presence
        ratio_12 = counts.get(type1, 0) / counts.get(type2, 0.1)  # Use 0.1 to avoid division by zero
        ratio_21 = counts.get(type2, 0) / counts.get(type1, 0.1)  # Use 0.1 to avoid division by zero
        ratios[f"{type1}/{type2}"] = ratio_12
        ratios[f"{type2}/{type1}"] = ratio_21
    return ratios
def compute_density(df, area, phenotype_col = 'celltypes'):
    counts = df[phenotype_col].value_counts()
    return counts / area

def compute_distribution_uniformity(df, grid_size, phenotype_col = 'celltypes'):
    # Discretize the space into a grid
    df['grid_x'] = pd.cut(df['x'], bins=int(df['x'].max()/grid_size))
    df['grid_y'] = pd.cut(df['y'], bins=int(df['y'].max()/grid_size))
    grid_counts = df.groupby(['grid_x', 'grid_y', 'celltypes']).size().unstack(fill_value=0)
    # Calculate standard deviation of counts in each grid cell
    return grid_counts.std()

def compute_spatial_entropy(df, phenotype_col = 'celltypes'):
    counts = df[phenotype_col].value_counts()
    probabilities = counts / counts.sum()
    return entropy(probabilities)

def compute_area_fraction(df, total_area, phenotype_col = 'celltypes'):
    counts = df[phenotype_col].value_counts()
    # Assuming each point represents an equal unit of area
    return (counts * (1 / counts.sum()) * total_area)

def compute_geometric_features(df, phenotype_col='celltypes', area=None, grid_size=None, total_area=None):
        # Initialize an empty DataFrame
    results_df = pd.DataFrame(index=[0])

    # Compute each metric and store directly into DataFrame
    results_df['Abundance'] = [compute_abundance(df, phenotype_col).to_dict()]
    results_df['Proportions'] = [compute_proportions(df, phenotype_col).to_dict()]
    results_df['Ratios'] = [compute_ratio(df, phenotype_col)]
    if area:
        results_df['Density'] = [compute_density(df, area, phenotype_col).to_dict()]
    if grid_size:
        results_df['Uniformity'] = [compute_distribution_uniformity(df, grid_size, phenotype_col).std().to_dict()]
    results_df['Spatial Entropy'] = [compute_spatial_entropy(df, phenotype_col)]
    if total_area:
        results_df['Area Fraction'] = [compute_area_fraction(df, total_area, phenotype_col).to_dict()]

    return results_df