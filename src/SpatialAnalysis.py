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
        break
    # Remove axis and grid
    plt.axis('off')
    # Set labels and title
    plt.show()