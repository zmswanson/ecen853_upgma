# %%
import numpy as np

# %%
distance_matrix = np.array([[0, 2.5, 10.44, 4.12, 11.75],
                            [2.5, 0, 12.5, 6.4, 13.93],
                            [10.44, 12.5, 0, 6.48, 1.41],
                            [4.12, 6.4, 6.48, 0, 7.35],
                            [11.75, 13.93, 1.41, 7.35, 0]])

distance_matrix[np.tril_indices_from(distance_matrix)] = np.inf

# (label, number of taxa)
labels = [("1", 1), ("2", 1), ("3", 1), ("4", 1), ("5", 1)]

# %%
def find_min_distance(distance_matrix):
    min_distance = np.min(distance_matrix)
    min_distance_index = np.unravel_index(np.argmin(distance_matrix), distance_matrix.shape)
    return min_distance, min_distance_index

# %% copy the upper triangle to the lower triangle to make it symmetric and to make it
# easier to calculate the new distances
def prep_for_new_distances(distance_matrix):
    distance_matrix = np.triu(distance_matrix) + np.triu(distance_matrix, 1).T
    return distance_matrix

# %% calculate the new distances using the average linkage
# new_distance = (old_distance_1 * n1 + old_distance_2 * n2) / (n1 + n2)
# n1 and n2 are the number of taxa in the two clusters
# old_distance_1 and old_distance_2 are the distances between the two clusters and the other clusters
def calculate_new_distances(distance_matrix, min_distance_index, labels):
    new_distances = [(distance_matrix[min_distance_index[0], i] * labels[min_distance_index[0]][1] +
                      distance_matrix[min_distance_index[1], i] * labels[min_distance_index[1]][1]) /
                     (labels[min_distance_index[0]][1] + labels[min_distance_index[1]][1])
                     for i in range(distance_matrix.shape[0])]
    return new_distances

# %% append the new distances as a new column and then append a new row of np.inf
def append_new_distances(distance_matrix, new_distances):
    distance_matrix = np.column_stack((distance_matrix, new_distances))
    distance_matrix = np.vstack((distance_matrix, np.inf * np.ones(len(distance_matrix) + 1)))
    return distance_matrix

# %% find the label with the most taxa and put it first
def create_new_label(labels, min_distance_index):
    min_taxa_idx = 0
    max_taxa_idx = 1

    if labels[min_distance_index[0]][1] < labels[min_distance_index[1]][1]:
        min_taxa_idx = 1
        max_taxa_idx = 0

    new_label = (f"({labels[min_distance_index[min_taxa_idx]][0]}, {labels[min_distance_index[max_taxa_idx]][0]})",
                 labels[min_distance_index[min_taxa_idx]][1] + labels[min_distance_index[max_taxa_idx]][1])
    return new_label

# %% find and remove the larger index first, so that the smaller index is still valid
# find the larger index
def find_larger_index(min_distance_index):
    smaller_index = min(min_distance_index)
    larger_index = max(min_distance_index)
    return smaller_index, larger_index

def update_distance_matrix(distance_matrix, min_distance_index):
    smaller_index, larger_index = find_larger_index(min_distance_index)

    distance_matrix = np.delete(distance_matrix, larger_index, axis=0)
    distance_matrix = np.delete(distance_matrix, larger_index, axis=1)

    distance_matrix = np.delete(distance_matrix, smaller_index, axis=0)
    distance_matrix = np.delete(distance_matrix, smaller_index, axis=1)

    distance_matrix[np.tril_indices_from(distance_matrix)] = np.inf
    return distance_matrix

# %%
def update_labels(labels, new_label, min_distance_index):
    smaller_index, larger_index = find_larger_index(min_distance_index)
    del labels[larger_index]
    del labels[smaller_index]
    labels.append(new_label)
    return labels


# %% 
while (len(labels) > 1):
    min_distance, min_distance_index = find_min_distance(distance_matrix)
    distance_matrix = prep_for_new_distances(distance_matrix)
    new_distances = calculate_new_distances(distance_matrix, min_distance_index, labels)
    distance_matrix = append_new_distances(distance_matrix, new_distances)
    new_label = create_new_label(labels, min_distance_index)
    print(new_label[0])
    labels = update_labels(labels, new_label, min_distance_index)
    distance_matrix = update_distance_matrix(distance_matrix, min_distance_index)
# %%
import numpy as np
from upgma import UPGMA
distance_matrix = np.array([[0, 2.5, 10.44, 4.12, 11.75, 2.5],
                            [2.5, 0, 12.5, 6.4, 13.93, 3.5],
                            [10.44, 12.5, 0, 6.48, 1.41, 4.5],
                            [4.12, 6.4, 6.48, 0, 7.35, 5.5],
                            [11.75, 13.93, 1.41, 7.35, 0, 6.5],
                            [2.5, 3.5, 4.5, 5.5, 6.5, 0]])

labels = ["p53", "mdm2", "bcl2", "cyclinE", "caspase8", "dummy"]

upgma = UPGMA(distance_matrix, labels)
upgma.compute_clusters()

# %%