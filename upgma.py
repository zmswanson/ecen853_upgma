import numpy as np

class UPGMA:
    def __init__(self, distance_matrix: np.array, labels):
        """
        Initialize the UPGMA class with the distance matrix and the labels of the taxa.

        Parameter
        """
        self.distance_matrix = distance_matrix
        self.distance_matrix[np.tril_indices_from(self.distance_matrix)] = np.inf

        # (label, number of taxa)
        self.labels = [(f"{i+1}", 1) for i in range(len(labels))]

        for i in range(len(labels)):
            print(f"{self.labels[i][0]}: {labels[i]}")
        print()
        
    def __str__(self):
        return f"UPGMA(distance_matrix={self.distance_matrix}, labels={self.labels})"
    
    def __repr__(self):
        return f"UPGMA(distance_matrix={self.distance_matrix}, labels={self.labels})"

    # %%
    def find_min_distance(self):
        min_distance = np.min(self.distance_matrix)
        min_distance_index = np.unravel_index(
            np.argmin(self.distance_matrix), self.distance_matrix.shape
        )
        return min_distance, min_distance_index

    # %% copy the upper triangle to the lower triangle to make it symmetric and to make it
    # easier to calculate the new distances
    def prep_for_new_distances(self):
        self.distance_matrix = np.triu(self.distance_matrix) + np.triu(self.distance_matrix, 1).T

    # %% calculate the new distances using the average linkage
    # new_distance = (old_distance_1 * n1 + old_distance_2 * n2) / (n1 + n2)
    # n1 and n2 are the number of taxa in the two clusters
    # old_distance_1 and old_distance_2 are the distances between the two clusters and the other clusters
    def calculate_new_distances(self, min_distance_index):
        new_distances = [(self.distance_matrix[min_distance_index[0], i] * 
                            self.labels[min_distance_index[0]][1] +
                          self.distance_matrix[min_distance_index[1], i] * 
                            self.labels[min_distance_index[1]][1]) /
                        (self.labels[min_distance_index[0]][1] + 
                            self.labels[min_distance_index[1]][1]) 
                        for i in range(self.distance_matrix.shape[0])]
        return new_distances

    # %% append the new distances as a new column and then append a new row of np.inf
    def append_new_distances(self, new_distances):
        self.distance_matrix = np.column_stack((self.distance_matrix, new_distances))
        self.distance_matrix = np.vstack(
            (self.distance_matrix, np.inf * np.ones(len(self.distance_matrix) + 1))
        )

    # %% find the label with the most taxa and put it first
    def create_new_label(self, min_distance_index):
        min_taxa_idx = 0
        max_taxa_idx = 1

        if self.labels[min_distance_index[0]][1] < self.labels[min_distance_index[1]][1]:
            min_taxa_idx = 1
            max_taxa_idx = 0

        new_label = (f"({self.labels[min_distance_index[min_taxa_idx]][0]}, " + 
                     f"{self.labels[min_distance_index[max_taxa_idx]][0]})",
                     self.labels[min_distance_index[min_taxa_idx]][1] +
                     self.labels[min_distance_index[max_taxa_idx]][1])
        return new_label

    # %% find and remove the larger index first, so that the smaller index is still valid
    # find the larger index
    def find_larger_index(self, min_distance_index):
        smaller_index = min(min_distance_index)
        larger_index = max(min_distance_index)
        return smaller_index, larger_index

    def update_distance_matrix(self, min_distance_index):
        smaller_index, larger_index = self.find_larger_index(min_distance_index)

        self.distance_matrix = np.delete(self.distance_matrix, larger_index, axis=0)
        self.distance_matrix = np.delete(self.distance_matrix, larger_index, axis=1)

        self.distance_matrix = np.delete(self.distance_matrix, smaller_index, axis=0)
        self.distance_matrix = np.delete(self.distance_matrix, smaller_index, axis=1)

        self.distance_matrix[np.tril_indices_from(self.distance_matrix)] = np.inf

    def update_labels(self, new_label, min_distance_index):
        smaller_index, larger_index = self.find_larger_index(min_distance_index)
        del self.labels[larger_index]
        del self.labels[smaller_index]
        self.labels.append(new_label)

    def compute_clusters(self, verbose=False):
        while (len(self.labels) > 1):
            min_distance, min_distance_index = self.find_min_distance()
            self.prep_for_new_distances()
            new_distances = self.calculate_new_distances(min_distance_index)
            self.append_new_distances(new_distances)
            new_label = self.create_new_label(min_distance_index)
            if verbose:
                print(new_label[0])
            self.update_labels(new_label, min_distance_index)
            self.update_distance_matrix(min_distance_index)