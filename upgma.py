import numpy as np

class UPGMA:
    def __init__(self, distance_matrix: np.array, labels=None):
        """
        Initialize the UPGMA class with the distance matrix and the labels of the taxa.

        Parameters:
        distance_matrix (np.array): The distance matrix between the taxa. Must be upper triangular.
        labels (list): (Optional) The labels of the taxa. Must have same order as distance matrix.
        """
        self.distance_matrix = distance_matrix
        self.distance_matrix[np.tril_indices_from(self.distance_matrix)] = np.inf

        # (label, number of taxa)
        if labels is None:
            self.labels = [(f"taxon{i+1}", 1) for i in range(distance_matrix.shape[0])]
        else:
            self.labels = [(f"{i+1}", 1) for i in range(len(labels))]

        for i in range(len(labels)):
            print(f"{self.labels[i][0]}: {labels[i]}")
        print()
        
    def __str__(self):
        return f"UPGMA(distance_matrix={self.distance_matrix}, labels={self.labels})"
    
    def __repr__(self):
        return f"UPGMA(distance_matrix={self.distance_matrix}, labels={self.labels})"

    def find_min_distance(self):
        """
        Find the minimum distance in the distance matrix and return the distance and the index of
        the minimum distance.
        
        Returns:
        min_distance (float): The minimum distance in the distance matrix.
        min_distance_index (tuple): The index of the minimum distance in the distance matrix.
        """
        min_distance = np.min(self.distance_matrix)
        min_distance_index = np.unravel_index(
            np.argmin(self.distance_matrix), self.distance_matrix.shape
        )
        return min_distance, min_distance_index

    def prep_for_new_distances(self):
        """
        Copy the upper triangle to the lower triangle to make the distance matrix symmetric and to
        make it easier to calculate the new distances using average linkage.

        Returns:
        None
        """
        self.distance_matrix = np.triu(self.distance_matrix) + np.triu(self.distance_matrix, 1).T

    def calculate_new_distances(self, min_distance_index):
        """
        Calculate the new distances using the average linkage method. The new distances are
        calculated, new_distance = (old_distance_1 * n1 + old_distance_2 * n2) / (n1 + n2), where
        n1 and n2 are the number of taxa in the two clusters and old_distance_1 and old_distance_2
        are the distances between the two clusters and the other clusters.

        Parameters:
        min_distance_index (tuple): The index of the minimum distance in the distance matrix.

        Returns:
        new_distances (list): The new distances between the new cluster and the other clusters.
        """
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
        """
        Append the new distances as a new column in the distance matrix and then append a new row
        of np.inf to the distance matrix to keep it symmetric.

        Parameters:
        new_distances (list): The new distances between the new cluster and the other clusters.

        Returns:
        None
        """
        self.distance_matrix = np.column_stack((self.distance_matrix, new_distances))
        self.distance_matrix = np.vstack(
            (self.distance_matrix, np.inf * np.ones(len(self.distance_matrix) + 1))
        )

    def create_new_label(self, min_distance_index):
        """
        Generate a new Newick label with the associated number of tax children for the new cluster.
        For example, ('((1, 4), 6)', 3) means that the new cluster is a cluster of 3 tax children
        where the 1:4 cluster is clustered with 6.

        Parameters:
        min_distance_index (tuple): The index of the minimum distance in the distance matrix.

        Returns:
        new_label (tuple): The new Newick label with the associated number of tax children.

        Note:
        The new label is generated such that cluster with most tax is the first element (left side)
        of the new label.
        """
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

    def find_larger_index(self, min_distance_index):
        """
        Find the larger index of the two indices in the min_distance_index tuple.
        
        Parameters:
        min_distance_index (tuple): The index of the minimum distance in the distance matrix.
        
        Returns:
        smaller_index (int): The smaller index of the two indices in the min_distance_index tuple.
        larger_index (int): The larger index of the two indices in the min_distance_index tuple.
        """
        smaller_index = min(min_distance_index)
        larger_index = max(min_distance_index)
        return smaller_index, larger_index

    def update_distance_matrix(self, min_distance_index, new_distances):
        """
        Update the distance matrix by appending the new calculated distances, deleting the larger 
        index row and column first, and then the smaller index row and column. Sets the lower
        triangle of the distance matrix to np.inf for finding the next minimum distance.

        Parameters:
        min_distance_index (tuple): The index of the minimum distance in the distance matrix.
        new_distances (list): The new distances between the new cluster and the other clusters.

        Returns:
        None
        """
        # find the smaller and larger index so we can delete the larger index first
        smaller_index, larger_index = self.find_larger_index(min_distance_index)
        self.append_new_distances(new_distances)

        self.distance_matrix = np.delete(self.distance_matrix, larger_index, axis=0)
        self.distance_matrix = np.delete(self.distance_matrix, larger_index, axis=1)

        self.distance_matrix = np.delete(self.distance_matrix, smaller_index, axis=0)
        self.distance_matrix = np.delete(self.distance_matrix, smaller_index, axis=1)

        self.distance_matrix[np.tril_indices_from(self.distance_matrix)] = np.inf

    def update_labels(self, min_distance_index, new_label):
        """
        Updates the labels by deleting the larger index first and then the smaller index. Appends
        the new label to the labels list.

        Parameters:
        min_distance_index (tuple): The index of the minimum distance in the distance matrix.
        new_label (tuple): The new Newick label with the associated number of tax children.

        Returns:
        None
        """
        smaller_index, larger_index = self.find_larger_index(min_distance_index)
        del self.labels[larger_index]
        del self.labels[smaller_index]
        self.labels.append(new_label)

    def compute_clusters(self, verbose=False):
        """
        Compute the clusters using the UPGMA algorithm. The algorithm finds the minimum distance
        in the distance matrix, calculates the new distances using the average linkage method,
        appends the new distances to the distance matrix, creates a new label for the new cluster,
        updates the labels, and updates the distance matrix. This process is repeated until there
        is only one cluster left.
        
        Parameters:
        verbose (bool): (Optional) Print the new label at each iteration.
        
        Returns:
        None
        """
        while (len(self.labels) > 1):
            min_distance, min_distance_index = self.find_min_distance()
            self.prep_for_new_distances()
            new_distances = self.calculate_new_distances(min_distance_index)
            new_label = self.create_new_label(min_distance_index)
            if verbose:
                print(new_label[0])
            self.update_labels(min_distance_index, new_label)
            self.update_distance_matrix(min_distance_index, new_distances)