import argparse
import os
import sys
from upgma import UPGMA
import numpy as np

def main():
    parser = argparse.ArgumentParser(
        description="Command line interface for testing the UPMA algorithm. Assumes upper-" \
                    "triangular distance matrix. If the distance matrix is not upper-triangular, " \
                    "the lower triangle will be ignored and unexpected results may occur. The " \
                    "distance matrix file must be comma-separated numbers. The optional labels "\
                    "file must be comma-separated strings. Relative paths (to the script) will " \
                    "be used for these if the absolute path is not provided. An example is " \
                    "provided to test the algorithm based on the LN7 example if you use -e.")
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-e", "--example", action='store_true', 
                        help="Select to run the provided example (LN7 example).")
    group.add_argument("-e2", "--example2", action='store_true',
                        help="Select to run the second provided example (borrowed from " \
                            "http://www.slimsuite.unsw.edu.au/teaching/upgma/).")
    group.add_argument("-d", "--distance_matrix", type=str,
                         help="The path to the distance matrix as a file")
    parser.add_argument("-l", "--labels", type=str, help="(Optional) The path to the labels as a file")
    # parser.add_argument("-o", "--output", type=str,
    #                     help="(Optional) The output file if you want to save the results. " \
    #                          "Otherwise, the results will be printed to the console.")
    args = parser.parse_args()

    distance_matrix = None
    dm_path = ""
    labels = None
    labels_path = ""

    if args.example:
        print("Running example 1. See UNL ECEN-853 lecture notes 7 for more details on the data "\
              "and expected results.")
        print()

        distance_matrix = np.array([[0, 2.5, 10.44, 4.12, 11.75],
                                    [2.5, 0, 12.5, 6.4, 13.93],
                                    [10.44, 12.5, 0, 6.48, 1.41],
                                    [4.12, 6.4, 6.48, 0, 7.35],
                                    [11.75, 13.93, 1.41, 7.35, 0]])
        
        labels = ["p53", "mdm2", "bcl2", "cyclinE", "caspase8"]
    elif args.example2:
        print("Running example 2. See http://www.slimsuite.unsw.edu.au/teaching/upgma/ for more "\
              "details on the data and expected results.")
        print()
        distance_matrix = np.array([[0, 19, 27, 8, 33, 18, 13],
                                    [19, 0, 31, 18, 36, 1, 13],
                                    [27, 31, 0, 26, 41, 32, 29],
                                    [8, 18, 26, 0, 31, 17, 14],
                                    [33, 36, 41, 31, 0, 35, 28],
                                    [18, 1, 32, 17, 35, 0, 12],
                                    [13, 13, 29, 14, 28, 12, 0]], dtype=float)
        labels = ["A-turtle","B-man","C-tuna","D-chicken","E-moth","F-monkey","G-dog"]

    else:
        if not args.distance_matrix:
            print("Please provide a distance matrix file.")
            sys.exit(1)
        else:
            dm_path = args.distance_matrix
            if not os.path.exists(dm_path):
                print(f"File {dm_path} does not exist.")
                sys.exit(1)
            distance_matrix = np.genfromtxt(dm_path, delimiter=',', dtype=float, filling_values=0)
        
        labels = [f"taxon{i+1}" for i in range(distance_matrix.shape[0])]

        if args.labels:
            labels_path = args.labels
            if not os.path.exists(labels_path):
                print(f"File {labels_path} does not exist. Assuming generic labels.")
            else:
                tmp_labels = np.loadtxt(labels_path, dtype=str, delimiter=',')
                
                if tmp_labels.shape[0] != distance_matrix.shape[0]:
                    print("The number of labels does not match the number of taxa in the "\
                          "distance matrix. Assuming generic labels.")
                else:
                    labels = tmp_labels.tolist()

    upgma = UPGMA(distance_matrix, labels)
    upgma.compute_clusters(verbose=True)

if __name__ == "__main__":
    main()