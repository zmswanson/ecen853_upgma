# ECEN-853: Computational and Systems Biology Homework 6

## Description
Python class to implement the unweighted pair-group method with arithmetic mean (UPGMA) hierarchical
clustering algorithm. The UPGMA class produces a Newick formatted representation of the clusters
given an upper-triangular distance matrix. The class is wrapped in a command-line interface (CLI)
that allows the user to select from one of two examples or to specify a file containing the distance
matrix as comma-seperated values. The clustering results are printed to the terminal at each
iteration of the algorithm.

## Installation
1. Clone the repository or copy the zipped file to your desktop.
2. If you do not have Python installed, please install the latest version from
   https://www.python.org/downloads/.
3. (Optional) If you want to work from a virtual environment, please setup and activate your virtual
   environment using the virtual environment tool of your choice. For example,
   ```shell
   # Linux and Mac
   python3 -m pip install virtualenv
   python3 -m virtualenv venv
   source ./venv/bin/activate
   
   # Windows
   python -m pip install virtualenv
   python -m virtualenv venv
   venv\Scripts\activate
   ```
4. Install the required dependencies by running the following command:
    ```shell
    pip install -r requirements.txt
    ```

## Execution
From a command terminal, execute the following with the appropriate command line arguments:
```shell
python3 main.py <arguments>
```

For help while using the CLI, please execute the following from a terminal:
```shell
python3 main.py -h
```

The command line arguments are listed below. Note that the user must specify one of the examples or
provide a distance_matrix file with the -d option to get a result. The CLI will return an error with
a reminder to select an example or provide data otherwise.

* -e, --example

    - Executes the UPGMA algorithm on a example based on data and labels from lecture notes seven.
    Doesn't require any additional input from the user.

* -e2, --example2

    - Executes the UPGMA algorithm on a example based on data and labels from 
    http://www.slimsuite.unsw.edu.au/teaching/upgma/. Doesn't require additional input from user.

* -d, --distance_matrix

    - Relative or complete path to a file containing comma-separated numbers representing the
    distance matrix. The matrix values must be upper triangular. It is not necessary to provide the
    values on the diagonal and lower triangle, as long as the correct number of commas is provided.
    See example_distance.txt for an example without values in the lower triangle and
    example2_distance.txt for an example with.

* -l, --labels

    - (Optional) Relative or complete path to a file containing labels corresponding to what each of
    the rows in the distance_matrix represent. See example_labels.txt for an example.

If a label file is provided or an example is selected, the script will print a mapping of the
numbers used to the labels. This allows for clean output in the Newick format. The script will
proceed to print the clusters generated at each iteration. The final line is the complete
hierachical cluster. As an example, the example 1 output should be:

```shell
$ python3 main.py -e
Running example 1. See UNL ECEN-853 lecture notes 7 for more details on the data and expected results.

1: p53
2: mdm2
3: bcl2
4: cyclinE
5: caspase8

(3, 5)
(1, 2)
((1, 2), 4)
(((1, 2), 4), (3, 5))
```

## License
This project is licensed under the [MIT License](https://opensource.org/licenses/MIT). You are free
to use, modify, and distribute the code in this repository for any purpose, including commercial
applications, as long as you include the original copyright notice and the license terms.