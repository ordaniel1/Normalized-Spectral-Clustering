# Normalized-Spectral-Clustering

## Introduction

Final project 
This project implements a version of the normalized spectral clustering algorithm. The algorithm takes a set of n points in Rd and performs the following steps:
1. Form the weighted adjacency matrix W from the input points.
2. Compute the normalized graph Laplacian Lnorm.
3. Determine the number of clusters k and obtain the first k eigenvectors of Lnorm.
4. Cluster the points using the K-means algorithm on the rows of the matrix containing the eigenvectors.
5. Assign each original point to a cluster based on the K-means result.


## Running the Code

### Python Interface

1. **Install Dependencies:** Make sure you have NumPy installed (`pip install numpy`).
2. **Run the Program:**
   ```bash
   python spkmeans.py <k> <goal> <filename>
-<k>: Number of clusters (if 0, use eigengap heuristic).
-<goal>: Choose from spk, wam, ddg, lnorm, or jacobi.
-<filename>: Path to the file containing N observations (.txt or .csv).

### C Implementation

1. **Compile the Code:**
      ```bash
   gcc -ansi -Wall -Wextra -Werror -pedantic-errors spkmeans.c -lm -o spkmeans
2. **Run the Program:**
   ```bash
   ./spkmeans <k> <goal> <filename>
-<k>: Number of clusters (if 0, use eigengap heuristic).
-<goal>: Choose from spk, wam, ddg, lnorm, or jacobi.
-<filename>: Path to the file containing N observations (.txt or .csv).

## Output

- For `spk`, the program outputs the calculated final centroids from the K-means algorithm.
- For other cases, the program outputs the required matrix.

## Files

1. `spkmeans.py`: Python interface.
2. `spkmeans.h`: C header file.
3. `spkmeans.c`: C implementation.
4. `spkmeansmodule.c`: Python C API wrapper.
5. `setup.py`: Setup file for building the extension.

## Building the Extension

```bash
python setup.py build_ext --inplace


