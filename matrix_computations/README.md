# Matrix Computations

This directory contains the implementation and related utilities for performing various matrix computations. These computations are designed to be efficient and useful for numerical analysis, scientific computing, and other applications requiring matrix operations.

## Features

- **Matrix Operations**: Includes basic matrix operations such as addition, subtraction, multiplication, and scalar operations.
- **Decompositions**: Implements matrix decompositions like LU, QR, and Cholesky.
- **Solvers**: Provides methods for solving systems of linear equations.
- **Eigenvalues and Eigenvectors**: Contains utilities to compute eigenvalues and eigenvectors of matrices.
- **Matrix Utilities**: Additional helper functions for matrix manipulation (e.g., transposition, inversion).

## Directory Structure

matrix_computations/  
├── core/ # Core matrix computation logic  
├── utils/ # Helper functions and utilities  
├── tests/ # Unit tests for all implemented functionalities 
└── examples/ # Example scripts showcasing how to use the module  

## Usage

1. **Checking that dependencies are installed**:

$ type cmake  
$ type git  
$ type cpp  

**Building**: The script Helper.sh will 'make clean' and also contains
a cheat sheet of the commands that I use to build

My build of Harmonic Balance using Matrix Computations is odd.  The build directories can only be in one location:

$ mkdir build_mc  
$ cp matrix_computations/Helper.sh build_mc 
$ cd build_mc/  
$ cat Helper.sh  


2. **Importing the Module**:
   Matrix computations is used in Harmonic Balance.  See HarmonicBalance/CMakeList.txt


Contributing

I do not provide my email address with this repository.  My intent is to welcome contributions from a few birds of a feather.
Steps to Contribute:

    Fork the repository.
    Clone your fork:
    bash

git clone https://github.com/DavidMDay/Nomad.git

Create a new branch:
bash

    git checkout -b feature/matrix-enhancements

    Commit your changes and push to your fork.
    Open a pull request to the main repository.

License

This project is licensed under the MIT License. See the LICENSE file for details.

Thank you for using the matrix_computations module!
