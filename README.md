# Triplet-Correlation-Function

This project provides a set of functions implemented in Cython for calculating triplet correlation functions, along with analysis scripts in Python. These tools are designed to analyze trajectories generated from molecular dynamics simulations.

The triplet correlation function is based on the methodology described by McNeil, Madden, Haymet and Rice:  

    --The Journal of Chemical Physics 78, 388 (1983); doi: 10.1063/1.444514
    
    With further testing by Dhabal et al.
    
    -- Phys. Chem. Chem. Phys., 2017,19, 3265-3278
    
    -- J. Chem. Phys. 7 November 2014; 141 (17): 174504. https://doi.org/10.1063/1.4898755
    
 

## Introduction

The goal of this project is to provide a convenient and efficient solution for calculating triplet correlation functions in molecular dynamics simulations. The Cython implementation allows for fast execution, while the Python analysis scripts provide flexibility and ease of use.

## Installation

To use the functions and analysis scripts provided in this repository, follow these steps:

1. Clone the repository:

git clone [https://github.com/andresordorica/Triplet-Correlation-Function/tree/main ](https://github.com/andresordorica/Triplet-Correlation-Function.git)

2. Install the required dependencies:

pip install -r requirements.txt

3. Build the Cython modules:

pip install cython 

python setup.py build_ext --inplace

4. Ensure that the required simulation trajectories are available in the appropriate format.

## Usage

The project consists of two main components:

- Cython functions: The Cython implementation provides high-performance functions for calculating the triplet correlation functions. These functions can be imported and utilized in your own Python scripts.

- Python analysis scripts: The provided Python scripts demonstrate how to utilize the Cython functions to analyze molecular dynamics trajectories. These scripts can be customized and extended to suit your specific analysis needs.

## Example

To get started, refer to the Jupyter Notebook tutorial `example.ipynb`. This tutorial provides step-by-step instructions on how to load trajectory data, calculate triplet correlation functions using the Cython functions, and visualize the results.

## Contributing

Contributions are welcome! If you have any suggestions, bug reports, or feature requests, please open an issue or submit a pull request. For major changes, please discuss them in advance.
