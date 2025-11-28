# PFFM

- **`Dcb2D-ModeI.cpp` and `Dcb2D-ModeI.h`** contain C++ code implementing a finite-element method (FEM) model of a double cantilever beam, which uses the phase-field fracture model for interface fracture propagation.

- The diffuse phase-field damage interacts with the surrounding bulk material, which artificially increases the apparent interface fracture toughness (assuming the bulk material has a higher fracture toughness). This effect can be mitigated by using an *effective fracture toughness* value. 

- **`BVPSolverForGceff.m`** is a MATLAB script that pre-computes the effective fracture toughness for PFFM simulations using various methods.

The associated academic journal article is available [here](https://doi.org/10.1016/j.engfracmech.2025.111546).


## Requirements

To compile the C++ code, you must include the [Eigen 3.4.0 header-only C++ library](https://gitlab.com/libeigen/eigen/-/releases/3.4.0), which provides linear algebra functionalities.
More information about Eigen can be found on the [Eigen website](https://eigen.tuxfamily.org/index.php?title=Main_Page).

## Build

The C++ code has been successfully compiled on:
- **Windows 11** using **MinGW-W64 14.2.0**
- **Ubuntu 20.04.6 LTS** using **GCC 9.4.0**

You may compile the C++ code using the provided `Makefile`:

```bash
make modeI
```
Ensure that the Eigen library path is correctly set in the `Makefile` or environment if it is not installed in a standard location.

**Note:**  
For best performance, the following compile flags for Eigen are recommended:  
```
-O3 -fno-math-errno -DNDEBUG -march=native
```
These flags enable high optimisation, disable unnecessary math error checking, and allow Eigen to fully exploit your hardware’s capabilities.

## Usage

Model parameters (geometry, material properties, solver settings, etc.) can be modified directly in the source code near the start of the `main` function. After changing these parameters, the code must be recompiled before running.

## License

This project is licensed under the [MIT License](./LICENSE).  
See the [LICENSE](./LICENSE) file for details.

## Citation

C. A. Fear, S. Wang, C. M. Harvey (2025), [Effective fracture toughness in phase-field models for interface fracture](https://doi.org/10.1016/j.engfracmech.2025.111546), Engineering Fracture Mechanics (328), 111546.
