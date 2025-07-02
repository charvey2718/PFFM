# PFFM

**PFFM** is a double cantilever beam phase-field model for interface fracture. The associated academic journal article is in review. Details will be added here in due course.

## Requirements

To compile this code, you must include the [Eigen 3.4.0 header-only library](https://gitlab.com/libeigen/eigen/-/releases/3.4.0).  
More information about Eigen can be found on the [Eigen website](https://eigen.tuxfamily.org/index.php?title=Main_Page).

## Build

The code has been successfully compiled on:
- **Windows 11** using **MinGW-W64 14.2.0**
- **Ubuntu 20.04.6 LTS** using **GCC 9.4.0**

You may compile the code using the provided `Makefile`:

```bash
make modeI
```
Ensure that the Eigen library path is correctly set in the Makefile or environment if it is not installed in a standard location.

**Note:**  
For best performance, the following compile flags for Eigen are recommended:  
```
-O3 -fno-math-errno -DNDEBUG -march=native
```
These flags enable high optimisation, disable unnecessary math error checking, and allow Eigen to fully exploit your hardware’s capabilities.

## Usage

Model parameters (geometry, material properties, solver settings, etc.) can be modified directly in the source code near the start of the `main` function. After changing these parameters, the code must be recompiled before running.

## License

## Citation

Once the article is published, citation info with be provided here so that the work can be properly referenced.
