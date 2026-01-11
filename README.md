# FEM Study

Small C++ project implementing a 2D Galerkin finite element method for Poisson-type problems. Intended for experimenting with assembly strategies, matrix storage formats, and linear solvers.

## Features

- Linear P1 triangular elements (stiffness matrix and load vector assembly)
- Dirichlet boundary conditions
- Gmsh `.msh` input, VTK output
- Multiple matrix backends: dense, banded, stack-based, sparse
- Linear solvers:
  - Direct: custom LDLᵀ, Eigen (SimplicialLDLT, SparseLU)
  - Iterative: Conjugate Gradient (CPU and OpenCL)

## Code layout

- `src/FEM*.cpp` — FEM assembly and problem setup
- `src/MatrixImpl/` — matrix formats and solvers
- `src/time_series.cpp` — simple benchmarking

## Build notes

Solver backend is selected at compile time via:
`FEM_EIGEN_IMPL`, `FEM_CONJ_GRAD_IMPL`, or `FEM_NAIVE_IMPL`.
