# Quantum-Well-Solver-CPP


# 1D & 2D Infinite Potential Well: Numerical Schrödinger Solvers

**Authors:** Yuliana Guerrero Uchima, Dayana Andrea Henao Arbelaez

## Introduction
This project provides a comprehensive numerical framework for solving both the **Time-Independent (TISE)** and **Time-Dependent (TDSE)** Schrödinger equations for an infinite potential well in one and two dimensions. 

The implementation focuses on precision and stability, utilizing finite difference methods to model quantum particle confinement within rigid boundaries.

### Time-Independent Schrödinger Equation (TISE)
For an infinite potential well, the TISE is defined as:
$$-\frac{\hbar^2}{2m} \nabla^2 \psi(\mathbf{r}) = E \psi(\mathbf{r})$$
Subject to **Dirichlet boundary conditions**: $\psi(\mathbf{r}) = 0$ at the well boundaries.

### Time-Dependent Schrödinger Equation (TDSE)
The temporal evolution of the quantum state is governed by:
$$i\hbar \frac{\partial}{\partial t} \psi(\mathbf{r}, t) = \left[ -\frac{\hbar^2}{2m} \nabla^2 + V(\mathbf{r}) \right] \psi(\mathbf{r}, t)$$

## Numerical Method: Crank-Nicolson Scheme
To ensure unconditional stability and preserve the unitarity of the time-evolution operator, we implement the **Crank-Nicolson method**.

In one dimension, the discretized evolution of the wavefunction $\Psi$ at time step $n$ is:
$$\left(I + \frac{i \Delta t}{2 \hbar} H\right) \Psi^{n+1} = \left(I - \frac{i \Delta t}{2 \hbar} H\right) \Psi^n$$
Where:
- $H$ is the discretized Hamiltonian matrix.
- $\Delta t$ is the temporal step.
- $I$ is the identity matrix.

### 2D Extension
In two dimensions, the Laplacian operator is approximated using a five-point stencil finite difference scheme:
$$\nabla^2 \Psi \approx \frac{\Psi_{i+1,j} - 2\Psi_{i,j} + \Psi_{i-1,j}}{\Delta x^2} + \frac{\Psi_{i,j+1} - 2\Psi_{i,j} + \Psi_{i,j-1}}{\Delta y^2}$$
The resulting system leads to a large, sparse matrix equation solved via numerical linear algebra techniques (e.g., LU decomposition or iterative solvers).

## Implementation Details
The solver is built in **C++** following a **Modular and Object-Oriented Programming (OOP)** approach to ensure code reusability and efficiency.

### Spatial Discretization
We discretize a rectangular domain where the potential $V$ is zero inside and infinite at the boundaries.

#### 1D Grid Setup:
```cpp
// 1D Mesh: N internal points, L well length
double dx = L / (N + 1);
for (int i = 0; i < N; i++) {
    double x = (i + 1) * dx;
    // Boundary conditions: ψ(0) = ψ(L) = 0
}
```
#### 2D Grid Setup:
```cpp
// 2D Mesh: Nx × Ny internal points
double dx = Lx / (Nx + 1);
double dy = Ly / (Ny + 1);
for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
        double x = (i + 1) * dx;
        double y = (j + 1) * dy;
        // Boundary conditions: ψ at boundaries = 0
    }
}
```
