
# Project_PC

This repository contains two small parallel programming projects used for experiments and coursework:

- `Heat_eq_MPI` — an MPI implementation of a 2D heat equation solver (C++).
- `TSP_OpenMP` — example code for a Travelling Salesman Problem (sequential notebook and an OpenMP C program).

**Contents**

- `Heat_eq_MPI/`
	- `heat_mpi.cpp` — MPI C++ implementation of the heat equation solver.
	- `Results/` — sample program outputs (e.g. `output_1p.txt`, `output_2p.txt`, ...).

- `TSP_OpenMP/`
	- `SM_PROJECT.c` — C program using OpenMP for parallelization.
	- `Pc_sequential_.ipynb` — Jupyter notebook with a sequential TSP reference/experiments.

**Requirements**

- A C/C++ compiler with OpenMP support (e.g. `gcc` / `g++`).
- An MPI implementation (e.g. OpenMPI or MPICH) for the MPI project.


**Build & Run**

Heat_eq_MPI

- Build:

```bash
cd Heat_eq_MPI
mpicxx -O2 heat_mpi.cpp -o heat_mpi
```

- Run (example with 4 processes):

```bash
mpirun -np 4 ./heat_mpi
```

Note: command-line arguments and exact runtime options depend on the implementation in `heat_mpi.cpp`. Check the source header or comments for any required input parameters.

TSP_OpenMP

- Build the OpenMP C program:

```bash
cd TSP_OpenMP
gcc -fopenmp -O2 SM_PROJECT.c -o sm_project
```

- Run:

```bash
./sm_project
```

- Open the sequential notebook:

```bash
cd TSP_OpenMP
jupyter notebook Pc_sequential_.ipynb
```

**Results**

- Example outputs for the heat solver are in `Heat_eq_MPI/Results/` (e.g. `output_1p.txt`, `output_2p.txt`, `output_4p.txt`, ...).


---

**Contributors & Supervision**

- Contributors:

	- Sarah Alhalees
    - Raya Abu Aljamal 
	- Zahrah saleh
	- Dalia Babain
	- Alanoud Almakadi

- Supervision:

	- Dr. Aisha Blfgeeh


