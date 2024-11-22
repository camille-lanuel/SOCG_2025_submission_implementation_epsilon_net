This repository contains all that is needed to run the code relative to the SOCG 2025 submission $\varepsilon$*-Net Algorithm Implementation on Hyperbolic Surfaces*. It also contains the complete results of the benchmarks cited in the paper.

**It can be used only for review purpose. All credit go to the (anonymous) authors of this project.**

# Folder structure

- cgal: our version of CGAL (see [CGAL installation](#cgal-installation)).

- code: our code presented in the paper
  
  - include: main code
    
    - Anchored_hyperbolic_surface_triangulation_2: class containing the data structure and the $\varepsilon$-net algorithm
    
    - AHST2_epsilon_net_combinatorics: tweak of Anchored_hyperbolic_surface_triangulation_2 with additional instructions in some methods for the combinatorics benchmarks
    
    - Hyperbolic_Dirichlet_domain_2: code to compute a Dirichlet domain from fundamental domain whose vertices are the same point on the surface
  
  - demos: two demos, used to produce the figures seen in the paper in Sections 8 & 9 (see [Execution of the demos](#execution-of-the-demos))
  
  - benchmarks: code we used to run the benchmarks presented in the paper

- results: the complete benchmarks results presented in the paper. The raw results are in the .txt files. The spreadsheet "results_summary.ods" contains a compilation of the results (for example, the computation of the average of each indicator).

# Running the code

## CGAL installation

We provide our version of CGAL in the "cgal" folder. It is based on an obsolete version of Loïc Dubois' CGAL branch (in the Hyperbolic_surface_triangulation_2 folder). If you already have CGAL installed on your computer, please use our version to compile our code as it won't work with your version. You do not need to build CGAL to run our code, as it is a header only library. However, you must have the [essential third party dependencies](https://doc.cgal.org/latest/Manual/thirdparty.html) and Qt6 installed.

## Compilation of the demos

To compile the demos:

```bash
cd path/to/CGAL_2024_submission/code/demo/build
CGAL_DIR=../../../cgal cmake ..
make
```

## Execution of the demos

### epsilon_net_demo

Usage:

```bash
./epsilon_net_demo [arg1] [arg2]
```

- `arg1` (default: `0.1`): value of $\varepsilon$

- `arg2` (default: random): seed for the surface generation (must  be an integer)

Runs the `epsilon_net` method with $\varepsilon$=`arg1` on a surface generated with seed `arg2`. Opens a window with a drawing of a lift of the obtained Delaunay triangulation, as shown in Section 8 of the paper. Prints the execution time of the algorithm, and the characteristics of the combinatorial map.

### dirichlet_and_epsilon_net_demo

Usage:

```bash
./dirichlet_and_epsilon_net_demo [arg1] [arg2]
```

- `arg1` (default: `0.1`): value of $\varepsilon$

- `arg2` (default: random): seed for the surface generation (must be an integer)

Runs the `epsilon_net` method with $\varepsilon$=`arg1` on a surface generated with seed `arg2`. Opens a window with a drawing the obtained Delaunay triangulation drawn following a Breadth First Search algorithm, and a Dirichlet domain (see Section 9 of the paper).
