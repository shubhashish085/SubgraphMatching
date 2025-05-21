# Subgraph Matching Algorithm
The Problem is NP-Hard. So levering the shared memory implementation to find the exact number of matches of particular query graph in the data graph.
Specifically **Parallel (Shared Memory) Subgraph Matching Algorithm** in C++ using OpenMP.

---

# Overview

Subgraph matching checks if a smaller query graph can be mapped to a part of a larger data graph while preserving structure and optionally node/edge labels. This implementation includes:

- Filtering unpromising candidates
- Breaking Automorphism
- Optimized Set-Intersection
- Stack-based DFS strategy
- The paper : https://arifuzzaman.faculty.unlv.edu/paper/HPEC24_subgraph.pdf
---

# File Description
Automorphism.h/.cpp - Breaking Automorphism
backtracking.h/.cpp - Different type of traversals of the graph
commandparser.h/.cpp - Parsing the command from the command line
Enumeration.h/.cpp - Sequential Implementation of subgraph matching enumeration
ParallelEnumeration.h/.cpp - Parallel Implementation of subgraph matching enumeration
FilterVertices.h/.cpp - Filtering algorithms for pruning
matchingcommand.h/.cpp - Matching commands from the command line
graph.h/.cpp - Graph Data Structures
PruningConstraints.h/.cpp - The cycle-based pruning constraints
GeneratingFilterPlan.h/.cpp - Filtering Plan for execution
types.h - Defining the types
util.h/utilities.cpp/wtime.h - These are all about the utilities and time
StudyPerformance.cpp - The file including the main function



# File
# Query Graph
t 3 3  ----- t #vertices #edges
v 0 0 2 ----- v id label_of_vertex degree
v 1 0 2
v 2 0 2
e 0 1 ------ e first_vertex_of_edge second_vertex_of_edge
e 1 2
e 2 0

# Data Graph
0 5 ----- first_vertex_of_edge second_vertex_of_edge

# Run
mkdir build
cmake ..
make
./SubgraphMatching.out -q /path/to/query -d /path/to/data -output /path/to/output





