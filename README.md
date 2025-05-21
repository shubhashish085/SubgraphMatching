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
Automorphism.h/.cpp - Breaking Automorphism<br />
backtracking.h/.cpp - Different type of traversals of the graph<br />
commandparser.h/.cpp - Parsing the command from the command line<br />
Enumeration.h/.cpp - Sequential Implementation of subgraph matching enumeration<br />
ParallelEnumeration.h/.cpp - Parallel Implementation of subgraph matching enumeration<br />
FilterVertices.h/.cpp - Filtering algorithms for pruning<br />
matchingcommand.h/.cpp - Matching commands from the command line<br />
graph.h/.cpp - Graph Data Structures<br />
PruningConstraints.h/.cpp - The cycle-based pruning constraints<br />
GeneratingFilterPlan.h/.cpp - Filtering Plan for execution<br />
types.h - Defining the types<br />
util.h/utilities.cpp/wtime.h - These are all about the utilities and time<br />
StudyPerformance.cpp - The file including the main function<br />



# File
# Query Graph
t 3 3  ----- t #vertices #edges<br />
v 0 0 2 ----- v id label_of_vertex degree<br />
v 1 0 2<br />
v 2 0 2<br />
e 0 1 ------ e first_vertex_of_edge second_vertex_of_edge<br />
e 1 2<br />
e 2 0<br />

# Data Graph
0 5 ----- first_vertex_of_edge second_vertex_of_edge<br />

# Run
mkdir build<br />
cmake ..<br />
make<br />
./SubgraphMatching.out -q /path/to/query -d /path/to/data -output /path/to/output<br />





