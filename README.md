# par6cycle

This repository uses Intel Threading Building Blocks, which needs to be installed in order to run.

Parallel hashmap from https://github.com/greg7mdp/parallel-hashmap.

To run the parallel butterfly counting algorithm:
    g++ butterflyCount.cpp -O3 -ltbb
    ./a.out path_to_dataset
    
Dataset format:
    |E| |U| |V|
    u1 v1
    u2 v1
    u2 v2

Example dataset:
    3 2 2
    0 0
    0 1
    1 0
