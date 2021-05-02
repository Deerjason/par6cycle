# par6cycle

This repository uses Intel Threading Building Blocks, which needs to be installed in order to run.

Parallel hashmap from https://github.com/greg7mdp/parallel-hashmap.

To run the parallel 6-cycle counting algorithm:
<br />
g++ 6CycleCount.cpp -O3 -ltbb
<br />
./a.out path_to_dataset
    
Dataset format:
<br />
|E| |U| |V|
<br />
u1 v1
<br />
u2 v1
<br />
u2 v2

Example dataset:
<br />
3 2 2
<br />
0 0
<br />
0 1
<br />
1 0
