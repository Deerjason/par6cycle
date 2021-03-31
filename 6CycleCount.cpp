/* 
    Prints the induced 6-cycle count of simple bipartite graphs

    To run:
        g++ 6CycleCount.cpp -O3 -ltbb
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
*/

#include <chrono>
#include <iostream>
#include <vector>
#include "parallel_hashmap/phmap.h"
#include "parallel_hashmap/phmap_utils.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>

typedef std::vector<std::vector<int>> graph;

graph readGraph(const char *filename, int& nEdge, int& vLeft, int& vRight) {
    graph G;

    char *f;
    int size;
    struct stat s;
    int fd = open (filename, O_RDONLY);

    /* Get the size of the file. */
    int status = fstat (fd, & s);
    size = s.st_size;

    f = (char *) mmap (0, size, PROT_READ, MAP_PRIVATE, fd, 0);

    int headerEnd = 0;
    nEdge = 0, vLeft = 0, vRight = 0;
    bool nE = true, vL = false;
    char c = f[headerEnd];
    while (c != '\n') {
        if (isdigit(c)) {
            if (nE) {
                nEdge = nEdge * 10 + c - '0';
            }
            else if (vL) {
                vLeft = vLeft * 10 + c - '0';
            }
            else {
                vRight = vRight * 10 + c - '0';
            }
        }
        else {
            if (nE) {
                nE = false;
                vL = true;
            }
            else {
                vL = false;
            }
        }
        headerEnd += 1;
        c = f[headerEnd];
    }
    G.resize(vLeft + vRight);

    int u = 0, v = 0;
    bool left = true;
    int count = 0;
    for (int i = headerEnd + 1; i < size; i++) {
        c = f[i];

        if (isdigit(c)) {
            if (left) {
                u = u * 10 + c - '0';
            }
            else {
                v = v * 10 + c - '0';
            }
        }
        else {
            if (c == ' ') {
                left = false;
            }
            // if edge has been processed
            else if (!left) {
                left = true;
                G[u].push_back(v + vLeft);
                G[v + vLeft].push_back(u);
                u = 0; v = 0;           
            }
        }
    }

    return G;
}


/*
Proof of Correctness

We will prove this algorithm is correct by contradiction.

Say that there exists a set of nodes that is inaccurately counted as an induced 6-cycle.
Since our algorithm checks for inducedness given a set of 6 unique nodes, every node set counted must be an induced 6-cycle, contradicting our claim.

Say that there exists an induced 6-cycle that is counted twice.
Since there exists a specific ordering of wedge traversal in our algorithm, a wedge can't be processed twice.
This contradicts our claim of duplicity in induced 6-cycle counting.

Say that there exists an induced 6-cycle x = {u1, u2, u3, v1, v2, v3} that isn't counted.
Then x contains three wedges a = {u2, u3, v2}, b = {u1, u3, v3}, and c = {u2, u1, v1}, which the algorithm processes in that order.
Therefore, when c is processed, the algorithm counts the triangle of wedges a, b, and c as an induced 6-cycle, contradicting our claim.

Since all possible cases lead to contradiction, this algorithm returns the correct induced 6-cycle count.
*/

bool find (const std::vector<std::pair<int, int>>& W, const int& start, const int& end, const int& u, const int& v, int& m) {
    if (m != -1) {
        int idx = m;
        while (idx < end) {
            const std::pair<int, int>& curr = W[idx];
            if (curr.first != u)
                break;
            if (curr.second == v)
                return true;
            idx++;
        }
        idx = m - 1;
        while (idx >= start) {
            const std::pair<int, int>& curr = W[idx];
            if (curr.first != u)
                break;
            if (curr.second == v)
                return true;
            idx--;
        }
        return false;
    }

    int l = start;
    int r = end - 1;
    while (l <= r) {
        m = (l + r) / 2;

        const std::pair<int, int>& wedge = W[m];
        
        if (wedge.first < u)
            l = m + 1;
        else if (wedge.first > u)
            r = m - 1;
        else if (wedge.second == v)
            return true;
        else {
            int idx = m + 1;
            while (idx < end) {
                const std::pair<int, int>& curr = W[idx];
                if (curr.first != u)
                    break;
                if (curr.second == v)
                    return true;
                idx++;
            }
            idx = m - 1;
            while (idx >= start) {
                const std::pair<int, int>& curr = W[idx];
                if (curr.first != u)
                    break;
                if (curr.second == v)
                    return true;
                idx--;
            }
            return false;
        }
    }
    return false;
}

int getCount(const graph& G, const int& vLeft) {

    std::vector<int> partitions(vLeft);

    for (int u = 0; u < vLeft; ++u) {
        int val = 0;
        for (int v : G[u]) {
            for (int u2 : G[v]) {
                if (u2 < u) {
                    val++;
                }
            }
        }
        if (u == 0) {
            partitions[u] = val;
        }
        else {
            partitions[u] = partitions[u - 1] + val;
        }
    }

    std::vector<std::pair<int, int>> W(partitions[vLeft - 1]);

    std::vector<int> wedgeCnt(vLeft);

    int c, m, m2, idx = -1, u3 = -1;

    bool skip = false;

    int count = 0;

    for (int u1 = 0; u1 < vLeft; ++u1)
        for (int v1 : G[u1])
            for (int u2 : G[v1])
                if (u2 > u1) {
                    for (int w_idx = partitions[u2 - 1]; w_idx < partitions[u2 - 1] + wedgeCnt[u2]; w_idx++) {
                        const std::pair<int, int>& w = W[w_idx];
                        if (w.first == u3 && skip)
                            continue;
                        u3 = w.first;
                        m = -1;
                        skip = find(W, partitions[u1 - 1], partitions[u1 - 1] + wedgeCnt[u1], u3, v1, m);
                        if (!skip && !find(W, partitions[u1 - 1], partitions[u1 - 1] + wedgeCnt[u1], u3, w.second, m)) {
                            if (idx != u3) {
                                c = 0;
                                m2 = -1;
                                idx = m;
                                while (idx < partitions[u1 - 1] + wedgeCnt[u1]) {
                                    const std::pair<int, int>& curr = W[idx];
                                    if (curr.first != u3)
                                        break;
                                    if (!find(W, partitions[u2 - 1], partitions[u2 - 1] + wedgeCnt[u2], u3, curr.second, m2))
                                        c++;
                                    idx++;
                                }
                                idx = m - 1;
                                while (idx >= partitions[u1 - 1]) {
                                    const std::pair<int, int>& curr = W[idx];
                                    if (curr.first != u3)
                                        break;
                                    if (!find(W, partitions[u2 - 1], partitions[u2 - 1] + wedgeCnt[u2], u3, curr.second, m2))
                                        c++;
                                    idx--;
                                }
                                idx = u3;
                            }
                            count += c;
                        }
                            
                    }
                    W[partitions[u2 - 1] + wedgeCnt[u2]] = std::make_pair(u1, v1);
                    wedgeCnt[u2]++;
                    idx = -1;
                    u3 = -1;
                    skip = false;
                }
    return count;
}

auto get_time() {return std::chrono::high_resolution_clock::now(); }

int main(int argc, char *argv[]) {

    if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " path_to_dataset" << "\n";
        return 1;
	}

    char *filename = argv[1];

    int nEdge, vLeft, vRight;

    graph G = readGraph(filename, nEdge, vLeft, vRight);

    auto start = get_time();

    int B = getCount(G, vLeft);
    
    std::cout << "Number of induced 6 cycles: " << B << "\n";

    auto finish = get_time();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
    std::cout << "Elapsed time = " << duration.count() << " ms\n";

    return 0;
}
