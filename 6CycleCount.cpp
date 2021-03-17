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
#include "tbb/parallel_sort.h"

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>

typedef phmap::parallel_flat_hash_map<long long, int> wedgeMap;
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
    G.reserve(vLeft + vRight);

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

int getCount(const graph& G, const int& vLeft) {
    std::vector<phmap::flat_hash_map<int, std::vector<int>>> W;
    W.resize(vLeft);
    
    int count = 0;

    for (int u1 = 0; u1 < vLeft; ++u1)
        for (int v1 : G[u1])
            for (int u2 : G[v1])
                if (u2 > u1) {
                    for (auto const& w : W[u1]) {
                        int u3 = w.first;
                        if (std::find(W[u1][u3].begin(), W[u1][u3].end(), v1) == W[u1][u3].end())
                            for (int v2 : w.second)
                                if (std::find(W[u2][u3].begin(), W[u2][u3].end(), v2) == W[u2][u3].end()) {
                                    int c = 0;
                                    for (int v3 : W[u2][u3])
                                        if (std::find(W[u1][u3].begin(), W[u1][u3].end(), v3) == W[u1][u3].end())
                                            c++;
                                    count += c;
                                }
                    }
                    W[u2][u1].push_back(v1);
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
