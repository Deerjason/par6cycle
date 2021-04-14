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
#include "tbb/parallel_scan.h"
#include "tbb/parallel_sort.h"

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>

typedef std::vector<std::vector<int>> graph;

typedef std::vector<phmap::flat_hash_set<int>> edges;

edges readGraph(const char *filename, uint32_t& nEdge, uint32_t& vLeft, uint32_t& vRight, graph& newG) {

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

    graph G(vLeft + vRight);

    int minV = 0;

    int u = 0, v = 0;
    bool left = true;
    int count = 0;
    for (int i = headerEnd + 1; i < size; ++i) {
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

    // sorting in increasing order of degrees
    std::vector<int> idx(vLeft);

    tbb::parallel_for(tbb::blocked_range<int>(0, vLeft), [&](tbb::blocked_range<int> r) {
        for (int i = r.begin(); i < r.end(); ++i){
            idx[i] = i;
        }
    });
    
    tbb::parallel_sort(idx.begin(), idx.end(), [&G](int i1, int i2) {return G[i1].size() < G[i2].size();});
    
    std::vector<int> rank(vLeft);

    edges E(vLeft);
    
    tbb::parallel_for(tbb::blocked_range<int>(0, vLeft), [&](tbb::blocked_range<int> r) {
        for (int i = r.begin(); i < r.end(); ++i){
            rank[idx[i]] = i;
            E[i].reserve(G[i].size());
        }
    });

    newG.resize(vLeft + vRight);

    tbb::parallel_for(tbb::blocked_range<int>(0, vLeft + vRight), [&](tbb::blocked_range<int> r) {
        for (int u = r.begin(); u < r.end(); ++u){
            for (int v : G[u]) {
                if (u < vLeft) {
                    newG[rank[u]].push_back(v);
                    E[rank[u]].insert(v - vLeft);
                }
                else {
                    newG[u].push_back(rank[v]);
                }
            }
        }
    });

    tbb::parallel_for(tbb::blocked_range<int>(vLeft, vLeft + vRight), [&](tbb::blocked_range<int> r) {
        for (int i = r.begin(); i < r.end(); ++i){
            std::sort(newG[i].begin(), newG[i].end(), std::greater<int>());
        }
    });

    return E;
}

// finds location of a wedge with endpoint u in the vector of wedges W
int getm (const std::vector<std::pair<int, int>>& W, int start, int end, int u) {

    int l = start;
    int r = end - 1;
    while (l <= r) {
        int m = (l + r) / 2;

        int u2 = W[m].first;
        
        if (u2 < u)
            l = m + 1;
        else if (u2 > u)
            r = m - 1;
        else
            return m;
    }
    return -1;
}

struct Sum {
    uint64_t value;
    Sum() : value(0) {}
    Sum(Sum& s, tbb::split) {value = 0;}
    void operator()(const tbb::blocked_range<std::vector<uint64_t>::iterator>& r) {
        uint64_t temp = value;
        for(std::vector<uint64_t>::iterator it = r.begin(); it != r.end(); ++it) {
            temp += *it;
        }
        value = temp;
    }
    void join(Sum& rhs) {value += rhs.value;}
};

// returns number of induced 6 cycles
uint64_t getCount(const graph& G, const uint32_t nEdge, const uint32_t vLeft, const uint32_t vRight, const edges& E) {

    std::vector<int> partitions(vLeft);

    tbb::parallel_for(tbb::blocked_range<int>(0, vLeft - 1), [&](tbb::blocked_range<int> r) {
        for (int u = r.begin(); u < r.end(); ++u) {
            int val = 0;
            for (int v : G[u]) {
                for (int u2 : G[v]) {
                    if (u2 > u) {
                        ++val;
                    }
                }
            }
            partitions[u + 1] = val;
        }
    });

    tbb::parallel_scan(tbb::blocked_range<int>(0, vLeft), 0,
		[&](tbb::blocked_range<int> r, int sum, bool is_final_scan) {
			int tmp = sum;
			for (int u = r.begin(); u < r.end(); ++u) {
                tmp += partitions[u];
				if (is_final_scan) {
					partitions[u] = tmp;
				}
			}
			return tmp;
		},
		[](int a, int b) {
			return a + b;
	    }
    );

    std::vector<std::pair<int, int>> W(partitions[vLeft - 1]);

    tbb::parallel_for(tbb::blocked_range<int>(0, vLeft - 1), [&](tbb::blocked_range<int> r) {
        for (int u1 = r.begin(); u1 < r.end(); ++u1) {
            int i = 0;
            for (int v1 : G[u1])
                for (int u2 : G[v1])
                    if (u2 > u1) {
                        for (int w = partitions[u1]; w < partitions[u1] + i + 1; ++w) {
                            if (w == partitions[u1] + i) {
                                W[w] = std::make_pair(u2, v1);
                                break;
                            }
                            if (u2 < W[w].first) {
                                for (int w2 = partitions[u1] + i - 1; w2 >= w; --w2) {
                                    W[w2 + 1] = W[w2];
                                }
                                W[w] = std::make_pair(u2, v1);
                                break;
                            }
                        }
                        ++i;
                    }
                    else
                        break;
        }
    });

    std::vector<uint64_t> counts(vLeft - 1);

    tbb::parallel_for(tbb::blocked_range<int>(0, vLeft - 1), [&](tbb::blocked_range<int> r) {
        for (int u1 = r.begin(); u1 < r.end(); ++u1) {
            int c;
            uint64_t count = 0;
            // u1 -> v1 -> u2
            for (int w1_idx = partitions[u1]; w1_idx < partitions[u1 + 1]; ++w1_idx) {
                int idx = -1, u3 = -1;
                bool skip = false;
                int u2 = W[w1_idx].first;
                int v1 = W[w1_idx].second - vLeft;
                // u2 -> v2 -> u3
                for (int w2_idx = partitions[u2]; w2_idx < partitions[u2 + 1]; ++w2_idx) {
                    if (skip && u3 == W[w2_idx].first)
                        continue;
                    u3 = W[w2_idx].first;
                    int v2 = W[w2_idx].second - vLeft;
                    // u3 -> v1
                    skip = E[u3].contains(v1);
                    // u1 -> v2
                    if (!skip && !E[u1].contains(v2)) {
                        if (idx != u3) {
                            c = 0;
                            // u1 -> v3 -> u3
                            int m = getm(W, partitions[u1], partitions[u1 + 1], u3);
                            idx = m;
                            while (idx < partitions[u1 + 1]) {
                                const std::pair<int, int>& curr = W[idx];
                                if (curr.first != u3)
                                    break;
                                int v3 = curr.second - vLeft;
                                // u2 -> v3
                                if (v3 != v1 && !E[u2].contains(v3))
                                    ++c;
                                ++idx;
                            }
                            idx = m - 1;
                            while (idx >= partitions[u1 - 1]) {
                                const std::pair<int, int>& curr = W[idx];
                                if (curr.first != u3)
                                    break;
                                int v3 = curr.second - vLeft;
                                // u2 -> v3
                                if (v3 != v1 && !E[u2].contains(v3))
                                    ++c;
                                --idx;
                            }
                            idx = u3;
                        }
                        count += c;
                    }
                }
            }
            counts[u1] = count;
        }
    });

    Sum total;
    tbb::parallel_reduce(tbb::blocked_range<std::vector<uint64_t>::iterator>(counts.begin(), counts.end()), total);
    return total.value;
}

auto get_time() {return std::chrono::high_resolution_clock::now(); }

int main(int argc, char *argv[]) {

    if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " path_to_dataset" << "\n";
        return 1;
	}

    char *filename = argv[1];

    uint32_t nEdge, vLeft, vRight;

    graph G;

    edges E = readGraph(filename, nEdge, vLeft, vRight, G);

    auto start = get_time();

    uint64_t B = getCount(G, nEdge, vLeft, vRight, E);
    
    std::cout << "Number of induced 6 cycles: " << B << "\n";

    auto finish = get_time();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
    std::cout << "Elapsed time = " << duration.count() << " milliseconds\n";

    return 0;
}
