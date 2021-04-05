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

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>

typedef std::vector<std::vector<int>> graph;

typedef phmap::flat_hash_map<uint64_t, bool> edges;

edges readGraph(const char *filename, uint32_t& nEdge, uint32_t& vLeft, uint32_t& vRight, graph& G) {

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

    edges E(nEdge);

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
                E[(u + v) * (u + v + 1) / 2 + v] = true;
                u = 0; v = 0;
            }
        }
    }

    return E;
}


/*
Proof of Correctness

We will prove this algorithm is correct by contradiction.

Say that there exists a set of nodes that is inaccurately counted as an induced 6-cycle.
Since our algorithm checks for inducedness given a set of 6 unique nodes, every node set counted must be an induced 6-cycle, contradicting our claim.

Say that there exists an induced 6-cycle that is counted twice.
Since there exists a specific ordering of wedge traversal in our algorithm and each wedge is stored exactly once, a wedge can't be processed twice.
Since all traversed induced 6-cycles, which is represented by the set of unique nodes {u1, u2, u3, v1, v2, v3}, have the property u1 < u2 < u3, an induced 6-cycle can't be traversed twice.
This contradicts our claim of duplicity in induced 6-cycle counting.

Say that there exists an induced 6-cycle x = {u1, u2, u3, v1, v2, v3} that isn't counted.
Then x contains three wedges a = {u1, u2, v1}, b = {u2, u3, v2}, and c = {u1, u3, v3}, which the algorithm processes in that order.
Therefore, when a is processed, the algorithm counts the triangle of wedges a, b, and c as an induced 6-cycle, contradicting our claim.

Since all possible cases lead to contradiction, this algorithm returns the correct induced 6-cycle count.
*/

// finds location of a wedge with endpoint u in the vector of wedges W
void getm (const std::vector<std::pair<int, int>>& W, const int& start, const int& end, const int& u, int& m) {

    int l = start;
    int r = end - 1;
    while (l <= r) {
        m = (l + r) / 2;

        const std::pair<int, int>& wedge = W[m];
        
        if (wedge.first < u)
            l = m + 1;
        else if (wedge.first > u)
            r = m - 1;
        else
            return;
    }
    m = -1;
    return;
}

// returns number of induced 6 cycles
int getCount(const graph& G, const uint32_t& nEdge, const uint32_t& vLeft, const uint32_t& vRight, const edges& E) {

    std::vector<int> partitions(vLeft);

    tbb::parallel_for(tbb::blocked_range<int>(0, vLeft - 1), [&](tbb::blocked_range<int> r) {
        for (int u = r.begin(); u < r.end(); ++u) {
            int val = 0;
            for (const int& v : G[u]) {
                for (const int& u2 : G[v]) {
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
		[&](const int &a, const int &b) {
			return a + b;
		});

    std::vector<std::pair<int, int>> W(partitions[vLeft - 1]);

    tbb::parallel_for(tbb::blocked_range<int>(0, vLeft - 1), [&](tbb::blocked_range<int> r) {
        for (int u1 = r.begin(); u1 < r.end(); ++u1) {
            int i = 0;
            for (const int& v1 : G[u1])
                for (const int& u2 : G[v1])
                    if (u2 > u1) {
                        for (int w = partitions[u1]; w < partitions[u1] + i + 1; ++w) {
                            if (w == partitions[u1] + i) {
                                W[w] = std::make_pair(u2, v1);
                                break;
                            }
                            int u = W[w].first;
                            int v = W[w].second;
                            if (u2 < u) {
                                for (int w2 = partitions[u1] + i - 1; w2 >= w; --w2) {
                                    W[w2 + 1] = W[w2];
                                }
                                W[w] = std::make_pair(u2, v1);
                                break;
                            }
                        }
                        ++i;
                    }
        }
    });

    std::vector<int> counts(vLeft);

    tbb::parallel_for(tbb::blocked_range<int>(0, vLeft - 1), [&](tbb::blocked_range<int> r) {
        for (int u1 = r.begin(); u1 < r.end(); ++u1) {
            int c, m, idx = -1, u3 = -1;
            bool skip = false;
            int count = 0;
            // u1 -> v1 -> u2
            for (int w1_idx = partitions[u1]; w1_idx < partitions[u1 + 1]; ++w1_idx) {
                int u2 = W[w1_idx].first;
                int v1 = W[w1_idx].second - vLeft;
                // u2 -> v2 -> u3
                for (int w2_idx = partitions[u2]; w2_idx < partitions[u2 + 1]; ++w2_idx) {
                    if (skip && u3 == W[w2_idx].first)
                        continue;
                    u3 = W[w2_idx].first;
                    int v2 = W[w2_idx].second - vLeft;
                    skip = E.contains((u3 + v1) * (u3 + v1 + 1) / 2 + v1);
                    // u3 -> v1
                    // u1 -> v2
                    if (!skip && !E.contains((u1 + v2) * (u1 + v2 + 1) / 2 + v2)) {
                        if (idx != u3) {
                            c = 0;
                            // u1 -> v3 -> u3
                            getm(W, partitions[u1], partitions[u1 + 1], u3, m);
                            idx = m;
                            while (idx < partitions[u1 + 1]) {
                                const std::pair<int, int>& curr = W[idx];
                                if (curr.first != u3)
                                    break;
                                int v3 = curr.second - vLeft;
                                // u2 -> v3
                                if (v3 != v1 && !E.contains((u2 + v3) * (u2 + v3 + 1) / 2 + v3))
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
                                if (v3 != v1 && !E.contains((u2 + v3) * (u2 + v3 + 1) / 2 + v3))
                                    ++c;
                                --idx;
                            }
                            idx = u3;
                        }
                        count += c;
                    }
                }
                idx = -1;
                skip = false;
            }
            counts[u1] = count;
        }
    });

    return tbb::parallel_reduce( 
            tbb::blocked_range<int>(0,counts.size()),
            0.0,
            [&](tbb::blocked_range<int> r, int running_total)
            {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    running_total += counts[i];
                }

                return running_total;
            }, std::plus<int>());
}

auto get_time() {return std::chrono::high_resolution_clock::now(); }

int main(int argc, char *argv[]) {

    if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " path_to_dataset" << "\n";
        return 1;
	}

    std::cout << argv[1] << std::endl;

    char *filename = argv[1];

    uint32_t nEdge, vLeft, vRight;

    graph G;

    edges E = readGraph(filename, nEdge, vLeft, vRight, G);

    auto start = get_time();

    int B = getCount(G, nEdge, vLeft, vRight, E);
    
    std::cout << "Number of induced 6 cycles: " << B << "\n";

    auto finish = get_time();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
    std::cout << "Elapsed time = " << duration.count() << " milliseconds\n";

    return 0;
}
