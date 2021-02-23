// g++ count.cpp -O3 -ltbb

#include <chrono>
#include <iostream>
#include <vector>
#include "parallel_hashmap/phmap.h"
#include "parallel_hashmap/phmap_utils.h"
#include "tbb/parallel_for.h"
#include <tbb/parallel_reduce.h>
#include "tbb/parallel_sort.h"
#include "tbb/spin_mutex.h"

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>

typedef phmap::parallel_flat_hash_map<long long, int, phmap::priv::hash_default_hash<long long>, phmap::priv::hash_default_eq<long long>, phmap::priv::Allocator<phmap::priv::Pair<const long long, int>>, 4, tbb::spin_mutex> wedgeMap;
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

std::vector<int> preProcessing(graph& G, const int& nNodes) {

    // sorting in decreasing order of degrees
    std::vector<int> idx(nNodes);

    tbb::parallel_for(tbb::blocked_range<int>(0, nNodes), [&](tbb::blocked_range<int> r) {
        for (int i = r.begin(); i < r.end(); ++i){
            idx[i] = i;
        }
    });
    
    tbb::parallel_sort(idx.begin(), idx.end(), [&G](int i1, int i2) {return G[i1].size() > G[i2].size();});
    
    std::vector<int> rank(nNodes);
    
    tbb::parallel_for(tbb::blocked_range<int>(0, nNodes), [&](tbb::blocked_range<int> r) {
        for (int i = r.begin(); i < r.end(); ++i){
            rank[idx[i]] = i;
        }
    });

    tbb::parallel_for(tbb::blocked_range<int>(0, nNodes), [&](tbb::blocked_range<int> r) {
        for (int i = r.begin(); i < r.end(); ++i){
            tbb::parallel_sort(G[i].begin(), G[i].end(), [&rank](int i1, int i2) {return rank[i1] > rank[i2];});
        }
    });

    return rank;
}

std::vector<wedgeMap> getWedges(const graph& G, const std::vector<int>& rank, const int& nNodes) {

    int partitions = nNodes / 64;
    
    std::vector<wedgeMap> wedgeMapList(partitions);

    int partitionSize = nNodes / partitions;

    tbb::parallel_for(tbb::blocked_range<int>(0, partitions), [&](tbb::blocked_range<int> b) {
        for (int s = b.begin(); s < b.end(); ++s) {
            int start = partitionSize * s;
            tbb::parallel_for(tbb::blocked_range<int>(start, (s == partitions - 1) ? nNodes : start + partitionSize), [&](tbb::blocked_range<int> r) {
                for (long long u1 = r.begin(); u1 < r.end(); ++u1)
                    for (int v : G[u1])
                        if (rank[v] > rank[u1])
                            for (long long u2 : G[v])
                                if (rank[u2] > rank[u1])
                                    wedgeMapList[s][(u1 + u2) * (u1 + u2 + 1) / 2 + u2] += 1;
                                else
                                    break;
                        else
                            break;
            });
        }
    });

    return wedgeMapList;
}

int countEWedges(const std::vector<wedgeMap>& wedgeMapList) {

    int partitions = wedgeMapList.size();

    std::vector<int> wedgeCounts(partitions);
    tbb::parallel_for(tbb::blocked_range<int>(0, partitions), [&](tbb::blocked_range<int> r) {
        for (int s = r.begin(); s < r.end(); ++s) {
            wedgeMap W = wedgeMapList[s];
            int B = 0;
            for (auto const& w : W) {
                int n = w.second;
                B += n * (n - 1) / 2;
            }
            wedgeCounts[s] = B;
        }
    });

    return tbb::parallel_reduce( 
            tbb::blocked_range<int>(0, partitions),
            0,
            [&](tbb::blocked_range<int> r, int count)
            {
                for (int i = r.begin(); i < r.end(); ++i) {
                    count += wedgeCounts[i];
                }

                return count;
            }, std::plus<int>() );
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
    
    std::vector<int> rank = preProcessing(G, vLeft + vRight);

    auto start = get_time();

    std::vector<wedgeMap> W = getWedges(G, rank, vLeft + vRight);

    int B = countEWedges(W);
    
    std::cout << "Number of butterflies: " << B << "\n";

    auto finish = get_time();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
    std::cout << "Elapsed time = " << duration.count() << " ms\n";

    return 0;
}
