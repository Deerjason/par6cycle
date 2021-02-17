// g++ count.cpp -O3 -ltbb

#include <chrono>
#include <iostream>
#include <vector>
#include <numeric>
#include "parallel_hashmap/phmap.h"
#include "parallel_hashmap/phmap_utils.h"
#include "tbb/parallel_for.h"

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>

using namespace std;

vector<vector<int>> readGraph(char *filename, int* nEdge, int* vLeft, int* vRight) {
    vector<vector<int>> G;

    char *f;
    int size;
    struct stat s;
    int fd = open (filename, O_RDONLY);

    /* Get the size of the file. */
    int status = fstat (fd, & s);
    size = s.st_size;

    f = (char *) mmap (0, size, PROT_READ, MAP_PRIVATE, fd, 0);

    int headerEnd = 0;
    *nEdge = 0, *vLeft = 0, *vRight = 0;
    bool nE = true, vL = false;
    char c = f[headerEnd];
    while (c != '\n') {
        if (isdigit(c)) {
            if (nE) {
                *nEdge = *nEdge * 10 + c - '0';
            }
            else if (vL) {
                *vLeft = *vLeft * 10 + c - '0';
            }
            else {
                *vRight = *vRight * 10 + c - '0';
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
    G.resize(*vLeft + *vRight);

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
                G[u].push_back(v + *vLeft);
                G[v + *vLeft].push_back(u);
                u = 0; v = 0;           
            }
        }
    }

    return G;
}

vector<int> preProcessing(vector<vector<int>>* G) {
    int nNodes = (*G).size();

    // sorting in decreasing order of degrees
    vector<int> idx(nNodes);
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(),
       [&G](int i1, int i2) {return (*G)[i1].size() > (*G)[i2].size();});
    
    vector<int> rank(nNodes);
    tbb::parallel_for(tbb::blocked_range<int>(0, nNodes), [&](tbb::blocked_range<int> r) {
        for (int i = r.begin(); i < r.end(); ++i){
            rank[idx[i]] = i;
        }
    });

    return rank;
}

phmap::flat_hash_map<tuple<int, int>, vector<int>> getWedges(vector<vector<int>>* G, vector<int>* rank) {
    int nNodes = (*G).size();
    
    phmap::flat_hash_map<tuple<int, int>, vector<int>> W;

    for (int u1 = 0; u1 < nNodes; ++u1){
        for (int i = 0; i < (*G)[u1].size(); ++i) {
            int v = (*G)[u1][i];
            if ((*rank)[v] > (*rank)[u1]) {
                for (int j = 0; j < (*G)[v].size(); ++j) {
                    int u2 = (*G)[v][j];
                    if ((*rank)[u2] > (*rank)[u1])
                        W[make_tuple(u1, u2)].push_back(v);
                }
            }
        }
    }
    return W;
}

int nChooseK(int n, int k){
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    int result = n;
    for(int i = 2; i <= k; ++i) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}

int countEWedges(phmap::flat_hash_map<tuple<int, int>, vector<int>> *W) {
    int B = 0;
    for (auto const& w : *W) {
        B += nChooseK(w.second.size(), 2);
    }
    return B;
}

auto get_time() {return chrono::high_resolution_clock::now(); }

int main(int argc, char *argv[]) {

    if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " path_to_dataset" << "\n";
        return 1;
	}

    auto start = get_time();

    char *filename = argv[1];

    int nEdge, vLeft, vRight;

    vector<vector<int>> G = readGraph(filename, &nEdge, &vLeft, &vRight);
    
    vector<int> rank = preProcessing(&G);

    phmap::flat_hash_map<tuple<int, int>, vector<int>> W = getWedges(&G, &rank);

    int B = countEWedges(&W);
    
    std::cout << "Number of butterflies: " << B << "\n";

    auto finish = get_time();
    auto duration = chrono::duration_cast<chrono::milliseconds>(finish-start);
    std::cout << "Elapsed time = " << duration.count() << " ms\n";

    return 0;
}
