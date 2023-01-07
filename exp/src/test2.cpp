#include <bits/stdc++.h>
#include "top.h"
#include <omp.h>

using namespace std;


int main(int argc, char const *argv[]) {
    init_commandLine(argc, argv);
    Graph G(graphFilePath, DIRECTED_G);
    G.set_diffusion_model(IC);
    vector<bi_node> seeds;

    auto A_batch = AP_from_file(ApFilePath);
    int k = 10;

    for (auto &A : A_batch) {
        printf("AP size :%zu\n", A.size());
        printf("\ttime0: %.3f\n", method_random(G, k, A, seeds));
        printf("\tseed0 size:%zu\n", seeds.size());
        seeds.clear();
    }
    return 0;
}