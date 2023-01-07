#include <bits/stdc++.h>
#include "top.h"

using namespace std;

int main(int argc, char const *argv[]) {
    printf("Start--");
    double cur = clock();
    init_commandLine(argc, argv);
    Graph G(graphFilePath, DIRECTED_G);
    G.set_diffusion_model(IC, 15);
    printf("read time = %.3f n=%ld m=%ld\n", time_by(cur), G.n, G.m);
    printf("Evaluating influence in [0.99,1.01]*EPT with prob. 99.9%%.\n");
    vector<bi_node> seeds;

    auto A_batch = AP_from_file(ApFilePath);
    vector<int64> k_batch = {10};

    for (auto &A : A_batch) {
        for (auto k : k_batch) {
            printf("***************************\napsize: %zu k:%ld\n\n", A.size(), k);
            printf("RR-OPIM:\n\ttime = %.3f\n", method_FOPIM(G, k, A, seeds, "RR"));
            printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
            seeds.clear();
            printf("RR-OPIM+:\n\ttime = %.3f\n", method_FOPIM(G, k, A, seeds, "RR+"));
            printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
            seeds.clear();
            printf("MG-OPIM:\n\ttime = %.3f\n", method_FOPIM(G, k, A, seeds, "MG"));
            printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
            seeds.clear();
            printf("DT-OPIM:\n\ttime = %.3f\n", method_FOPIM(G, k, A, seeds, "DT"));
            printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
            seeds.clear();
//            printf("OPIM-C:\n\ttime = %.3f\n", method_local_OPIM(G, k, A, seeds));
//            printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
//            seeds.clear();
            printf("Degree:\n\ttime = %.3f\n", method_local_Degree(G, k, A, seeds));
            printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
            seeds.clear();
            printf("PageRank:\n\ttime = %.3f\n", method_local_PageRank(G, k, A, seeds));
            printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
            seeds.clear();
            printf("***************************\n");
        }
    }
    return 0;
}