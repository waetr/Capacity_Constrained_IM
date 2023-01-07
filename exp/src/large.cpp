#include <bits/stdc++.h>
#include "top.h"

using namespace std;

int main(int argc, char const *argv[]) {
    double cur = clock();
    init_commandLine(argc, argv);
    Graph G(graphFilePath, DIRECTED_G);
    G.set_diffusion_model(IC, 15);
    printf("(eps = 0.1) read time = %.3f n=%ld m=%ld\n", time_by(cur), G.n, G.m);
    vector<bi_node> seeds;

    vector<int64> A_size = {15000, 30000, 60000, 150000};

    int64 k = 10;

    for (auto ap_size : A_size) {
        vector<int64> A;
        generate_ap(G, A, ap_size);
        printf("apsize: %zu overlap: %.3f\n", A.size(), estimate_neighbor_overlap(G, A));


        double start0 = clock();
        printf("time = %.3f\n", method_RR_FOPIM(G, k, A, seeds, 0.1, 1.0 / G.n, false));
        printf("vanilla D-Prob quality = %.3f size = %zu\n", effic_inf(G, seeds, A),
               seeds.size());
        ASSERT_SEED(seeds);
        seeds.clear();
        printf("time = %.3f\n", method_RR_FOPIM(G, k, A, seeds, 0.1, 1.0 / G.n, true));
        printf("tighten D-Prob quality = %.3f size = %zu\n", effic_inf(G, seeds, A),
               seeds.size());
        ASSERT_SEED(seeds);
        seeds.clear();

        printf("time = %.3f\n", method_MG_FOPIM(G, k, A, seeds, 0.1, 1.0 / G.n, false));
        printf("vanilla Greedy quality = %.3f size = %zu\n", effic_inf(G, seeds, A),
               seeds.size());
        ASSERT_SEED(seeds);
        seeds.clear();
        printf("time = %.3f\n", method_DT_FOPIM(G, k, A, seeds, 0.1, 1.0 / G.n, 0.05));
        printf("Threshold quality = %.3f size = %zu\n", effic_inf(G, seeds, A),
               seeds.size());
        ASSERT_SEED(seeds);
        seeds.clear();

        cur = clock();
        method_local_Degree(G, k, A, seeds);
        printf("time = %.3f\nDegree quality = %.3f size = %zu\n", time_by(cur),
               effic_inf(G, seeds, A), seeds.size());
        ASSERT_SEED(seeds);
        seeds.clear();

        cur = clock();
        method_local_PageRank(G, k, A, seeds);

        printf("time = %.3f\nPGrank quality = %.3f size = %zu\n", time_by(cur),
               effic_inf(G, seeds, A), seeds.size());
        ASSERT_SEED(seeds);
        seeds.clear();

    }
    return 0;
}