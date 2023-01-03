#include <bits/stdc++.h>
#include "top.h"

using namespace std;

int main(int argc, char const *argv[]) {
    double cur = clock();
    init_commandLine(argc, argv);
    Graph G(graphFilePath, DIRECTED_G);
    G.set_diffusion_model(IC, 15);
    printf("(eps = 0.1) read time = %.3f n=%ld m=%ld\n", time_by(cur), G.n, G.m);
    vector<int64> A;
    vector<bi_node> seeds;

    int64 k = 10;

    vector<int64> AP_size = {50,100,200,500};
    for (auto apsize: AP_size) {
        generate_ap(G, A, apsize);
        printf("overlap: %.3f\n", estimate_neighbor_overlap(G, A));

        RR_OPIM_Main(G, A, k, 0.1, 1.0 / G.n, seeds, false);
        printf("vanilla D-Prob quality = %.3f size = %zu\n", effic_inf(G, seeds, A),
               seeds.size());
        ASSERT_SEED(seeds);
        seeds.clear();
        RR_OPIM_Main(G, A, k, 0.1, 1.0 / G.n, seeds, true);
        printf("tighten D-Prob quality = %.3f size = %zu\n", effic_inf(G, seeds, A),
               seeds.size());
        ASSERT_SEED(seeds);
        seeds.clear();

        Greedy_OPIM_Main(G, A, k, 0.1, 1.0 / G.n, seeds, false);
        printf("vanilla Greedy quality = %.3f size = %zu\n", effic_inf(G, seeds, A),
               seeds.size());
        ASSERT_SEED(seeds);
        seeds.clear();
        Greedy_OPIM_Main(G, A, k, 0.1, 1.0 / G.n, seeds, true);
        printf("tighten Greedy quality = %.3f size = %zu\n", effic_inf(G, seeds, A),
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

//        cur = clock();
//        method_random(G, k, A, seeds);
//        printf("time = %.3f\nRandom quality = %.3f size = %zu\n", time_by(cur),
//               effic_inf(G, seeds, A), seeds.size());
//        ASSERT_SEED(seeds);
//        seeds.clear();
    }
    return 0;
}