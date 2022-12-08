#include <bits/stdc++.h>
#include "top.h"

using namespace std;

int main(int argc, char const *argv[]) {
    double cur = clock();
    init_commandLine(argc, argv);
    Graph G(graphFilePath, DIRECTED_G);
    G.set_diffusion_model(IC, 15);
    printf("(eps = 0.05) read time = %.3f\n", time_by(cur));
    vector<int64> A;
    vector<bi_node> seeds;

    int64 k = 10;

    vector<int64> A_size = {1000};
    for (auto apsize: A_size) {
        printf("d = %ld ", apsize);
        cur = clock();
        generate_ap_by_degree(G, A, apsize);
        printf("overlap: %.3f\n", estimate_neighbor_overlap(G, A));
        RRContainer R(G, A, true);
        R.resize(G, 1000000);
        printf("generate time = %.3f\n", time_by(cur));
        vector<bi_node> degree_seeds;
        method_greedy_Degree(G, k, A, degree_seeds);

        RR_OPIM_Main(G, A, k, 0.05, 1.0 / G.n, seeds, false);
        printf("vanilla D-Prob quality = %.3f size = %zu\n", 1.0 * R.self_inf_cal(G, seeds) / R.numOfRRsets() * G.n,
               seeds.size());
        ASSERT_SEED(seeds);
        seeds.clear();
        RR_OPIM_Main(G, A, k, 0.05, 1.0 / G.n, seeds, true);
        printf("tighten D-Prob quality = %.3f size = %zu\n", 1.0 * R.self_inf_cal(G, seeds) / R.numOfRRsets() * G.n,
               seeds.size());
        ASSERT_SEED(seeds);
        seeds.clear();

        Greedy_OPIM_Main(G, A, k, 0.05, 1.0 / G.n, seeds, false);
        printf("vanilla Greedy quality = %.3f size = %zu\n", 1.0 * R.self_inf_cal(G, seeds) / R.numOfRRsets() * G.n,
               seeds.size());
        ASSERT_SEED(seeds);
        seeds.clear();
        Greedy_OPIM_Main(G, A, k, 0.05, 1.0 / G.n, seeds, true);
        printf("tighten Greedy quality = %.3f size = %zu\n", 1.0 * R.self_inf_cal(G, seeds) / R.numOfRRsets() * G.n,
               seeds.size());
        ASSERT_SEED(seeds);
        seeds.clear();

        cur = clock();
        method_greedy_Degree(G, k, A, seeds);
        printf("time = %.3f\nDegree quality = %.3f size = %zu\n", time_by(cur),
               1.0 * R.self_inf_cal(G, seeds) / R.numOfRRsets() * G.n, seeds.size());
        ASSERT_SEED(seeds);
        seeds.clear();

        cur = clock();
        method_greedy_PageRank(G, k, A, seeds);
        printf("time = %.3f\nPGrank quality = %.3f size = %zu\n", time_by(cur),
               1.0 * R.self_inf_cal(G, seeds) / R.numOfRRsets() * G.n, seeds.size());
        ASSERT_SEED(seeds);
        seeds.clear();

        cur = clock();
        method_random(G, k, A, seeds);
        printf("time = %.3f\nRandom quality = %.3f size = %zu\n", time_by(cur),
               1.0 * R.self_inf_cal(G, seeds) / R.numOfRRsets() * G.n, seeds.size());
        ASSERT_SEED(seeds);
        seeds.clear();
    }
    return 0;
}