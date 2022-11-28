#include <bits/stdc++.h>
#include "top.h"

using namespace std;

int main(int argc, char const *argv[]) {
    double cur = clock();
    init_commandLine(argc, argv);
    Graph G(graphFilePath, DIRECTED_G);
    G.set_diffusion_model(IC, 15);
    printf("read time = %.3f\n", time_by(cur));
    cur = clock();
    vector<int64> A;
    vector<bi_node> seeds;
    generate_ap_by_degree(G, A, 100);
    int64 k = 10;

    RRContainer R(G, A, true);
    R.resize(G, 1000000);
    printf("gen time = %.3f\n", time_by(cur));


    printf("start!\n");
    RR_OPIM_Main(G, A, k, 0.05, 1.0 / G.n, seeds, false);
    printf("vanilla D-Prob quality = %.3f\n", 1.0 * R.self_inf_cal(G, seeds) / R.numOfRRsets() * G.n);
    seeds.clear();
    RR_OPIM_Main(G, A, k, 0.05, 1.0 / G.n, seeds, true);
    printf("tighten D-Prob quality = %.3f\n", 1.0 * R.self_inf_cal(G, seeds) / R.numOfRRsets() * G.n);
    seeds.clear();

    Greedy_OPIM_Main(G, A, k, 0.05, 1.0 / G.n, seeds, false);
    printf("vanilla Greedy quality = %.3f\n", 1.0 * R.self_inf_cal(G, seeds) / R.numOfRRsets() * G.n);
    seeds.clear();
    Greedy_OPIM_Main(G, A, k, 0.05, 1.0 / G.n, seeds, true);
    printf("tighten Greedy quality = %.3f\n", 1.0 * R.self_inf_cal(G, seeds) / R.numOfRRsets() * G.n);
    seeds.clear();

    cur = clock();
    method_greedy_Degree(G, k, A, seeds);
    printf("time = %.3f\nDegree quality = %.3f\n", time_by(cur),
           1.0 * R.self_inf_cal(G, seeds) / R.numOfRRsets() * G.n);
    seeds.clear();

    cur = clock();
    method_greedy_PageRank(G, k, A, seeds);
    printf("time = %.3f\nPGrank quality = %.3f\n", time_by(cur),
           1.0 * R.self_inf_cal(G, seeds) / R.numOfRRsets() * G.n);
    seeds.clear();

    cur = clock();
    method_random(G, k, A, seeds);
    printf("time = %.3f\nRandom quality = %.3f\n", time_by(cur),
           1.0 * R.self_inf_cal(G, seeds) / R.numOfRRsets() * G.n);
    seeds.clear();
    return 0;
}