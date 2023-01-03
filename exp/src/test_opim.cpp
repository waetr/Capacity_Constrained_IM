#include <bits/stdc++.h>
#include "top.h"

using namespace std;

int main(int argc, char const *argv[]) {
    double cur = clock();
    init_commandLine(argc, argv);
    Graph G(graphFilePath, DIRECTED_G);
    G.set_diffusion_model(IC);
    printf("(eps = 0.05) read time = %.3f n=%ld m=%ld\n", time_by(cur), G.n, G.m);
    vector<int64> S, empty;

    int64 k = 50;
    OPIM_CG(G,k,0.01,1.0/G.n,S);
    //printf("inf: %.3f\n", FI_simulation_new(G, S, empty));

    return 0;
}