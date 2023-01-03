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
    vector<int64> AP_size = {50, 100, 200, 500};
    generate_ap(G, A, AP_size[AP_size.size() - 1]);
    shuffle(A.begin(), A.end(), mt19937(random_device()()));
    for (auto apsize: AP_size) {
        printf("\nd = %ld\n",apsize);
        vector<int64> A0;
        for(int i = 0; i < apsize; i++) A0.emplace_back(A[i]);
        print_set(A0);
    }
    return 0;
}