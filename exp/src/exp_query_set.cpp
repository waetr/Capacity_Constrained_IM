#include <bits/stdc++.h>
#include "top.h"

using namespace std;

int main(int argc, char const *argv[]) {
    double cur = clock();
    freopen("../data/orkut.ap", "w", stdout);
    init_commandLine(argc, argv);
    Graph G(graphFilePath, DIRECTED_G);
    G.set_diffusion_model(IC, 15);
    vector<int64> A;
    vector<int64> AP_size = {15000, 30000, 60000, 150000};
    generate_ap(G, A, AP_size[AP_size.size() - 1]);
    shuffle(A.begin(), A.end(), mt19937(random_device()()));
    for (auto apsize: AP_size) {
        vector<int64> A0;
        for (int i = 0; i < apsize; i++) A0.emplace_back(A[i]);
        for (int i = 0; i < A0.size(); i++) printf("%ld%c", A0[i], (i == A0.size() - 1) ? '\n' : ',');
    }
    return 0;
}