#include <bits/stdc++.h>
#include "top.h"

using namespace std;

int main(int argc, char const *argv[]) {
    if (argc != 2) {
        printf("Usage: [dataset_name]");
        return 0;
    }
    string dataset_name = argv[1];
    string graphFilePath = "../data/" + dataset_name + ".txt";
    string ApFilePath = "../data/" + dataset_name + ".ap";
    freopen(ApFilePath.c_str(), "w", stdout);
    Graph G(graphFilePath);
    G.set_diffusion_model(IC);
    vector<int64> A;
    vector<int64> AP_size = {1,2,5,12};
    generate_ap(G, A, AP_size[AP_size.size() - 1]);
    shuffle(A.begin(), A.end(), mt19937(random_device()()));
    for (auto apsize: AP_size) {
        vector<int64> A0;
        for (int i = 0; i < apsize; i++) A0.emplace_back(A[i]);
        printf("%.3f ", estimate_neighbor_overlap(G,A0));
        for (int i = 0; i < A0.size(); i++) printf("%ld%c", A0[i], (i == A0.size() - 1) ? '\n' : ',');
    }
    return 0;
}