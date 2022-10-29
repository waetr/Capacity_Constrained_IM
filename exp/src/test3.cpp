//
// Created by asuka on 2022/10/29.
//

#include "top.h"
using namespace std;

int main(int argc, char const *argv[]) {
    init_commandLine(argc, argv);
    double cur = clock();
    Graph G(graphFilePath, DIRECTED_G);
    G.set_diffusion_model(IC, 15);
    cout << "open graph time = " << time_by(cur) << endl;

    vector<int64> A, seeds;
    generate_seed(G, A, 1);
    double avgSeedDegree = 0;

    vector<int64> k_ = {1, 2, 5, 10};

    for(auto k : k_) {
        cur = clock();
        advanced_CELF_method(G, k, A, seeds, avgSeedDegree);
        cout << "time: " << time_by(cur) << endl;
        cout << "seed: " << FI_simulation_new(G, seeds, A) << endl;
        seeds.clear();

        cur = clock();
        Thresholding_CELF(G, k, A, seeds, 0.05, avgSeedDegree);
        cout << "time: " << time_by(cur) << endl;
        cout << "seed: " << FI_simulation_new(G, seeds, A) << endl;
        seeds.clear();

        cur = clock();
        dprob_CELF(G, k, A, seeds);
        cout << "time: " << time_by(cur) << endl;
        cout << "seed: " << FI_simulation_new(G, seeds, A) << endl;
        seeds.clear();
    }

    return 0;
}