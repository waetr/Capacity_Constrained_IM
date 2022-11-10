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

    double avgSeedDegree = 0;

    vector<int64> k_ = {5, 10};


    for (auto k : k_) {
        double time[5], res[5]{}, sizes[5]{};
        for (int i = 0; i < 3; i++) {
            generate_seed(G, A, 2);
 
//            cur = clock();
//            advanced_CELF_method(G, k, A, seeds, avgSeedDegree);
//            time[0] += time_by(cur) / 3.0;
//            res[0] += FI_simulation_new(G, seeds, A, 10000) / 3.0;
//            sizes[0] += seeds.size() / 3.0;
//            seeds.clear();
//
//            cur = clock();
//            Thresholding_CELF(G, k, A, seeds, 0.05, avgSeedDegree);
//            time[1] += time_by(cur) / 3.0;
//            res[1] += FI_simulation_new(G, seeds, A, 10000) / 3.0;
//            sizes[1] += seeds.size() / 3.0;
//            seeds.clear();

            cur = clock();
            dprob_CELF(G, k, A, seeds);
            time[2] += time_by(cur) / 3.0;
            res[2] += FI_simulation_new(G, seeds, A, 10000) / 3.0;
            sizes[2] += seeds.size() / 3.0;
            seeds.clear();
        }
        cout << "k = " << k << ":" << endl;
        cout << "result of maxcover-CELF:  " << res[0] << " time: " << time[0] << " size: " << sizes[0] << endl;
        cout << "result of threshold-CELF: " << res[1] << " time: " << time[1] << " size: " << sizes[1] << endl;
        cout << "result of d-prob-CELF:    " << res[2] << " time: " << time[2] << " size: " << sizes[2] << endl;
    }

    return 0;
}