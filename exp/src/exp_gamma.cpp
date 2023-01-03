#include <bits/stdc++.h>
#include "top.h"
#include <omp.h>

using namespace std;

const int MAX_THREADS = 128;
const int exp_round = 3;

int main(int argc, char const *argv[]) {
    double start_time = omp_get_wtime();
    double cur = clock();
    init_commandLine(argc, argv);
    Graph G(graphFilePath, DIRECTED_G);
    G.set_diffusion_model(IC);
    printf("read time = %.3f n = %ld m = %ld\n", time_by(cur), G.n, G.m);

    auto A_batch = AP_from_file(ApFilePath);
    double candidate_exc_v[exp_round][G.n], one_v_[exp_round][G.n];

    for (auto &AP: A_batch) {
        set<int64> AP_set(AP.begin(), AP.end());
        set<int64> candidate_set;
        for (auto ap : AP) {
            for (auto &e: G.g[ap])
                if (AP_set.find(e.v) == AP_set.end())
                    candidate_set.insert(e.v);
        }
        vector<int64> candidate;
        candidate.assign(candidate_set.begin(), candidate_set.end());
        printf("**********AP size = %zu\n", AP.size());

#pragma omp parallel for default(none) schedule(guided) shared(G, AP, candidate, candidate_exc_v, one_v_)
        for (int j = 0; j < candidate.size(); j++) {
            int64 v = candidate[j];
            set<int64> candidate_exc_one_set(candidate.begin(), candidate.end());
            candidate_exc_one_set.erase(v);
            vector<int64> candidate_exc_one;
            candidate_exc_one.assign(candidate_exc_one_set.begin(), candidate_exc_one_set.end());
            vector<int64> one_v = {v};
            for (int i = 0; i < exp_round; i++) {
                candidate_exc_v[i][v] = FI_simulation_new(G, candidate_exc_one, AP);
                one_v_[i][v] = FI_simulation_new(G, one_v, AP);
            }
            candidate_exc_one_set.insert(v);
        }

        vector<double> gamma(exp_round, 1.0);
        for (int i = 0; i < exp_round; i++) {
            auto candidate_eval = FI_simulation_new(G, candidate, AP);
            for (auto v : candidate)
                gamma[i] = min(gamma[i], (candidate_eval - candidate_exc_v[i][v]) / one_v_[i][v]);
            gamma[i] = 1.0 - gamma[i];
        }
        printf("Gamma:%.3f (SD: %.3f) ", average(gamma), SD(gamma));
        printvec(gamma);
        printf("**********\n");
    }
    double end_time = omp_get_wtime();
    printf("time with omp = %.3f\n", end_time - start_time);
    return 0;
}