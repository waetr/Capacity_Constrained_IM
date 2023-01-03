#include "top.h"

using namespace std;
int exp_round = 3;

double average(vector<double> &a) {
    if (a.size() != exp_round) return -1;
    double res = 0;
    for (double e : a) res += e;
    return res / a.size();
}

double SD(vector<double> &a) {
    if (a.size() != exp_round) return -1;
    double res = 0;
    double mean = average(a);
    for (double e : a) res += (e - mean) * (e - mean);
    return sqrt(res / (a.size() - 1.0));
}

void printvec(vector<double> &a) {
    printf("[");
    for (double e : a) printf("%.3f,", e);
    printf("]\n");
}

int main(int argc, char const *argv[]) {
    init_commandLine(argc, argv);
    double init_cur = clock();
    Graph G(graphFilePath, DIRECTED_G);
    G.set_diffusion_model(IC, 15);
    vector<int64> A[3];
    vector<bi_node> seeds;
    printf("open graph time = %.3f n = %ld m = %ld\n", time_by(init_cur), G.n, G.m);
    vector<int64> apsize_ = {5}, k_ = {5};

    for (auto apsize : apsize_) {
        for (auto k : k_) {
            printf("**********\nd = %ld, k = %ld\n", apsize, k);
            vector<double> times, res, sizes, OP;
            for (int i = 0; i < exp_round; i++) {
                generate_ap_by_degree(G, A[i], apsize);
                OP.emplace_back(estimate_neighbor_overlap(G, A[i]));
            }
            printf("Overlap ratio:%.3f (SD: %.3f) ", average(OP), SD(OP));
            printvec(OP);
            printf("**********\n");
//            times.clear(), res.clear(), sizes.clear();
//            for (int i = 0; i < exp_round; i++) {
//                auto x = method_DProb_CELF(G, k, A[i], seeds);
//                if (x < 0) break;
//                times.emplace_back(x);
//                res.emplace_back(FI_simulation_binode(G, seeds, A[i], 10000));
//                sizes.emplace_back(seeds.size());
//                seeds.clear();
//            }
//            printf("D-Prob CELF:\n\toverall spread: %.3f (SD: %.3f) ", average(res), SD(res));
//            printvec(res);
//            printf("\tsize: %.3f (SD: %.3f) ", average(sizes), SD(sizes));
//            printvec(sizes);
//            printf("\ttime: %.3f (SD: %.3f) ", average(times), SD(times));
//            printvec(times);

//            times.clear(), res.clear(), sizes.clear();
//            for (int i = 0; i < exp_round; i++) {
//                auto x = method_local_Degree(G, k, A[i], seeds);
//                if (x < 0) break;
//                times.emplace_back(x);
//                res.emplace_back(FI_simulation_binode(G, seeds, A[i], 10000));
//                sizes.emplace_back(seeds.size());
//                seeds.clear();
//            }
//            printf("degree:\n\toverall spread: %.3f (SD: %.3f) ", average(res), SD(res));
//            printvec(res);
//            printf("\tsize: %.3f (SD: %.3f) ", average(sizes), SD(sizes));
//            printvec(sizes);
//            printf("\ttime: %.3f (SD: %.3f) ", average(times), SD(times));
//            printvec(times);

//            times.clear(), res.clear(), sizes.clear();
//            for (int i = 0; i < exp_round; i++) {
//                auto x = method_local_PageRank(G, k, A[i], seeds);
//                if (x < 0) break;
//                times.emplace_back(x);
//                res.emplace_back(FI_simulation_binode(G, seeds, A[i], 10000));
//                sizes.emplace_back(seeds.size());
//                seeds.clear();
//            }
//            printf("PageRank:\n\toverall spread: %.3f (SD: %.3f) ", average(res), SD(res));
//            printvec(res);
//            printf("\tsize: %.3f (SD: %.3f) ", average(sizes), SD(sizes));
//            printvec(sizes);
//            printf("\ttime: %.3f (SD: %.3f) ", average(times), SD(times));
//            printvec(times);

//            times.clear(), res.clear(), sizes.clear();
//            for (int i = 0; i < exp_round; i++) {
//                auto x = method_local_CELF(G, k, A[i], seeds);
//                if (x < 0) break;
//                times.emplace_back(x);
//                res.emplace_back(FI_simulation_binode(G, seeds, A[i], 10000));
//                sizes.emplace_back(seeds.size());
//                seeds.clear();
//            }
//            printf("local CELF:\n\toverall spread: %.3f (SD: %.3f) ", average(res), SD(res));
//            printvec(res);
//            printf("\tsize: %.3f (SD: %.3f) ", average(sizes), SD(sizes));
//            printvec(sizes);
//            printf("\ttime: %.3f (SD: %.3f) ", average(times), SD(times));
//            printvec(times);

//            times.clear(), res.clear(), sizes.clear();
//            for (int i = 0; i < exp_round; i++) {
//                auto x = method_greedy_vanilla(G, k, A[i], seeds);
//                if (x < 0) break;
//                times.emplace_back(x);
//                res.emplace_back(FI_simulation_binode(G, seeds, A[i], 10000));
//                sizes.emplace_back(seeds.size());
//                seeds.clear();
//            }
//            printf("greedy CELF:\n\toverall spread: %.3f (SD: %.3f) ", average(res), SD(res));
//            printvec(res);
//            printf("\tsize: %.3f (SD: %.3f) ", average(sizes), SD(sizes));
//            printvec(sizes);
//            printf("\ttime: %.3f (SD: %.3f) ", average(times), SD(times));
//            printvec(times);

            times.clear(), res.clear(), sizes.clear();
            for (int i = 0; i < exp_round; i++) {
                auto x = method_Threshold_vanilla(G, k, A[i], seeds, 0.5);
                if (x < 0) break;
                times.emplace_back(x);
                res.emplace_back(FI_simulation_binode(G, seeds, A[i], 10000));
                sizes.emplace_back(seeds.size());
                seeds.clear();
            }
            printf("Threshold CELF 0.5:\n\toverall spread: %.3f (SD: %.3f) ", average(res), SD(res));
            printvec(res);
            printf("\tsize: %.3f (SD: %.3f) ", average(sizes), SD(sizes));
            printvec(sizes);
            printf("\ttime: %.3f (SD: %.3f) ", average(times), SD(times));
            printvec(times);

            times.clear(), res.clear(), sizes.clear();
            for (int i = 0; i < exp_round; i++) {
                auto x = method_Threshold_vanilla(G, k, A[i], seeds, 0.2);
                if (x < 0) break;
                times.emplace_back(x);
                res.emplace_back(FI_simulation_binode(G, seeds, A[i], 10000));
                sizes.emplace_back(seeds.size());
                seeds.clear();
            }
            printf("Threshold CELF 0.2:\n\toverall spread: %.3f (SD: %.3f) ", average(res), SD(res));
            printvec(res);
            printf("\tsize: %.3f (SD: %.3f) ", average(sizes), SD(sizes));
            printvec(sizes);
            printf("\ttime: %.3f (SD: %.3f) ", average(times), SD(times));
            printvec(times);

            times.clear(), res.clear(), sizes.clear();
            for (int i = 0; i < exp_round; i++) {
                auto x = method_Threshold_vanilla(G, k, A[i], seeds, 0.1);
                if (x < 0) break;
                times.emplace_back(x);
                res.emplace_back(FI_simulation_binode(G, seeds, A[i], 10000));
                sizes.emplace_back(seeds.size());
                seeds.clear();
            }
            printf("Threshold CELF 0.1:\n\toverall spread: %.3f (SD: %.3f) ", average(res), SD(res));
            printvec(res);
            printf("\tsize: %.3f (SD: %.3f) ", average(sizes), SD(sizes));
            printvec(sizes);
            printf("\ttime: %.3f (SD: %.3f) ", average(times), SD(times));
            printvec(times);

            times.clear(), res.clear(), sizes.clear();
            for (int i = 0; i < exp_round; i++) {
                auto x = method_Threshold_vanilla(G, k, A[i], seeds, 0.05);
                if (x < 0) break;
                times.emplace_back(x);
                res.emplace_back(FI_simulation_binode(G, seeds, A[i], 10000));
                sizes.emplace_back(seeds.size());
                seeds.clear();
            }
            printf("Threshold CELF 0.04:\n\toverall spread: %.3f (SD: %.3f) ", average(res), SD(res));
            printvec(res);
            printf("\tsize: %.3f (SD: %.3f) ", average(sizes), SD(sizes));
            printvec(sizes);
            printf("\ttime: %.3f (SD: %.3f) ", average(times), SD(times));
            printvec(times);

            times.clear(), res.clear(), sizes.clear();
            for (int i = 0; i < exp_round; i++) {
                auto x = method_Threshold_vanilla(G, k, A[i], seeds, 0.02);
                if (x < 0) break;
                times.emplace_back(x);
                res.emplace_back(FI_simulation_binode(G, seeds, A[i], 10000));
                sizes.emplace_back(seeds.size());
                seeds.clear();
            }
            printf("Threshold CELF 0.02:\n\toverall spread: %.3f (SD: %.3f) ", average(res), SD(res));
            printvec(res);
            printf("\tsize: %.3f (SD: %.3f) ", average(sizes), SD(sizes));
            printvec(sizes);
            printf("\ttime: %.3f (SD: %.3f) ", average(times), SD(times));
            printvec(times);
        }
    }
    return 0;
}