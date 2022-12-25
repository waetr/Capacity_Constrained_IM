//
// Created by asuka on 2022/10/29.
//

#include "top.h"

using namespace std;

double average(vector<double> &a) {
    double res = 0;
    for (double e : a) res += e;
    return res / a.size();
}

double SD(vector<double> &a) {
    double res = 0;
    double mean = average(a);
    for (double e : a) res += (e - mean) * (e - mean);
    return sqrt(res / (a.size() - 1.0));
}

int main(int argc, char const *argv[]) {
    init_commandLine(argc, argv);
    double cur = clock();
    Graph G(graphFilePath, DIRECTED_G);
    G.set_diffusion_model(IC, 15);
    cout << "open graph time = " << time_by(cur) << endl;
    cout << "total prob: " << G.sum_non_single_p << endl;

    vector<int64> empty;


    vector<double> a100[10], b100[10], a1000[10], b1000[10], a10000[10], b10000[10];
    vector<int64> s100[10], s1000[10], s10000[10];
    for (int i = 0; i < 10; i++) {
        generate_ap_by_degree(G, s100[i], 10);
        generate_ap_by_degree(G, s1000[i], 100);
        generate_ap_by_degree(G, s10000[i], 1000);
    }
    for (int i = 0; i < 10; i++) {
        RRContainer R1(G, empty, true);
        R1.resize(G, 100000);
        RRContainer R2(G, empty, true);
        R2.resize_with_IIS(G, 100000);
        for (int j = 0; j < 10; j++) {
            a100[j].emplace_back(1.0 * R1.self_inf_cal(G, s100[j]) / R1.numOfRRsets() * G.n);
            b100[j].emplace_back(R2.influence_IIS(G, s100[j]));
            a1000[j].emplace_back(1.0 * R1.self_inf_cal(G, s1000[j]) / R1.numOfRRsets() * G.n);
            b1000[j].emplace_back(R2.influence_IIS(G, s1000[j]));
            a10000[j].emplace_back(1.0 * R1.self_inf_cal(G, s10000[j]) / R1.numOfRRsets() * G.n);
            b10000[j].emplace_back(R2.influence_IIS(G, s10000[j]));
        }
        printf("count:%d\n", i);
    }
    vector<double> A100, B100, A1000, B1000, A10000, B10000;
    for (int i = 0; i < 10; i++) {
        A100.emplace_back(SD(a100[i]));
        B100.emplace_back(SD(b100[i]));
        A1000.emplace_back(SD(a1000[i]));
        B1000.emplace_back(SD(b1000[i]));
        A10000.emplace_back(SD(a10000[i]));
        B10000.emplace_back(SD(b10000[i]));
    }
    printf("Size 100   | Average error of RIS: %.3f (SD: %.3f)\n", average(A100), SD(A100));
    printf("Size 100   | Average error of IIS: %.3f (SD: %.3f)\n", average(B100), SD(B100));
    printf("Size 1000  | Average error of RIS: %.3f (SD: %.3f)\n", average(A1000), SD(A1000));
    printf("Size 1000  | Average error of IIS: %.3f (SD: %.3f)\n", average(B1000), SD(B1000));
    printf("Size 10000 | Average error of RIS: %.3f (SD: %.3f)\n", average(A10000), SD(A10000));
    printf("Size 10000 | Average error of IIS: %.3f (SD: %.3f)\n", average(B10000), SD(B10000));
    return 0;
}