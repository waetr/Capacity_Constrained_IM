#include <bits/stdc++.h>
#include "top.h"

using namespace std;

double calculate1(Graph &G, vector<int64> &A) {
    int64 k = 10;
    double sum_log = 0;
    for (auto ap:A) sum_log += logcnk(G.deg_out[ap], k);
    return sum_log;
}

double evaluation(Graph &G, vector<int64> &A) {
    double sum_log = 0;
    for (auto ap:A) sum_log += 2 * G.deg_out[ap] * log(2) - 0.5 * log(3.1415926 * G.deg_out[ap]);
    return sum_log;
}

double evaluation1(Graph &G, vector<int64> &A) {
    int64 k = 10;
    int64 sum_pp = 0;
    for (auto ap:A) sum_pp += G.deg_out[ap];
    double sum_log = logcnk(sum_pp, k * A.size());
    return sum_log;
}

double evaluation2(Graph &G, vector<int64> &A) {
    int64 k = 10;
    int64 sum_pp = 0;
    for (auto ap:A) sum_pp += G.deg_out[ap];
    double sum_log = log(sum_pp) * k * A.size();
    return sum_log;
}

int main(int argc, char const *argv[]) {
    double cur = clock();
    init_commandLine(argc, argv);
    Graph G(graphFilePath, DIRECTED_G);
    G.set_diffusion_model(IC, 15);
    vector<int64> A;
    vector<int64> AP_size = {800, 1600, 3200, 8000};
    for (auto apsize: AP_size) {
        generate_ap(G, A, apsize);
        cout << "cal:" << calculate1(G, A) << endl;
        cout << "eva:" << evaluation(G, A) << endl;
        cout << "eva1:" << evaluation1(G, A) << endl;
        cout << "eva2:" << evaluation2(G, A) << endl << endl;
    }
    return 0;
}