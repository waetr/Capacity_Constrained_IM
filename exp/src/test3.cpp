//
// Created by asuka on 2022/10/29.
//

#include "top.h"

using namespace std;

void generate_ap_by_what(Graph &graph, std::vector<int64> &A, int64 size) {
    A.clear();
    std::vector<double> p(graph.n);
    double a = -16.1, b = -0.5486, c = 31.08;
    for (int i = 0; i < graph.n; i++) {
        p[i] = graph.deg_out[i] == 0 ? 0 : a * pow(graph.deg_out[i], b) + c;
        if (i > 0) p[i] += p[i - 1];
    }
    for (int i = 0; i < size; i++) {
        int64 v = choose_from_distributionP(p);
        while (std::find(A.begin(), A.end(), v) != A.end()) v = choose_from_distributionP(p);
        A.emplace_back(v);
    }
}

int main(int argc, char const *argv[]) {
    init_commandLine(argc, argv);
    double cur = clock();
    Graph G(graphFilePath, DIRECTED_G);
    G.set_diffusion_model(IC, 15);
    cout << "open graph time = " << time_by(cur) << endl;

    double a = 0, b = 0, c = 0;

    for (int i = 1; i <= 10; i++) {
        vector<int64> A;
        generate_ap(G, A, 100);
        a += estimate_neighbor_overlap(G, A);

        generate_ap_by_degree(G, A, 100);
        b += estimate_neighbor_overlap(G, A);

        generate_ap_by_what(G, A, 100);
        c += estimate_neighbor_overlap(G, A);
    }
    printf("%.3f %.3f %.3f\n", a / 10, b / 10, c / 10);

    return 0;
}