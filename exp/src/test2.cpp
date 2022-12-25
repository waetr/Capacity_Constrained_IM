#include "top.h"
#include "map"

using namespace std;

double estimate_seeds(Graph &graph, vector<int64> &ap, int64 k, vector<bi_node> &seeds) {
    map<int64, int64> h;
    for (auto u : ap)
        h[u] = min(k, graph.deg_out[u]);
    for (auto v : seeds)
        h[v.second]--;
    double res = 0;
    for (auto u : ap)
        res += 1.0 * h[u] / ap.size();
    return res;
}

int main(int argc, char const *argv[]) {
    init_commandLine(argc, argv);
    double cur = clock();
    Graph G(graphFilePath, DIRECTED_G);
    G.set_diffusion_model(IC, 15);
    cout << "open graph time = " << time_by(cur) << endl;

    vector<int64> A;
    vector<bi_node> seeds;
    vector<int64> A_size = {10};
    vector<int64> k_size = {10};
    for (auto apsize: A_size) {
        for(auto k : k_size) {
            printf("d = %ld ", apsize);
            generate_ap(G, A, apsize);
            method_local_Degree(G, k, A, seeds);
            printf("quality: %.3f size: %zu\n", FI_simulation_binode(G, seeds, A), seeds.size());
            seeds.clear();

            method_local_OPIM(G, k, A, seeds);
            printf("quality: %.3f size: %zu\n", FI_simulation_binode(G, seeds, A), seeds.size());
            seeds.clear();

            RR_OPIM_Main(G, A, k, 0.05, 1.0 / G.n, seeds, false);
            printf("quality: %.3f size: %zu\n", FI_simulation_binode(G, seeds, A), seeds.size());
            seeds.clear();
        }
    }


    stdFileOut.close();
    return 0;
}