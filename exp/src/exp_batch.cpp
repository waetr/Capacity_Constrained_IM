#include "IMs.h"
#include "OPIM_new.h"

using namespace std;

int main(int argc, char const *argv[]) {
    //blog-catalog 10 0.1 0
    if (argc != 5) {
        printf("Usage: [dataset_name] k eps greedy_option");
        return 0;
    }
    string dataset_name = argv[1];
    string graphFilePath = "../data/" + dataset_name + ".txt";
    string ApFilePath = "../data/" + dataset_name + ".ap";
    auto args_k = (int64)atoi(argv[2]);
    auto args_eps = (double)atof(argv[3]);
    auto args_greedy = atoi(argv[4]);
    printf("Start--\n");
    double cur = clock();
    Graph G(graphFilePath);
    G.set_diffusion_model(IC); // Only support IC
    printf("read time = %.3f n=%ld m=%ld\n", time_by(cur), G.n, G.m);
    printf("Evaluating influence in [0.99,1.01]*EPT with prob. 99.9%%.\n");
    vector<bi_node> seeds;

    auto A_batch = AP_from_file(ApFilePath);

    for (auto &A : A_batch) {
        if (A != A_batch[0]) continue; // |A| == 1% of |V|
        auto k = args_k;
        auto eps = args_eps;
        printf("***************************\napsize: %zu k:%ld eps:%.2f overlap:%.3f\n\n", A.size(), k,
               eps, estimate_neighbor_overlap(G, A));
        printf("RR-OPIM+:\n\ttime = %.3f\n", method_FOPIM(G, k, A, seeds, eps, "RR+"));
        printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
        seeds.clear();
        printf("RR-OPIM:\n\ttime = %.3f\n", method_FOPIM(G, k, A, seeds, eps, "RR"));
        printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
        seeds.clear();
        printf("MG-OPIM:\n\ttime = %.3f\n", method_FOPIM(G, k, A, seeds, eps, "MG"));
        printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
        seeds.clear();
        printf("OPIM-C:\n\ttime = %.3f\n", method_local_OPIM(G, k, A, eps, seeds));
        printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
        seeds.clear();
        printf("IMM:\n\ttime = %.3f\n", method_local_IMM(G, k, A, eps, seeds));
        printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
        seeds.clear();
        printf("Degree:\n\ttime = %.3f\n", method_local_Degree(G, k, A, seeds));
        printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
        seeds.clear();
        printf("PageRank:\n\ttime = %.3f\n", method_local_PageRank(G, k, A, seeds));
        printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
        seeds.clear();
        if (args_greedy == 1) {
            printf("MG-Greedy:\n\ttime = %.3f\n", method_MG_CELF(G, k, A, seeds));
            printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
            seeds.clear();
            printf("RR-Greedy:\n\ttime = %.3f\n", method_RR_CELF(G, k, A, seeds));
            printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
            seeds.clear();
        } else if (args_greedy == 2) {
            printf("MG-Greedy:\n\ttime = %.3f\n", method_MG_vanilla(G, k, A, seeds));
            printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
            seeds.clear();
            printf("RR-Greedy:\n\ttime = %.3f\n", method_RR_vanilla(G, k, A, seeds));
            printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
            seeds.clear();
        }
        printf("***************************\n");
    }
    return 0;
}