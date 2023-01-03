#include <bits/stdc++.h>
#include "top.h"

using namespace std;

int main(int argc, char const *argv[]) {
    double cur = clock();
    init_commandLine(argc, argv);
    Graph G(graphFilePath, DIRECTED_G);
    G.set_diffusion_model(IC, 15);
    printf("(eps = 0.1) read time = %.3f n=%ld m=%ld\n", time_by(cur), G.n, G.m);
    vector<int64> A, empty;
    A = {261,1478,1511,2204,2778,3503,3703,4432,6543,7216};

    double startTime = omp_get_wtime();
    cout << MC_simulation_p(G, A, empty) << endl;
    double endTime = omp_get_wtime();
    printf("time without OpenMP: %.3f\n", endTime - startTime);
    return 0;
}