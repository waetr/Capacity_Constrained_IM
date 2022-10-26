#include "top.h"

using namespace std;

int main(int argc, char const *argv[]) {
    init_commandLine(argc, argv);
    Graph G(graphFilePath, DIRECTED_G);
    cout << "!\n";
    G.set_diffusion_model(IC, 15);

    cout << argv[1] << endl;
    cout << "nodes: " << G.n << " edges: " << G.m << endl;
    int64 res = 0;
    for(int i = 0; i < G.n; i++) if(G.deg_out[i] < 5) res++;
    cout << "nodes that out-deg < 5: " <<  res << " | " << 100.0 * res / G.n << "%\n";
    return 0;
}