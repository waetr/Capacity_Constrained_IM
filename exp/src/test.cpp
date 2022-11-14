#include <bits/stdc++.h>
#include "top.h"
using namespace std;

int main(int argc, char const *argv[]) {
    init_commandLine(argc, argv);
    Graph G(graphFilePath, DIRECTED_G);
    G.set_diffusion_model(IC, 15);
    printf("n = %ld m = %ld\n", G.n, G.m);
    for(int64 i = 0; i < G.n; i++){
        if(!G.deg_in[i] && !G.deg_out[i]) printf("%ld %ld %ld\n", G.deg_out[i], G.deg_in[i], i);
    }
}