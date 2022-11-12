#include <bits/stdc++.h>
#include "top.h"
using namespace std;

int main(int argc, char const *argv[]) {
    init_commandLine(argc, argv);
    map<pair<int64,int64>, bool> fuck;
    Graph G(graphFilePath, DIRECTED_G);
    G.set_diffusion_model(IC, 15);
    for(int64 i = 0; i < G.n; i++)
        for(auto v:G.g[i]) fuck[make_pair(i, v.v)] = 1;
    for(int64 i = 0; i < G.n; i++)
        for(auto v:G.g[i]) {
            if(!fuck[make_pair(v.v, i)]) printf("%ld %ld\n", v.v, i);
        }
}