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

    stdFileOut.open("../output/result_.out");


    stdFileOut << G.n << " " << G.m << endl;

    vector<int64> A;
    vector<bi_node> seeds;

    vector<int64> apsize_ = {50, 100, 200}, k_ = {10}, rrsize_ = {200000, 1000000};

    int count = 0;
    for (auto rrsize : rrsize_) {
        for (auto k : k_) {
            for (auto apsize : apsize_) {
                double apAvgDegree = 0, overlap_ratio = 0, RRGenTime = 0, RRGenSize = 0;
                double sRatio[5]{}, res[5]{}, sizes[5]{};
                for (int i = 0; i < 3; i++) {
                    cout << "count: " << ++count << "\n";
                    generate_ap_by_degree(G, A, apsize);
                    CandidateNeigh candidate_source(G, A, k);
                    CandidateNeigh candidate;

                    for (auto u : A) apAvgDegree += (double) G.deg_out[u] / A.size() / 3.0;
                    overlap_ratio += estimate_neighbor_overlap(G, A) / 3.0;
                    RRContainer R(G, A, true);
                    cur = clock();
                    R.resize(G, rrsize);
                    RRGenTime += time_by(cur) / 3.0, RRGenSize += R.sizeOfRRsets() / 3.0;

                    candidate.assign(G.n, candidate_source);
                    auto x = 1.0 * IMMSelection(G, A, k, seeds, candidate, R) * G.n / R.numOfRRsets();
                    sRatio[0] += estimate_seeds(G, A, k, seeds) / 3.0;
                    res[0] += x / 3.0, sizes[0] += seeds.size() / 3.0;
                    seeds.clear();

                    candidate.assign(G.n, candidate_source);
                    x = 1.0 * ThresholdSelection(G, A, k, seeds, candidate, 0.05, R) * G.n / R.numOfRRsets();
                    sRatio[1] += estimate_seeds(G, A, k, seeds) / 3.0;
                    res[1] += x / 3.0, sizes[1] += seeds.size() / 3.0;
                    seeds.clear();

                    x = 1.0 * prob_determined(G, A, k, seeds, R) * G.n / R.numOfRRsets();
                    sRatio[2] += estimate_seeds(G, A, k, seeds) / 3.0;
                    res[2] += x / 3.0, sizes[2] += seeds.size() / 3.0;
                    seeds.clear();

//                    x = 1.0 * prob(G, A, k, seeds, 1.0, R) * G.n / R.numOfRRsets();
//                    sRatio[3] += estimate_seeds(G, A, k, seeds) / 3.0;
//                    res[3] += x / 3.0, sizes[3] += seeds.size() / 3.0;
//                    seeds.clear();
                }

                stdFileOut << "apsize = " << apsize << " k = " << k << " rrsize = " << rrsize << std::endl;

                stdFileOut << "overlap ratio = " << overlap_ratio << endl;
                stdFileOut << "ap avg degree = " << apAvgDegree << endl;

                stdFileOut << "RR set generation time = " << RRGenTime << " total size = " << RRGenSize
                           << endl << endl;

                stdFileOut << "IMM             = " << fixed << setprecision(2) << res[0] << "\tsize = " << fixed
                           << setprecision(2) << sizes[0];
                stdFileOut << "\tsquander ratio = " << sRatio[0] << endl;

                stdFileOut << "Threshold       = " << fixed << setprecision(2) << res[1] << "\tsize = " << fixed
                           << setprecision(2) << sizes[1];
                stdFileOut << "\tsquander ratio = " << sRatio[1] << endl;

                stdFileOut << "determined prob = " << fixed << setprecision(2) << res[2] << "\tsize = " << fixed
                           << setprecision(2) << sizes[2];
                stdFileOut << "\tsquander ratio = " << sRatio[2] << endl;

                cout << "apsize = " << apsize << " k = " << k << " rrsize = " << rrsize << " done!\n";
                stdFileOut << endl;
            }
        }
    }

    stdFileOut.close();
    return 0;
}