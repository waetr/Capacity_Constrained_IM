//
// Created by asuka on 2022/7/24.
//

#ifndef SIMULATION_H_IMM_H
#define SIMULATION_H_IMM_H

#include "graphs.h"

/*!
 * @brief temporary array for IMM_NodeSelection
 *
 * RISetCovered : marking the RI sets that has been covered by S
 * nodeRemain : check if an node is remained in the heap
 * covered_num_tmp : a copy of covered_num
 */

/*!
 * @brief Selection phase of IMM : Select a set S of size k that covers the maximum RI sets in R
 * @param graph : the graph
 * @param candidate : candidate node set in which seed set can choose nodes
 * @param k : the size of S
 * @param S : returns S as an passed parameter
 * @param RRI : the set of RR sets
 * @return : the fraction of RI sets in R that are covered by S
 */
double IMMNodeSelection(Graph &graph, std::vector<int64> &candidate, int64 k, std::vector<int64> &S, RRContainer &RRI) {
    S.clear();
    std::vector<bool> RISetCovered(RRI.R.size(), false);
    std::vector<bool> nodeRemain(graph.n, false);
    for (int64 i : candidate) nodeRemain[i] = true;
    auto *coveredNum_tmp = new int64[graph.n];
    memcpy(coveredNum_tmp, RRI.coveredNum, graph.n * sizeof(int64));
    std::priority_queue<std::pair<int64, int64>> Q;
    for (int64 i : candidate) Q.push(std::make_pair(coveredNum_tmp[i], i));
    int64 influence = 0;

    while (S.size() < k && !Q.empty()) {
        int64 value = Q.top().first;
        int64 maxInd = Q.top().second;
        Q.pop();
        if (value > coveredNum_tmp[maxInd]) {
            Q.push(std::make_pair(coveredNum_tmp[maxInd], maxInd));
            continue;
        }
        influence += coveredNum_tmp[maxInd];
        S.emplace_back(maxInd);
        nodeRemain[maxInd] = false;
        for (int64 RIIndex : RRI.covered[maxInd]) {
            if (RISetCovered[RIIndex]) continue;
            for (int64 u : RRI.R[RIIndex]) {
                if (nodeRemain[u]) coveredNum_tmp[u]--;
            }
            RISetCovered[RIIndex] = true;
        }
    }
    delete[] coveredNum_tmp;
    return (double) influence / RRI.R.size();
}

/*!
 * @brief Sampling phase of IMM : generate sufficient RI sets into R.
 * @param graph : the graph
 * @param candidate : candidate node set in which seed set can choose nodes
 * @param k : the size of the seed set
 * @param eps : argument related to accuracy.
 * @param iota : argument related to accuracy.
 * @param RRI : the set of RR sets
 */
void IMMSampling(Graph &graph, std::vector<int64> &candidate, int64 k, double eps, double iota, RRContainer &RRI) {
    double epsilon_prime = eps * sqrt(2);
    double LB = 1;
    std::vector<int64> S_tmp;
    auto End = (int) (log2(graph.n) + 1e-9 - 1);
    for (int i = 1; i <= End; i++) {
        auto ci = (int64) ((2.0 + 2.0 / 3.0 * epsilon_prime) *
                           (iota * log(graph.n) + logcnk(graph.n, k) + log(log2(graph.n))) /
                           sqr(epsilon_prime) *
                           pow(2.0, i));
        if (ci > (int64) 200000000) break;
        RRI.resize(graph, ci);

        double ept = IMMNodeSelection(graph, candidate, k, S_tmp, RRI);
        if (ept > 1.0 / pow(2.0, i)) {
            LB = ept * graph.n / (1.0 + epsilon_prime);
            break;
        }
    }
    double e = exp(1);
    double alpha = sqrt(iota * log(graph.n) + log(2));
    double beta = sqrt((1.0 - 1.0 / e) * (logcnk(graph.n, k) + iota * log(graph.n) + log(2)));
    auto C = (int64) (2.0 * graph.n * sqr((1.0 - 1.0 / e) * alpha + beta) / LB / sqr(eps));
    RRI.resize(graph, std::min(C, (int64) 200000000));
    if (verbose_flag) {
        std::cout << "\tfinal C = " << C << std::endl;
    }
}

/*!
 * @brief Use IMM to find k nodes in candidate with most influence.
 * @param graph : the graph
 * @param candidate : candidate node set in which seed set can choose nodes
 * @param k : the size of the seed set
 * @param eps : approximation argument. default as 0.5.
 * @param iota : argument related to failure probability. default as 1.
 * @param S : returns final S as a passed parameter
 * @param RRI : the set of RR sets
 */
void IMM(Graph &G, std::vector<int64> &candidate, int64 k, double eps, double iota, std::vector<int64> &S, RRContainer &RRI) {
    double iota_new = iota * (1.0 + log(2) / log(G.n));
    IMMSampling(G, candidate, k, eps, iota_new, RRI);
    IMMNodeSelection(G, candidate, k, S, RRI);
}

/*!
 * @brief Encapsulated operations for Option 2 using IM solver : IMM
 * @param graph : the graph
 * @param k : the number in the problem definition
 * @param A : the active participant set A
 * @param seeds : returns the seed set S = {S_1, S_2, ..., S_n}
 */
void IMM_method(Graph &graph, int64 k, std::vector<int64> &A, std::vector<int64> &seeds) {
    double cur = clock();
    std::set<int64> seeds_reorder;
    for (int64 u : A) {
        std::vector<int64> neighbours, one_seed;
        for (auto &edge : graph.g[u]) {
            if (find(A.begin(), A.end(), edge.v) == A.end()) {
                neighbours.emplace_back(edge.v);
            }
        }
        RRContainer RRI(graph, A, true);
        IMM(graph, neighbours, k, 0.5, 1, one_seed, RRI);
        for (int64 w : one_seed)
            seeds_reorder.insert(w);
    }
    for (int64 w : seeds_reorder) seeds.emplace_back(w);
    seeds_reorder.clear();
    if (verbose_flag) printf("IMM method done. total time = %.3f\n", time_by(cur));
}


/*!
 * @brief Selection phase of advanced IMM.
 * @param graph : the graph
 * @param A : active participant set
 * @param k : the size of S
 * @param S : returns S as an passed parameter
 * @param RRI : the set of RR sets
 * @return : the fraction of RI sets in R that are covered by S
 */
double IMMNodeSelection_advanced(Graph &graph, std::vector<int64> &A, int64 k, std::vector<int64> &S, RRContainer &RRI) {
    S.clear();
    CandidateNeigh candidate(graph, A, k);
    std::vector<bool> RISetCovered(RRI.R.size(), false);
    std::vector<bool> nodeRemain(graph.n, false);
    for (int64 i : candidate.N) nodeRemain[i] = true;
    auto *coveredNum_tmp = new int64[graph.n];
    memcpy(coveredNum_tmp, RRI.coveredNum, graph.n * sizeof(int64));
    std::priority_queue<std::pair<int64, int64>> Q;
    for (int64 i : candidate.N) Q.push(std::make_pair(coveredNum_tmp[i], i));
    int64 influence = 0;

    while (!Q.empty()) {
        int64 value = Q.top().first;
        int64 maxInd = Q.top().second;
        Q.pop();

        int64 u0 = candidate.source_participant(maxInd);
        if (u0 == -1) continue;

        if (value > coveredNum_tmp[maxInd]) {
            Q.push(std::make_pair(coveredNum_tmp[maxInd], maxInd));
            continue;
        }
        influence += coveredNum_tmp[maxInd];
        S.emplace_back(maxInd);
        candidate.choose(u0);
        nodeRemain[maxInd] = false;
        for (int64 RIIndex : RRI.covered[maxInd]) {
            if (RISetCovered[RIIndex]) continue;
            for (int64 u : RRI.R[RIIndex]) {
                if (nodeRemain[u]) coveredNum_tmp[u]--;
            }
            RISetCovered[RIIndex] = true;
        }
    }
    delete[] coveredNum_tmp;
    return (double) influence / RRI.R.size();
}

/*!
 * @brief Sampling phase of IMM : generate sufficient RI sets into R.
 * @param graph : the graph
 * @param candidate : candidate node set in which seed set can choose nodes
 * @param k : the size of the seed set
 * @param eps : argument related to accuracy.
 * @param iota : argument related to accuracy.
 * @param RRI : the set of RR sets
 */
void IMMSampling_advanced(Graph &graph, std::vector<int64> &A, int64 k, double eps, double iota, RRContainer &RRI) {
    int64 kA = 0;
    for (int64 u : A) kA += std::min(k, (int64)graph.g[u].size());
    kA = std::min(kA, (int64) graph.n);
    double epsilon_prime = eps * sqrt(2);
    double LB = 1;
    std::vector<int64> S_tmp;
    auto End = (int) (log2(graph.n) + 1e-9 - 1);
    for (int i = 1; i <= End; i++) {
        auto ci = (int64) ((2.0 + 2.0 / 3.0 * epsilon_prime) *
                           (iota * log(graph.n) + logcnk(graph.n, kA) + log(log2(graph.n))) *
                           pow(2.0, i) / sqr(epsilon_prime));
        if (verbose_flag) std::cout << "\tci = " << ci << std::endl;
        RRI.resize(graph, ci);

        double ept = IMMNodeSelection_advanced(graph, A, k, S_tmp, RRI);
        if (ept > 1.0 / pow(2.0, i)) {
            LB = ept * graph.n / (1.0 + epsilon_prime);
            break;
        }
    }
    double e = exp(1);
    double alpha = sqrt(iota * log(graph.n) + log(2));
    double beta = sqrt((1.0 - 1.0 / e) * (logcnk(graph.n, kA) + iota * log(graph.n) + log(2)));
    auto C = (int64) (2.0 * graph.n * sqr((1.0 - 1.0 / e) * alpha + beta) / LB / sqr(eps));
    RRI.resize(graph, C);
    if (verbose_flag) std::cout << "\tfinal C = " << C << std::endl;
}

/*!
 * @brief the same as IMM.
 */
void IMM_advanced(Graph &G, std::vector<int64> &A, int64 k, double eps, double iota, std::vector<int64> &S, RRContainer &RRI) {
    double iota_new = iota * (1.0 + log(2) / log(G.n));
    IMMSampling_advanced(G, A, k, eps, iota_new, RRI);
    IMMNodeSelection_advanced(G, A, k, S, RRI);
}

/*!
 * @brief Encapsulated operations for Greedy-IMM
 * @param graph : the graph
 * @param k : the number in the problem definition
 * @param A : the active participant set A
 * @param seeds : returns the seed set S = {S_1, S_2, ..., S_n}
 */
void advanced_IMM_method(Graph &graph, int64 k, std::vector<int64> &A, std::vector<int64> &seeds) {
    RRContainer RRI(graph, A, true);
    double cur = clock();
    IMM_advanced(graph, A, k, 0.5, 1, seeds, RRI);
    if (verbose_flag) printf("IMM advanced done. total time = %.3f\n", time_by(cur));
}

#endif //SIMULATION_H_IMM_H
