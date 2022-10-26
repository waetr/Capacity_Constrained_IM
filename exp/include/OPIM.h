//
// Created by asuka on 2022/7/29.
//

#ifndef EXP_OPIM_H
#define EXP_OPIM_H

#include "graphs.h"

/*!
 * @brief Calculate the coverage of S in RRI, and return an upper/lower bound.
 * @param graph : the graph
 * @param candidate : candidate node set in which seed set can choose nodes
 * @param k : the size of S
 * @param S : returns final seed set S as an passed parameter
 * @param RRI : the set of RR sets
 * @return return the upper bound of coverage S with R1.
 */
int64 OPIMGreedy(Graph &graph, std::vector<int64> &candidate, int64 k, std::vector<int64> &S, RRContainer &RRI) {
    S.clear();
    std::vector<bool> RISetCovered(RRI.R.size(), false);
    for (int64 i : candidate) nodeRemain[i] = true;
    auto *coveredNum_tmp = new int64[graph.n];
    memcpy(coveredNum_tmp, RRI.coveredNum, graph.n * sizeof(int64));
    std::vector<int64> coveredNum_maxK;
    int64 coverage = 0, A_u = INT64_MAX;

    for (int64 i = 0; i <= k; i++) {
        if (i > 0) {
            int64 maxInd = -1;
            for (int64 v : candidate) {
                if (!nodeRemain[v]) continue;
                if (maxInd == -1 || coveredNum_tmp[v] > coveredNum_tmp[maxInd]) maxInd = v;
            }
            if (maxInd == -1) break;

            coverage += coveredNum_tmp[maxInd];
            S.emplace_back(maxInd);
            nodeRemain[maxInd] = false;
            for (auto RIIndex : RRI.covered[maxInd]) {
                if (RISetCovered[RIIndex]) continue;
                for (auto u : RRI.R[RIIndex]) {
                    if (nodeRemain[u]) coveredNum_tmp[u]--;
                }
                RISetCovered[RIIndex] = true;
            }
        }
        //find the max-k nodes with maximum marginal coverage
        int64 A_u_i = coverage;
        coveredNum_maxK.clear();
        for (int64 u : candidate) {
            if (nodeRemain[u])
                coveredNum_maxK.emplace_back(coveredNum_tmp[u]);
        }
        if (coveredNum_maxK.size() <= k) {
            A_u_i += accumulate(coveredNum_maxK.begin(), coveredNum_maxK.end(), (int64) 0);
        } else {
            nth_element(coveredNum_maxK.begin(), coveredNum_maxK.end() - k, coveredNum_maxK.end());
            for (int64 j = coveredNum_maxK.size() - k; j < coveredNum_maxK.size(); j++) {
                A_u_i += coveredNum_maxK[j];
            }
        }
        A_u = std::min(A_u, A_u_i);
    }
    delete[] coveredNum_tmp;
    for (int64 i : candidate) nodeRemain[i] = false;
    return A_u;
}

/*!
 * @brief Use OPIM-C+ to find k nodes in candidate with most influence.
 * @param graph : the graph
 * @param candidate : candidate node set in which seed set can choose nodes
 * @param k : the size of the seed set
 * @param eps : error threshold. The recommended value is 0.01.
 * @param delta : failure threshold. The recommended value is 1/n.
 * @param S : returns final S as a passed parameter
 */
void OPIM_C(Graph &graph, std::vector<int64> &A, std::vector<int64> &candidate, int64 k, double eps, double delta,
            std::vector<int64> &S) {
    std::vector<int64> S_tmp;
    double approx = 1.0 - 1.0 / exp(1);
    double C_0 = 2.0 * sqr(
            approx * sqrt(log(6.0 / delta)) + sqrt(approx * (logcnk(graph.n, k) + log(6.0 / delta))));
    RRContainer R1(graph, A, true), R2(graph, A, true);
    R1.resize(graph, (size_t) C_0);
    R2.resize(graph, (size_t) C_0);
    auto i_max = (int64) (log2(1.0 * graph.n / eps / eps / k) + 1 - 1e-9);
    double d0 = log(3.0 * i_max / delta);
    for (int64 i = 1; i <= i_max; i++) {
        auto upperC = (double) OPIMGreedy(graph, candidate, k, S, R1);
        auto lowerC = (double) R2.self_inf_cal(graph, S);
        double lower = sqr(sqrt(lowerC + 2.0 * d0 / 9.0) - sqrt(d0 / 2.0)) - d0 / 18.0;
        double upper = sqr(sqrt(upperC + d0 / 2.0) + sqrt(d0 / 2.0));
        double a0 = lower / upper;
        if (verbose_flag) std::cout << "a0 = " << a0 << std::endl;
        if (a0 > approx - eps || i == i_max) break;
        R1.resize(graph, R1.R.size() * 2ll);
        R2.resize(graph, R2.R.size() * 2ll);
    }
    if (verbose_flag) std::cout << "final RR sets = " << 2 * R1.R.size() << std::endl;
}

/*!
 * @brief Encapsulated operations for Option 2 using IM solver : OPIM
 * @param graph : the graph
 * @param k : the number in the problem definition
 * @param A : the active participant set A
 * @param seeds : returns the seed set S = {S_1, S_2, ..., S_n}
 */
void OPIM_method(Graph &graph, int64 k, std::vector<int64> &A, std::vector<int64> &seeds) {
    double cur = clock();
    std::set<int64> seeds_reorder;
    for (int64 u : A) {
        std::vector<int64> neighbours, one_seed;
        for (auto &edge : graph.g[u]) {
            if (find(A.begin(), A.end(), edge.v) == A.end()) {
                neighbours.emplace_back(edge.v);
            }
        }
        OPIM_C(graph, A, neighbours, k, 0.01, 1.0 / graph.n, one_seed);
        for (int64 w : one_seed)
            seeds_reorder.insert(w);
    }
    for (int64 w : seeds_reorder) seeds.emplace_back(w);
    seeds_reorder.clear();
    if (verbose_flag) printf("OPIM method done. total time = %.3f\n", time_by(cur));
}

/*!
 * @brief Calculate the coverage of S in RRI, and return an upper/lower bound.
 * @param graph : the graph
 * @param candidate : candidate node set in which seed set can choose nodes
 * @param k : the size of S
 * @param S : returns final seed set S as an passed parameter
 * @param RRI : the set of RR sets
 * @return return the upper bound of coverage S with R1.
 */
int64 OPIMGreedy_G(Graph &graph, std::vector<int64> &A, int64 k, std::vector<int64> &S, RRContainer &RRI) {
    int64 kA = 0;
    for (int64 u : A) kA += std::min(k, (int64) graph.g[u].size());
    kA = std::min(kA, (int64) graph.n);
    S.clear();
    CandidateNeigh candidate(graph, A, k);
    std::vector<bool> RISetCovered(RRI.R.size(), false);
    for (int64 i : candidate.N) nodeRemain[i] = true;
    auto *coveredNum_tmp = new int64[graph.n];
    memcpy(coveredNum_tmp, RRI.coveredNum, graph.n * sizeof(int64));
    std::vector<int64> coveredNum_maxK;
    int64 coverage = 0, A_u = INT64_MAX;

    for (int64 i = 0;; i++) {
        if (i > 0) {
            int64 maxInd = -1;
            for (int64 v : candidate.N) {
                if (!nodeRemain[v]) continue;
                int64 u = candidate.source_participant(v);
                if (u == -1) {
                    nodeRemain[v] = false;
                    continue;
                }
                if (maxInd == -1 || coveredNum_tmp[v] > coveredNum_tmp[maxInd]) maxInd = v;
            }
            if (maxInd == -1) break;
            int64 u0 = candidate.source_participant(maxInd);

            coverage += coveredNum_tmp[maxInd];
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
        //find the max-k nodes with maximum marginal coverage
        int64 A_u_i = coverage;
        coveredNum_maxK.clear();
        for (int64 u : candidate.N) {
            if (nodeRemain[u])
                coveredNum_maxK.emplace_back(coveredNum_tmp[u]);
        }
        if (coveredNum_maxK.size() >= kA) {
            nth_element(coveredNum_maxK.begin(), coveredNum_maxK.end() - kA, coveredNum_maxK.end());
        }
        size_t begin0 = coveredNum_maxK.size() <= kA ? 0 : coveredNum_maxK.size() - kA;
        for (int64 j = begin0; j < coveredNum_maxK.size(); j++) {
            A_u_i += coveredNum_maxK[j];
        }
        A_u = std::min(A_u, A_u_i);
    }
    delete[] coveredNum_tmp;
    for (int64 i : candidate.N) nodeRemain[i] = false;
    return A_u;
}

/*!
 * @brief Use modified OPIM-C+ to find k nodes in candidate with most influence.
 * @param graph : the graph
 * @param A : active participant
 * @param k : the size of the seed set
 * @param eps : error threshold. The recommended value is 0.01.
 * @param delta : failure threshold. The recommended value is 1/n.
 * @param S : returns final S as a passed parameter
 */
void OPIM_CG(Graph &graph, std::vector<int64> &A, int64 k, double eps, double delta, std::vector<int64> &S) {
    int64 kA = 0;
    for (int64 u : A) kA += std::min(k, (int64) graph.g[u].size());
    kA = std::min(kA, (int64) graph.n);
    std::vector<int64> S_tmp;
    double approx = 1.0 - 1.0 / exp(1);
    double C_0 = 2.0 * sqr(
            approx * sqrt(log(6.0 / delta)) + sqrt(approx * (logcnk(graph.n, kA) + log(6.0 / delta))));
    RRContainer R1(graph, A, true), R2(graph, A, true);
    R1.resize(graph, (size_t) C_0);
    R2.resize(graph, (size_t) C_0);
    auto i_max = (int64) (log2(1.0 * graph.n / eps / eps / kA) + 1 - 1e-9);
    double d0 = log(3.0 * i_max / delta);
    for (int64 i = 1; i <= i_max; i++) {
        auto upperC = (double) OPIMGreedy_G(graph, A, k, S, R1);
        auto lowerC = (double) R2.self_inf_cal(graph, S);
        double lower = sqr(sqrt(lowerC + 2.0 * d0 / 9.0) - sqrt(d0 / 2.0)) - d0 / 18.0;
        double upper = sqr(sqrt(upperC + d0 / 2.0) + sqrt(d0 / 2.0));
        double a0 = lower / upper;
        if (verbose_flag) std::cout << "a0 = " << a0 << std::endl;
        if (a0 > approx - eps || i == i_max) break;
        R1.resize(graph, R1.R.size() * 2ll);
        R2.resize(graph, R2.R.size() * 2ll);
    }
    if (verbose_flag) std::cout << "final RR sets = " << 2 * R1.R.size() << std::endl;
}

/*!
 * @brief Encapsulated operations for Greedy-OPIM
 * @param graph : the graph
 * @param k : the number in the problem definition
 * @param A : the active participant set A
 * @param seeds : returns the seed set S = {S_1, S_2, ..., S_n}
 */
void advanced_OPIM_method(Graph &graph, int64 k, std::vector<int64> &A, std::vector<int64> &seeds) {
    double cur = clock();
    OPIM_CG(graph, A, k, 0.01, 1.0 / graph.n, seeds);
    if (verbose_flag) printf("OPIM advanced done. total time = %.3f\n", time_by(cur));
}

#endif //EXP_OPIM_H
