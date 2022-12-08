//
// Created by asuka on 2022/7/29.
//

#ifndef EXP_OPIM_OLD_H
#define EXP_OPIM_OLD_H

#include "graphs.h"

/*!
 * @brief Calculate the coverage of S in RRI, and return an upper/lower bound.
 * @param graph : the graph
 * @param k : the size of S
 * @param S : returns final seed set S as an passed parameter
 * @param RRI : the set of RR sets
 * @return return the upper bound of coverage S with R1.
 */
int64 OPIMGreedy_G(Graph &graph, int64 k, std::vector<int64> &S, RRContainer &RRI) {
    S.clear();
    std::vector<bool> RISetCovered(RRI.R.size(), false);
    std::vector<bool> nodeRemain(graph.n, true);
    auto *coveredNum_tmp = new int64[graph.n];
    memcpy(coveredNum_tmp, RRI.coveredNum, graph.n * sizeof(int64));
    std::vector<int64> coveredNum_maxK;
    int64 coverage = 0, A_u = INT64_MAX;

    for (int64 i = 0; i < k; i++) {
        coveredNum_maxK.clear();
        int64 maxInd = -1;
        for (int64 j = 0; j < graph.n; ++j) {
            if (nodeRemain[j]) {
                coveredNum_maxK.emplace_back(coveredNum_tmp[j]);
                if (maxInd == -1 || coveredNum_tmp[j] > coveredNum_tmp[maxInd]) maxInd = j;
            }
        }
        coverage += coveredNum_tmp[maxInd];
        S.emplace_back(maxInd);
        nodeRemain[maxInd] = false;
        for (int64 RIIndex : RRI.covered[maxInd]) {
            if (RISetCovered[RIIndex]) continue;
            for (int64 u : RRI.R[RIIndex]) {
                if (nodeRemain[u]) coveredNum_tmp[u]--;
            }
            RISetCovered[RIIndex] = true;
        }
        //find the max-k nodes with maximum marginal coverage
        int64 A_u_i = coverage;
        std::sort(coveredNum_maxK.begin(), coveredNum_maxK.end(), std::greater<>());
        for (int64 j = 0; j < std::min((int64) coveredNum_maxK.size(), k); j++) {
            A_u_i += coveredNum_maxK[j];
        }
        A_u = std::min(A_u, A_u_i);
    }
    delete[] coveredNum_tmp;
    printf("A_U: %ld %ld\n", A_u, coverage);
    return A_u;
}

/*!
 * @brief Use modified OPIM-C+ to find k nodes in candidate with most influence.
 * @param graph : the graph
 * @param k : the size of the seed set
 * @param eps : error threshold. The recommended value is 0.01.
 * @param delta : failure threshold. The recommended value is 1/n.
 * @param S : returns final S as a passed parameter
 */
void OPIM_G(Graph &graph, int64 k, int64 numOfRRSets, double delta, std::vector<int64> &S) {
    std::vector<int64> empty_set(0);
    RRContainer R1(graph, empty_set, true), R2(graph, empty_set, true);
    R1.resize(graph, numOfRRSets / 2);
    R2.resize(graph, numOfRRSets / 2);
    std::cout << "RR set done!\n";
    double d0 = log(2.0 / delta);
    auto upperC = (double) OPIMGreedy_G(graph, k, S, R1);
    auto lowerC = (double) R2.self_inf_cal(graph, S);
    double lower = sqr(sqrt(lowerC + 2.0 * d0 / 9.0) - sqrt(d0 / 2.0)) - d0 / 18.0;
    double upper = sqr(sqrt(upperC + d0 / 2.0) + sqrt(d0 / 2.0));
    double a0 = lower / upper;
    std::cout << "a0 = " << a0 << std::endl;
}

/*!
 * @brief Use modified OPIM-C+ to find k nodes in candidate with most influence.
 * @param graph : the graph
 * @param k : the size of the seed set
 * @param eps : error threshold. The recommended value is 0.01.
 * @param delta : failure threshold. The recommended value is 1/n.
 * @param S : returns final S as a passed parameter
 */
void OPIM_CG(Graph &graph, int64 k, double eps, double delta, std::vector<int64> &S) {
    std::vector<int64> empty_set(0);
    double approx = 1.0 - 1.0 / exp(1);
    double C_0 = 2.0 * sqr(
            approx * sqrt(log(6.0 / delta)) + sqrt(approx * (logcnk(graph.n, k) + log(6.0 / delta))));
    RRContainer R1(graph, empty_set, true), R2(graph, empty_set, true);
    R1.resize(graph, (size_t) C_0);
    R2.resize(graph, (size_t) C_0);
    auto i_max = (int64) (log2(1.0 * graph.n / eps / eps / k) + 1 - 1e-9);
    double d0 = log(3.0 * i_max / delta);
    for (int64 i = 1; i <= i_max; i++) {
        S.clear();
        auto upperC = (double) OPIMGreedy_G(graph, k, S, R1);
        auto lowerC = (double) R2.self_inf_cal(graph, S);
        double lower = sqr(sqrt(lowerC + 2.0 * d0 / 9.0) - sqrt(d0 / 2.0)) - d0 / 18.0;
        double upper = sqr(sqrt(upperC + d0 / 2.0) + sqrt(d0 / 2.0));
        double a0 = lower / upper;
        std::cout << "a0 = " << a0 << std::endl;
        if (a0 > approx - eps || i == i_max) break;
        R1.resize(graph, R1.R.size() * 2ll);
        R2.resize(graph, R2.R.size() * 2ll);
    }
    std::cout << "final RR sets = " << 2 * R1.R.size() << std::endl;
}


#endif //EXP_OPIM_OLD_H
