#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
//
// Created by lenovo on 2022/7/16.
//

#ifndef EXP_SIMULATION_H
#define EXP_SIMULATION_H

#include "graphs.h"


/*!
 * @brief generate a random node set of graph.
 * @param graph : the graph
 * @param A : stores the node set. Suppose it is initialized as empty.
 * @param size : the size of the node set
 */
void generate_ap(Graph &graph, std::vector<int64> &A, int64 size = 1) {
    A.clear();
    std::uniform_int_distribution<int64> uniformIntDistribution(0, graph.n - 1);
    for (int64 i = 0; i < size; i++) {
        int64 v = uniformIntDistribution(mt19937engine);
        while (std::find(A.begin(), A.end(), v) != A.end()) v = uniformIntDistribution(mt19937engine);
        A.emplace_back(v);
    }
}

/*!
 * @brief generate a random node set of graph.
 * @param graph : the graph
 * @param A : stores the node set. Suppose it is initialized as empty.
 * @param size : the size of the node set
 */
void generate_ap_by_degree(Graph &graph, std::vector<int64> &A, int64 size = 1) {
    A.clear();
    std::uniform_int_distribution<int64> uniformIntDistribution(0, graph.n - 1);
    for (int64 i = 0; i < size; i++) {
        int64 v = uniformIntDistribution(mt19937engine);
        while (graph.deg_out[v] <= graph.m / graph.n || std::find(A.begin(), A.end(), v) != A.end())
            v = uniformIntDistribution(mt19937engine);
        A.emplace_back(v);
    }
}

/*!
 * @brief run MC simulation to evaluate the influence spread.
 * @param graph : the graph that define propagation models(IC)
 * @param S : the seed set
 * @param Ap : the active participant
 * @return the estimated value of influence spread
 */
double MC_simulation(Graph &graph, std::vector<int64> &S, std::vector<int64> &Ap) {
    bool *active = new bool[graph.n];
    std::vector<int64> new_active, A, new_ones;
    std::vector<bool> Ap_bitwise(graph.n, false);
    for (auto u : Ap) Ap_bitwise[u] = true;
    double res = 0;
    for (int64 i = 1; i <= MC_iteration_rounds; i++) {
        if (graph.diff_model == IC) {
            new_active = S, A = S;
            for (int64 w : S) active[w] = true;
            new_ones.clear();
            while (!new_active.empty()) {
                for (int64 u : new_active) {
                    for (auto &edge : graph.g[u]) {
                        int64 v = edge.v;
                        if (Ap_bitwise[v]) continue;
                        if (active[v]) continue;
                        bool success = (random_real() < edge.p);
                        if (success) new_ones.emplace_back(v), active[v] = true;
                    }
                }
                new_active = new_ones;
                for (int64 u : new_ones) A.emplace_back(u);
                new_ones.clear();
            }
            for (int64 u : A) active[u] = false;
            res += (double) A.size() / MC_iteration_rounds;
            A.clear();
        }
    }
    delete[] active;
    return res;
}

/*!
 * @brief Calculate the degree of neighbor overlap at active participant.
 * @param graph : the graph
 * @param seeds : the active participant set
 * @param iteration_rounds : The number of selection
 * @return the mean overlap ratio
 */
double estimate_neighbor_overlap(Graph &graph, std::vector<int64> &seeds) {
    auto *num = new int64[graph.n]();
    int64 tot = 0, overlap = 0;
    for (int64 u : seeds)
        for (auto &e : graph.g[u])
            num[e.v]++;
    for (int64 u : seeds) {
        for (auto &e : graph.g[u]) {
            if (num[e.v] > 0) tot++;
            if (num[e.v] > 1) overlap++;
            num[e.v] = 0;
        }
    }
    delete[] num;
    return (double) overlap / tot;
}

/*!
 * @brief New calculation of influence spread in BIM, in which the effect of ap is excluded
 * @param graph : the graph that define propagation models(IC-M)
 * @param S : the seed set
 * @param A : the active participant
 * @return the estimated value of influence spread
 */
double FI_simulation_new(Graph &graph, std::vector<int64> &S, std::vector<int64> &A, int64 it_rounds = -1) {
    RRContainer RRI(graph, A, false);
    double res = 0, cur = clock();
    int64 it_ = (it_rounds == -1) ? MC_iteration_rounds : it_rounds;
    std::vector<int64> RR;
    for (int i = 1; i <= it_; i++) {
        RRI.RI_Gen(graph, S, RR);
        res += RR.size();
        RR.clear();
    }
    if (verbose_flag) printf("\t\tresult=%.3f time=%.3f\n", res, time_by(cur));
    return res / it_;
}

/// Efficiently estimate the influence spread with sampling error epsilon within probability 1-delta
double effic_inf(Graph &graph, std::vector<bi_node> &S, std::vector<int64> &A) {
    const double delta = 1e-3, eps = 0.01, c = 2.0 * (exp(1.0) - 2.0);
    const double LambdaL = 1.0 + 2.0 * c * (1.0 + eps) * log(2.0 / delta) / (eps * eps);
    size_t numHyperEdge = 0, numCoverd = 0;
    std::vector<bool> vecBoolSeed(graph.n), exclusive(graph.n);
    for (auto seed : S) vecBoolSeed[seed.first] = true;
    for (auto ap : A) exclusive[ap] = true;
    std::uniform_int_distribution<int64> uniformIntDistribution(0, graph.n - 1);
    std::vector<int64> nodes;
    std::queue<int64> Q;
    while (numCoverd < LambdaL) {
        numHyperEdge++;
        const auto uStart = uniformIntDistribution(mt19937engine);
        if (exclusive[uStart]) {
            continue;
        }
        if (vecBoolSeed[uStart]) {
            // Stop, this sample is covered
            numCoverd++;
            continue;
        }
        exclusive[uStart] = true;
        Q.push(uStart);
        nodes.emplace_back(uStart);
        bool break_flag = false;
        while (!Q.empty()) {
            int64 u = Q.front();
            Q.pop();
            for (auto &edgeT : graph.gT[u]) {
                if (exclusive[edgeT.v]) continue;
                if (random_real() < edgeT.p) {
                    if (vecBoolSeed[edgeT.v]) {
                        numCoverd++;
                        break_flag = true;
                        break;
                    }
                    exclusive[edgeT.v] = true;
                    Q.push(edgeT.v);
                    nodes.emplace_back(edgeT.v);
                }
            }
            if (break_flag) break;
        }
        for (auto e : nodes) exclusive[e] = false;
        nodes.clear();
        Q = std::queue<int64>(); // clear the queue
    }
    return 1.0 * numCoverd * graph.n / numHyperEdge;
}

#endif //EXP_SIMULATION_H

#pragma clang diagnostic pop