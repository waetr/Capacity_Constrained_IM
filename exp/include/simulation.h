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
    /// @brief Marks the point that was activated in the MC simulation
    bool *active = new bool[graph.n];
    double cur = clock();
    int64 meet_time = 0; //Too large for int32!
    std::vector<int64> new_active, A, new_ones;
    std::vector<Edge> meet_nodes, meet_nodes_tmp;
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
                        if (success) new_ones.emplace_back(v), active[v] = true, meet_time++;
                    }
                }
                new_active = new_ones;
                for (int64 u : new_ones) A.emplace_back(u);
                new_ones.clear();
            }
            for (int64 u : A) active[u] = false;
            res += (double) A.size() / MC_iteration_rounds;
            A.clear();
        } else if (graph.diff_model == IC_M) {
            meet_nodes.clear();
            new_active = S;
            for (int64 w : S) active[w] = true;
            for (int64 spread_rounds = 0; spread_rounds < graph.deadline; spread_rounds++) {
                for (int64 u : new_active) {
                    for (auto &edge : graph.g[u]) {
                        if (Ap_bitwise[edge.v]) continue;
                        if (active[edge.v]) continue;
                        meet_nodes.emplace_back(edge);
                    }
                }
                for (int64 u : new_active) A.emplace_back(u);
                new_active.clear();
                if (meet_nodes.empty()) break;
                meet_nodes_tmp.clear();
                for (auto &edge : meet_nodes) {
                    if (!active[edge.v]) {
                        meet_time++;
                        bool meet_success = (random_real() < edge.m);
                        if (meet_success) {
                            bool activate_success = (random_real() < edge.p);
                            if (activate_success) {
                                new_active.emplace_back(edge.v);
                                active[edge.v] = true;
                            }
                        } else {
                            meet_nodes_tmp.emplace_back(edge);
                        }
                    }
                }
                meet_nodes = meet_nodes_tmp;
            }
            for (int64 u : new_active) A.emplace_back(u);
            for (int64 u : A) active[u] = false;
            res += (double) A.size() / MC_iteration_rounds;
            A.clear();
        }
    }
    if (verbose_flag) {
        std::cout << "\t\tresult=" << res << " time=" << time_by(cur) << " meet time=" << meet_time << std::endl;
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
    std::vector<int64> RR;
    double res = 0, cur = clock();
    int64 it_ = (it_rounds == -1) ? MC_iteration_rounds : it_rounds;
    for (int i = 1; i <= it_; i++) {
        RRI.RI_Gen(graph, S, RR);
        res += RR.size();
    }
    if (verbose_flag) printf("\t\tresult=%.3f time=%.3f\n", res, time_by(cur));
    return res / it_;
}

double FI_simulation_binode(Graph &graph, std::vector<bi_node> &S, std::vector<int64> &A, int64 it_rounds = -1) {
    RRContainer RRI(graph, A, false);
    std::vector<int64> RR, S_(S.size());
    for (int i = 0; i < S.size(); i++) S_[i] = S[i].first;
    double res = 0, cur = clock();
    int64 it_ = (it_rounds == -1) ? MC_iteration_rounds : it_rounds;
    for (int i = 1; i <= it_; i++) {
        RRI.RI_Gen(graph, S_, RR);
        res += RR.size();
    }
    return res / it_;
}

#endif //EXP_SIMULATION_H
