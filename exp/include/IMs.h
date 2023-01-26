// IM Solvers

#ifndef EXP_IMS_H
#define EXP_IMS_H

#include "simulation.h"
#include <map>


double method_random(Graph &graph, int64 k, std::vector<int64> &A, std::vector<bi_node> &seeds) {
    assert(seeds.empty());
    CandidateNeigh candidate(graph, A, k);
    std::vector<int64> set0 = candidate.N;
    auto start_time = std::chrono::high_resolution_clock::now();
    //S_ordered is ordered by pageRank
    std::shuffle(set0.begin(), set0.end(), std::mt19937(std::random_device()()));
    for (auto v : set0) {
        int64 u = candidate.source_participant(v);
        if (u == -1) continue;
        candidate.choose(u);
        seeds.emplace_back(v, u);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}


/*!
 * @brief Use power iteration method to calculate pagerank values of nodes in graph.
 * @remarks Note that if there are isolated nodes in the graph, the iteration will not stop and will lead to an error.
 * @param graph : the graph
 * @param pi : returns a size-n vector, containing pagerank values of all nodes
 * @param alpha : initialized usually as 0.15 or 0.2
 * @param l1_error : The precision that needs to be achieved
 */
void power_iteration(Graph &graph, std::vector<double> &pi, double alpha, double l1_error = 1e-9) {
    std::vector<double> residuals(graph.n, 1.0 / graph.n);
    std::vector<double> new_residuals(graph.n, 0);
    double r_sum = 1;
    while (r_sum > l1_error) {
        r_sum = 0;
        for (int64 u = 0; u < graph.n; u++) {
            for (auto &edge : graph.g[u]) {
                int v = edge.v;
                new_residuals[v] += residuals[u] / graph.deg_out[u];
            }
        }
        for (int64 u = 0; u < graph.n; u++) {
            new_residuals[u] = (1.0 - alpha) * new_residuals[u] + (alpha * 1.0 / graph.n);
            r_sum += fabs(residuals[u] - new_residuals[u]);
            residuals[u] = new_residuals[u];
            new_residuals[u] = 0;
        }
    }
    pi = residuals;
}

/*!
 * @brief Encapsulated method using local-PageRank
 * @param graph : the graph
 * @param k : the number in the problem definition
 * @param A : the active participant set A
 * @param seeds : returns the seed set S (each element is a pair <node, AP>)
 */
double method_local_PageRank(Graph &graph, int64 k, std::vector<int64> &A, std::vector<bi_node> &seeds) {
    assert(seeds.empty());
    std::vector<double> pi(graph.n, 0);
    std::vector<std::pair<double, int64> > pg_rank;
    std::set<int64> seeds_reorder, A_ordered(A.begin(), A.end());

    auto start_time = std::chrono::high_resolution_clock::now();
    power_iteration(graph, pi, 0.2);
    for (int64 u : A) {
        pg_rank.clear();
        for (auto &edge : graph.g[u]) {
            if (A_ordered.find(edge.v) == A_ordered.end()) {
                pg_rank.emplace_back(std::make_pair(pi[edge.v], edge.v));
            }
        }
        sort(pg_rank.begin(), pg_rank.end(), std::greater<>());
        for (int64 i = 0; i < k && i < pg_rank.size(); i++) {
            int64 w = pg_rank[i].second;
            if (seeds_reorder.find(w) == seeds_reorder.end()) {
                seeds_reorder.insert(w);
                seeds.emplace_back(w, u);
            }
        }
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}

/*!
 * @brief Encapsulated method using local-degree
 * @param graph : the graph
 * @param k : the number in the problem definition
 * @param A : the active participant set A
 * @param seeds : returns the seed set S (each element is a pair <node, AP>)
 */
double method_local_Degree(Graph &graph, int64 k, std::vector<int64> &A, std::vector<bi_node> &seeds) {
    assert(seeds.empty());
    std::vector<std::pair<double, int64> > dg_rank;
    std::set<int64> seeds_reorder, A_ordered(A.begin(), A.end());

    auto start_time = std::chrono::high_resolution_clock::now();
    for (int64 u : A) {
        dg_rank.clear();
        for (auto &edge : graph.g[u]) {
            if (A_ordered.find(edge.v) == A_ordered.end()) {
                dg_rank.emplace_back(std::make_pair(graph.deg_out[edge.v], edge.v));
            }
        }
        sort(dg_rank.begin(), dg_rank.end(), std::greater<>());
        for (int64 i = 0; i < k && i < dg_rank.size(); i++) {
            int64 w = dg_rank[i].second;
            if (seeds_reorder.find(w) == seeds_reorder.end()) {
                seeds_reorder.insert(w);
                seeds.emplace_back(w, u);
            }
        }
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}

/*!
 * @brief Encapsulated method using greedy-CELF
 * @param graph : the graph
 * @param k : the number in the problem definition
 * @param A : the active participant set A
 * @param seeds : returns the seed set S (each element is a pair <node, AP>)
 */
double method_greedy_CELF(Graph &graph, int64 k, std::vector<int64> &A, std::vector<bi_node> &seeds) {
    assert(seeds.empty());
    /// initialization
    CandidateNeigh candidate(graph, A, k);
    typedef std::pair<double, std::pair<int64, int64> > node0;
    std::priority_queue<node0> Q;
    std::vector<int64> seeds_calc(1); //Add a temporary space
    double current_spread = 0;

    /// main
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int64 u : candidate.N) {
        seeds_calc[0] = u;
        Q.push(make_pair(MC_simulation(graph, seeds_calc, A), std::make_pair(u, 0)));
    }
    seeds_calc.clear(); //Clear the temporary space
    while (!Q.empty()) {
        node0 Tp = Q.top();
        int64 v = Tp.second.first;
        int64 it_round = Tp.second.second;
        double mg = Tp.first;
        Q.pop();
        int64 u = candidate.source_participant(v);
        if (u == -1) continue;
        if (it_round == seeds.size()) {
            candidate.choose(u);
            seeds.emplace_back(v, u);
            seeds_calc.emplace_back(v);
            current_spread += mg;
        } else {
            seeds_calc.emplace_back(v);
            Tp.first = MC_simulation(graph, seeds_calc, A) - current_spread;
            seeds_calc.pop_back();
            Tp.second.second = seeds.size();
            Q.push(Tp);
        }
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}

/*!
 * @brief Encapsulated method using greedy-vanilla
 * @param graph : the graph
 * @param k : the number in the problem definition
 * @param A : the active participant set A
 * @param seeds : returns the seed set S (each element is a pair <node, AP>)
 */
double method_greedy_vanilla(Graph &graph, int64 k, std::vector<int64> &A, std::vector<bi_node> &seeds) {
    assert(seeds.empty());
    /// initialization
    std::vector<bi_node> candidate;
    std::vector<int64> num(A.size(), 0);
    for (int i = 0; i < A.size(); i++) {
        for (auto &e: graph.g[A[i]]) {
            if (std::find(A.begin(), A.end(), e.v) == A.end()) candidate.emplace_back(e.v, i);
        }
    }
    std::vector<bool> node_selected(graph.n, false);
    std::vector<int64> seeds_calc(0); //Add a temporary space
    double current_spread = 0;
    /// main
    auto start_time = std::chrono::high_resolution_clock::now();
    while (true) {
        int64 u, v;
        double mg = -1;
        for (auto pair: candidate) {
            auto w = pair.first, ap_index = pair.second;
            if(node_selected[w] || num[ap_index] >= k) continue;
            seeds_calc.emplace_back(w);
            double mg0 = MC_simulation(graph, seeds_calc, A) - current_spread;
            if (mg0 > mg) v = w, u = ap_index, mg = mg0;
            seeds_calc.pop_back();
        }
        if (mg < 0) break;
        node_selected[v] = true;
        num[u] += 1;
        seeds.emplace_back(v, A[u]);
        seeds_calc.emplace_back(v);
        current_spread += mg;
    }
    //printf("thr: %d\n", r1);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}

/*!
 * @brief Encapsulated method using RR_vanilla
 * @param graph : the graph
 * @param k : the number in the problem definition
 * @param A : the active participant set A
 * @param seeds : returns the seed set S (each element is a pair <node, AP>)
 */
double method_RR_vanilla(Graph &graph, int64 k, std::vector<int64> &A, std::vector<bi_node> &seeds) {
    assert(seeds.empty());
    ///initialization
    int64 current_influence = 0, N_empty = 0;
    std::vector<int64> Ni_empty(A.size(), 0);
    std::vector<bool> selected(graph.n, false);
    std::vector<int64> seeds_calc(0); //Add a temporary space

    auto start_time = std::chrono::high_resolution_clock::now();
    for (auto ap : A) selected[ap] = true;
    while (N_empty < A.size()) {
        for (int i = 0; i < A.size(); i++) { ///N_numbers[i] == k + 1 means that N[i] is full
            if (Ni_empty[i] != k + 1 && Ni_empty[i] == k) {
                Ni_empty[i] = k + 1;
                N_empty++;
            }
            if (Ni_empty[i] == k + 1) continue;
            int64 v = -1;
            double v_inf = 0;
            for(auto &e : graph.g[A[i]]) {
                if (selected[e.v]) continue;
                seeds_calc.emplace_back(e.v);
                double v_inf_tmp = MC_simulation(graph, seeds_calc, A) - current_influence;
                if (v_inf_tmp > v_inf) v = e.v, v_inf = v_inf_tmp;
                seeds_calc.pop_back();
            }
            if (Ni_empty[i] != k + 1 && v == -1) {
                Ni_empty[i] = k + 1;
                N_empty++;
            }
            if (Ni_empty[i] == k + 1) continue;
            ///choose v
            seeds.emplace_back(v, A[i]);
            seeds_calc.emplace_back(v);
            current_influence += v_inf;
            selected[v] = true;
            Ni_empty[i]++;
        }
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}

/*!
 * @brief Encapsulated method using RR_CELF
 * @param graph : the graph
 * @param k : the number in the problem definition
 * @param A : the active participant set A
 * @param seeds : returns the seed set S (each element is a pair <node, AP>)
 */
double method_RR_CELF(Graph &graph, int64 k, std::vector<int64> &A, std::vector<bi_node> &seeds) {
    assert(seeds.empty());
    ///initialization
    typedef std::pair<double, std::pair<int64, int64> > node0;
    std::priority_queue<node0> Q[A.size()];
    std::map<int64, double> marginal_influence;
    std::set<int64> A_reorder(A.begin(), A.end());
    std::vector<int64> seeds_calc(1); //Add a temporary space
    int64 current_influence = 0, N_empty = 0;
    std::vector<int64> Ni_empty(A.size(), 0);
    std::vector<bool> selected(graph.n, false);

    /// main
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < A.size(); i++) {
        for (auto e : graph.g[A[i]])
            if (A_reorder.find(e.v) == A_reorder.end()) {
                auto tmp_mg = marginal_influence[e.v];
                if (tmp_mg != 0) {
                    Q[i].push(std::make_pair(tmp_mg, std::make_pair(e.v, 0)));
                } else {
                    seeds_calc[0] = e.v;
                    tmp_mg = MC_simulation(graph, seeds_calc, A);
                    marginal_influence[e.v] = tmp_mg;
                    Q[i].push(std::make_pair(tmp_mg, std::make_pair(e.v, 0)));
                }
            }
    }
    seeds_calc.clear(); //Clear the temporary space
    while (N_empty < A.size()) {
        for (int i = 0; i < A.size(); i++) { ///N_numbers[i] == k + 1 means that N[i] is full
            if (Ni_empty[i] != k + 1 && Ni_empty[i] == k) {
                Ni_empty[i] = k + 1;
                N_empty++;
            }
            if (Ni_empty[i] == k + 1) continue;
            /// skip all nodes that cannot be selected, since we must select one in this for-loop
            while (!Q[i].empty() && (selected[Q[i].top().second.first] || Q[i].top().second.second != seeds.size())) {
                node0 Tp = Q[i].top();
                Q[i].pop();
                if (selected[Tp.second.first]) continue;
                if (Tp.second.second != seeds.size()) {
                    seeds_calc.emplace_back(Tp.second.first);
                    Tp.first = MC_simulation(graph, seeds_calc, A) - current_influence;
                    seeds_calc.pop_back();
                    Tp.second.second = seeds.size();
                    Q[i].push(Tp);
                }
            }
            if (Ni_empty[i] != k + 1 && Q[i].empty()) {
                Ni_empty[i] = k + 1;
                N_empty++;
            }
            if (Ni_empty[i] == k + 1) continue;
            int64 v = Q[i].top().second.first;
            double mg = Q[i].top().first;
            Q[i].pop();
            ///choose v
            seeds.emplace_back(v, A[i]);
            seeds_calc.emplace_back(v);
            current_influence += mg;
            selected[v] = true;
            Ni_empty[i]++;
        }
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}

#endif //EXP_IMS_H
