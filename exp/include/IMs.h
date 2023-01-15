//
// Created by asuka on 2022/7/16.
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
 * @brief CELF algorithm is used to select k most influential nodes at a given candidate.
 * @param graph : the graph
 * @param k : the number of nodes to be selected
 * @param candidate : the candidate node set
 * @param A : the active participant set A
 * @param seeds : returns the most influential nodes set
 */
void
CELF_local(Graph &graph, int64 k, std::vector<int64> &candidate, std::vector<int64> &A, std::vector<int64> &seeds) {
    /// first double : magimal influence
    /// first node : index
    /// second int64 : iteration round
    /// initialization
    typedef std::pair<double, std::pair<int64, int64> > node0;
    std::priority_queue<node0> Q;
    seeds.resize(1);
    double current_spread = 0;
    /// main
    for (int64 u : candidate) {
        seeds[0] = u;
        Q.push(make_pair(FI_simulation_new(graph, seeds, A), std::make_pair(u, 0)));
    }
    seeds.clear();
    while (!Q.empty() && seeds.size() < k) {
        node0 u = Q.top();
        Q.pop();
        if (u.second.second == seeds.size()) {
            seeds.emplace_back(u.second.first);
            current_spread += u.first;
        } else {
            seeds.emplace_back(u.second.first);
            u.first = FI_simulation_new(graph, seeds, A) - current_spread;
            seeds.pop_back();
            u.second.second = seeds.size();
            Q.push(u);
        }
    }
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
 * @brief Encapsulated method using local-CELF
 * @param graph : the graph
 * @param k : the number in the problem definition
 * @param A : the active participant set A
 * @param seeds : returns the seed set S (each element is a pair <node, AP>)
 */
double method_local_CELF(Graph &graph, int64 k, std::vector<int64> &A, std::vector<bi_node> &seeds) {
    assert(seeds.empty());
    std::set<int64> seeds_reorder, A_ordered(A.begin(), A.end());

    auto start_time = std::chrono::high_resolution_clock::now();
    for (int64 u : A) {
        std::vector<int64> neighbours, one_seed;
        for (auto &edge : graph.g[u]) {
            if (A_ordered.find(edge.v) == A_ordered.end()) {
                neighbours.emplace_back(edge.v);
            }
        }
        CELF_local(graph, k, neighbours, A, one_seed);
        for (int64 w : one_seed)
            if (seeds_reorder.find(w) == seeds_reorder.end()) {
                seeds_reorder.insert(w);
                seeds.emplace_back(w, u);
            }
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}

/*!
 * @brief These global variables are used to assist the recursive functions.
 * No needs to be initialized.
 */
int64 node_selected_enum[MAX_NODE_SIZE];
int64 neighbour_selected_enum[MAX_NODE_SIZE];
int64 stack_kS_enum[MAX_NODE_SIZE], stack_kS_top_enum;

/*!
 * @brief A recursive method for enumerating all possible S.
 * @param graph : the graph
 * @param S : the active participant set
 * @param V_n : a container to store all possible S
 * @param k0 : number of neighbours k of each node
 * @param k_now : Integer internal parameters, should initialize as 0
 * @param i_now : Integer internal parameters, should initialize as 0
 * @param it : Iterator internal parameters, should initialize as S.begin()
 * @param is_new : Boolean internal parameters, should initialize as true
 */
void
select_neighbours(Graph &graph, std::vector<int64> &S, std::vector<std::vector<int64> > &V_n, int64 k0, int64 k_now,
                  int64 i_now,
                  std::vector<int64>::iterator it,
                  bool is_new) {
    if (it == S.end()) {
        std::vector<int64> set_tmp;
        for (int64 i = 1; i <= stack_kS_top_enum; i++) set_tmp.emplace_back(stack_kS_enum[i]);
        V_n.emplace_back(set_tmp);
        return;
    }
    int u = *it;
    if (is_new) {
        if (it == S.begin())
            for (auto w : S) node_selected_enum[w] = 2;
        neighbour_selected_enum[u] = 0;
        for (int64 i = 0; i < graph.g[u].size(); i++) {
            int64 v = graph.g[u][i].v;
            if (node_selected_enum[v] == 1) neighbour_selected_enum[u]++;
        }
    }
    if (k_now == k0 || neighbour_selected_enum[u] == graph.g[u].size()) {
        select_neighbours(graph, S, V_n, k0, 0, 0, ++it, true);
    } else {
        for (int64 i = i_now; i < graph.g[u].size(); i++) {
            int v = graph.g[u][i].v;
            if (!node_selected_enum[v]) {
                node_selected_enum[v] = 1;
                k_now++;
                neighbour_selected_enum[u]++;
                stack_kS_enum[++stack_kS_top_enum] = v;
                select_neighbours(graph, S, V_n, k0, k_now, i + 1, it, false);
                node_selected_enum[v] = 0;
                k_now--;
                neighbour_selected_enum[u]--;
                --stack_kS_top_enum;
            }
        }
    }
    if (is_new && it == S.begin()) for (int w : S) node_selected_enum[w] = 0;
}

/*!
 * @brief Encapsulated operations for Option 1
 * @param graph : the graph
 * @param k : the number in the problem definition
 * @param A : the active participant set A
 * @param seeds : returns the seed set S = {S_1, S_2, ..., S_n}
 */
void enumeration_method(Graph &graph, int64 k, std::vector<int64> &A, std::vector<int64> &seeds) {
    std::vector<std::vector<int64> > V_n;
    std::vector<double> value;
    select_neighbours(graph, A, V_n, k, 0, 0, A.begin(), true);
    if (verbose_flag) std::cout << "The size of the solution space : " << V_n.size() << std::endl;
    for (auto &S_n : V_n) {
        if (verbose_flag) print_set(S_n, "set S = ");
        double x = FI_simulation_new(graph, S_n, A);
        value.emplace_back(x);
        if (verbose_flag) printf(" simulation result = %.5f\n", x);
    }
    double max_value = 0;
    for (int64 i = 1; i < value.size(); i++) if (value[i] > value[max_value]) max_value = i;
    seeds = V_n[max_value];
}

/*!
 * @brief Encapsulated method using greedy-PageRank
 * @param graph : the graph
 * @param k : the number in the problem definition
 * @param A : the active participant set A
 * @param seeds : returns the seed set S (each element is a pair <node, AP>)
 */
double method_greedy_PageRank(Graph &graph, int64 k, std::vector<int64> &A, std::vector<bi_node> &seeds) {
    assert(seeds.empty());
    std::vector<double> pi(graph.n, 0);
    std::vector<std::pair<double, int64> > S_ordered;
    CandidateNeigh candidate(graph, A, k);

    auto start_time = std::chrono::high_resolution_clock::now();
    power_iteration(graph, pi, 0.2);
    for (int64 w : candidate.N) S_ordered.emplace_back(std::make_pair(pi[w], w));
    //S_ordered is ordered by pageRank
    sort(S_ordered.begin(), S_ordered.end(), std::greater<>());
    for (auto &i : S_ordered) {
        int64 v = i.second;
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
 * @brief Encapsulated method using greedy-degree
 * @param graph : the graph
 * @param k : the number in the problem definition
 * @param A : the active participant set A
 * @param seeds : returns the seed set S (each element is a pair <node, AP>)
 */
double method_greedy_Degree(Graph &graph, int64 k, std::vector<int64> &A, std::vector<bi_node> &seeds) {
    assert(seeds.empty());
    std::vector<std::pair<double, int64> > S_ordered;
    CandidateNeigh candidate(graph, A, k);

    auto start_time = std::chrono::high_resolution_clock::now();
    for (int64 w : candidate.N) S_ordered.emplace_back(std::make_pair(graph.deg_out[w], w));
    //S_ordered is ordered by pageRank
    sort(S_ordered.begin(), S_ordered.end(), std::greater<>());
    for (auto &i : S_ordered) {
        int64 v = i.second;
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
        Q.push(make_pair(FI_simulation_new(graph, seeds_calc, A), std::make_pair(u, 0)));
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
            Tp.first = FI_simulation_new(graph, seeds_calc, A) - current_spread;
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
 * @brief Encapsulated method using Threshold_CELF
 * @param graph : the graph
 * @param k : the number in the problem definition
 * @param A : the active participant set A
 * @param seeds : returns the seed set S (each element is a pair <node, AP>)
 */
double method_Threshold_CELF(Graph &graph, int64 k, std::vector<int64> &A, std::vector<bi_node> &seeds) {
    assert(seeds.empty());
    double epsilon = 0.05;
    std::vector<std::pair<double, int64>> mg(graph.n);
    double W = 0;
    CandidateNeigh candidate(graph, A, k);
    std::vector<int64> seeds_calc(1); //Add a temporary space
    double current_spread = 0;

    /// main
    auto start_time = std::chrono::high_resolution_clock::now();
    for (auto u : candidate.N) {
        seeds_calc[0] = u;
        mg[u] = std::make_pair(FI_simulation_new(graph, seeds_calc, A), 0);
        W = std::max(W, mg[u].first);
    }
    seeds_calc.clear();
    for (int64 T = 1 - log(1.0 * A.size() * k / epsilon) / log(1.0 - epsilon); T > 0; T--) {
        for (auto u : candidate.N) {
            if (mg[u].second == -1) continue; //which means u was selected as seed
            if (mg[u].first < W) continue; //if f(u|S)<w, then after recalculating it also satisfies.
            int64 v = candidate.source_participant(u);
            if (v == -1) {
                mg[u].second = -1;
                continue;
            }
            if (mg[u].second < seeds.size()) { //when f(u|S) is out-of-date, update it
                seeds_calc.emplace_back(u);
                mg[u] = std::make_pair(FI_simulation_new(graph, seeds_calc, A) - current_spread, seeds.size());
                seeds_calc.pop_back();
            }
            if (mg[u].first >= W) {
                seeds.emplace_back(u, v);
                seeds_calc.emplace_back(u);
                candidate.choose(v);
                current_spread += mg[u].first;
                mg[u].second = -1; //marked that u was selected
            }
        }
        W *= 1.0 - epsilon;
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
    CandidateNeigh candidate(graph, A, k);
    std::vector<bool> node_selected(graph.n, false);
    std::vector<int64> seeds_calc(0); //Add a temporary space
    double current_spread = 0;
    int r1 = 0;
    /// main
    auto start_time = std::chrono::high_resolution_clock::now();
    while (true) {
        int64 u, v;
        double mg = -1;
        for (auto w: candidate.N) {
            if(node_selected[w]) continue;
            int64 u0 = candidate.source_participant(w);
            if (u0 == -1) {
                node_selected[w] = true;
                continue;
            }
            seeds_calc.emplace_back(w);
            double mg0 = FI_simulation_new(graph, seeds_calc, A) - current_spread;
            r1++;
            if (mg0 > mg) v = w, u = u0, mg = mg0;
            seeds_calc.pop_back();
        }
        if (mg < 0) break;

        candidate.choose(u);
        node_selected[v] = true;
        seeds.emplace_back(v, u);
        seeds_calc.emplace_back(v);
        current_spread += mg;
    }
    //printf("thr: %d\n", r1);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}

/*!
 * @brief Encapsulated method using Threshold_vanilla
 * @param graph : the graph
 * @param k : the number in the problem definition
 * @param A : the active participant set A
 * @param seeds : returns the seed set S (each element is a pair <node, AP>)
 */
double method_Threshold_vanilla(Graph &graph, int64 k, std::vector<int64> &A, std::vector<bi_node> &seeds, double epsilon = 0.05) {
    assert(seeds.empty());
    double W = 0;
    CandidateNeigh candidate(graph, A, k);
    std::vector<bool> node_selected(graph.n, false);
    std::vector<double> single_mg(graph.n, -1);
    std::vector<int64> seeds_calc(1); //Add a temporary space
    double current_spread = 0;

    /// main
    auto start_time = std::chrono::high_resolution_clock::now();
    for (auto u : candidate.N) {
        seeds_calc[0] = u;
        single_mg[u] = FI_simulation_new(graph, seeds_calc, A);
        W = std::max(W, single_mg[u]);
    }
    seeds_calc.clear();
    for (int64 T = 1 - log(1.0 * A.size() * k / epsilon) / log(1.0 - epsilon); T > 0; T--) {
        for (auto u : candidate.N) {
            if (node_selected[u]) continue; //which means u was selected as seed
            int64 v = candidate.source_participant(u);
            if (v == -1) {
                node_selected[u] = true;
                continue;
            }
            double mg;
            if(single_mg[u] > 0) {
                mg = single_mg[u];
                single_mg[u] = -1;
            } else {
                seeds_calc.emplace_back(u);
                mg = FI_simulation_new(graph, seeds_calc, A) - current_spread;
                seeds_calc.pop_back();
            }

            if (mg >= W) {
                seeds.emplace_back(u, v);
                seeds_calc.emplace_back(u);
                candidate.choose(v);
                current_spread += mg;
                node_selected[u] = true;
            }
        }
        W *= 1.0 - epsilon;
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
                    tmp_mg = FI_simulation_new(graph, seeds_calc, A);
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
                    Tp.first = FI_simulation_new(graph, seeds_calc, A) - current_influence;
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
