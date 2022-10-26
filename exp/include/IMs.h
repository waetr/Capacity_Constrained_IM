//
// Created by asuka on 2022/7/16.
// IM Solvers

#ifndef EXP_IMS_H
#define EXP_IMS_H

#include "simulation.h"


/*!
 * @brief MG0[u] stores influence spread of {u}
 */
double MG0[MAX_NODE_SIZE];

/*!
 * @brief CELF algorithm is used to select k most influential nodes at a given candidate.
 * @param graph : the graph
 * @param k : the number of nodes to be selected
 * @param candidate : the candidate node set
 * @param A : the active participant set A
 * @param seeds : returns the most influential nodes set
 */
void CELF(Graph &graph, int64 k, std::vector<int64> &candidate, std::vector<int64> &A, std::vector<int64> &seeds){
    if (k >= candidate.size()) {
        seeds = candidate;
        if (verbose_flag) printf("Nodes are not exceeding k. All selected.\n");
        return;
    }
    double cur = clock();
    int64 r = 0;
    /// first double : magimal influence
    /// first node : index
    /// second int64 : iteration round
    typedef std::pair<double, std::pair<int64, int64> > node0;
    std::priority_queue<node0> Q;
    seeds.resize(1);
    for (int64 u : candidate) {
        seeds[0] = u;
        if (!local_mg) Q.push(make_pair(FI_simulation_new(graph, seeds, A), std::make_pair(u, 0)));
        else Q.push(make_pair(MG0[u], std::make_pair(u, 0)));
    }
    double current_spread = 0;
    seeds.clear();
    if (verbose_flag) printf("\tInitialization time = %.5f\n", time_by(cur));
    while (seeds.size() < k) {
        r++;
        node0 u = Q.top();
        Q.pop();
        if (u.second.second == seeds.size()) {
            if (verbose_flag) {
                std::cout << "\tnode = " << u.second.first << "\tround = " << r << "\ttime = " << time_by(cur) << std::endl;
            }
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
    if (verbose_flag) printf("CELF done. total time = %.3f\n", time_by(cur));
}

/*!
 * @brief Use power iteration method to calculate pagerank values of nodes in graph.
 * @remarks Note that if there are isolated nodes in the graph, the iteration will not stop causing an error.
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
 * @brief Encapsulated operations for Option 2 using IM solver : pagerank
 * @param graph : the graph
 * @param k : the number in the problem definition
 * @param A : the active participant set A
 * @param seeds : returns the seed set S = {S_1, S_2, ..., S_n}
 */
void pgrank_method(Graph &graph, int64 k, std::vector<int64> &A, std::vector<int64> &seeds) {
    std::vector<double> pi(graph.n, 0);
    power_iteration(graph, pi, 0.2);
    std::vector<std::pair<double, int64> > pg_rank;
    std::set<int64> seeds_reorder;
    for (int64 u : A) {
        pg_rank.clear();
        for (auto &edge : graph.g[u]) {
            int64 v = edge.v;
            if (find(A.begin(), A.end(), v) == A.end()) {
                pg_rank.emplace_back(std::make_pair(pi[v], v));
            }
        }
        sort(pg_rank.begin(), pg_rank.end(), std::greater<>());
        for (int64 i = 0; i < k && i < pg_rank.size(); i++)
            seeds_reorder.insert(pg_rank[i].second);
    }
    for (int64 w : seeds_reorder) seeds.emplace_back(w);
    seeds_reorder.clear();
    pg_rank.clear();
}

/*!
 * @brief Encapsulated operations for Option 2 using IM solver : degree
 * @param graph : the graph
 * @param k : the number in the problem definition
 * @param A : the active participant set A
 * @param seeds : returns the seed set S = {S_1, S_2, ..., S_n}
 */
void degree_method(Graph &graph, int64 k, std::vector<int64> &A, std::vector<int64> &seeds) {
    std::vector<std::pair<int64, int64> > degree_rank;
    std::set<int64> seeds_reorder;
    for (int64 u : A) {
        degree_rank.clear();
        for (auto &edge : graph.g[u]) {
            int64 v = edge.v;
            if (find(A.begin(), A.end(), v) == A.end()) {
                degree_rank.emplace_back(std::make_pair(graph.deg_out[v], v));
            }
        }
        sort(degree_rank.begin(), degree_rank.end(), std::greater<>());
        for (int64 i = 0; i < k && i < degree_rank.size(); i++)
            seeds_reorder.insert(degree_rank[i].second);
    }
    for (int64 w : seeds_reorder) seeds.emplace_back(w);
    seeds_reorder.clear();
    degree_rank.clear();
}

/*!
 * @brief Encapsulated operations for Option 2 using IM solver : CELF
 * @param graph : the graph
 * @param k : the number in the problem definition
 * @param A : the active participant set A
 * @param seeds : returns the seed set S = {S_1, S_2, ..., S_n}
 */
void CELF_method(Graph &graph, int64 k, std::vector<int64> &A, std::vector<int64> &seeds) {
    double cur = clock();
    std::set<int64> seeds_reorder;
    for (int64 u : A) {
        std::vector<int64> neighbours, one_seed;
        for (auto &edge : graph.g[u]) {
            if (find(A.begin(), A.end(), edge.v) == A.end()) {
                neighbours.emplace_back(edge.v);
            }
        }
        CELF(graph, k, neighbours, A, one_seed);
        for (int64 w : one_seed)
            seeds_reorder.insert(w);
    }
    for (int64 w : seeds_reorder) seeds.emplace_back(w);
    seeds_reorder.clear();
    if (verbose_flag) printf("CELF method done. total time = %.3f\n", time_by(cur));
}

/*!
 * @brief These global variables are used to assist the recursive functions.
 * No needs to be initialized.
 */
int64 node_selected[MAX_NODE_SIZE];
int64 neighbour_selected[MAX_NODE_SIZE];
int64 stack_kS[MAX_NODE_SIZE], stack_kS_top;

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
void select_neighbours(Graph &graph, std::vector<int64> &S, std::vector<std::vector<int64> > &V_n, int64 k0, int64 k_now, int64 i_now,
                       std::vector<int64>::iterator it,
                       bool is_new) {
    if (it == S.end()) {
        std::vector<int64> set_tmp;
        for (int64 i = 1; i <= stack_kS_top; i++) set_tmp.emplace_back(stack_kS[i]);
        V_n.emplace_back(set_tmp);
        return;
    }
    int u = *it;
    if (is_new) {
        if (it == S.begin())
            for (auto w : S) node_selected[w] = 2;
        neighbour_selected[u] = 0;
        for (int64 i = 0; i < graph.g[u].size(); i++) {
            int64 v = graph.g[u][i].v;
            if (node_selected[v] == 1) neighbour_selected[u]++;
        }
    }
    if (k_now == k0 || neighbour_selected[u] == graph.g[u].size()) {
        select_neighbours(graph, S, V_n, k0, 0, 0, ++it, true);
    } else {
        for (int64 i = i_now; i < graph.g[u].size(); i++) {
            int v = graph.g[u][i].v;
            if (!node_selected[v]) {
                node_selected[v] = 1;
                k_now++;
                neighbour_selected[u]++;
                stack_kS[++stack_kS_top] = v;
                select_neighbours(graph, S, V_n, k0, k_now, i + 1, it, false);
                node_selected[v] = 0;
                k_now--;
                neighbour_selected[u]--;
                --stack_kS_top;
            }
        }
    }
    if (is_new && it == S.begin()) for (int w : S) node_selected[w] = 0;
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
 * @brief Encapsulated operations for advanced version of IM solver : pagerank
 * @param graph : the graph
 * @param k : the number in the problem definition
 * @param A : the active participant set A
 * @param seeds : returns the seed set S = {S_1, S_2, ..., S_n}
 */
void advanced_pgrank_method(Graph &graph, int64 k, std::vector<int64> &A, std::vector<int64> &seeds, double &seedAvgDegree) {
    std::vector<double> pi(graph.n, 0);
    power_iteration(graph, pi, 0.2);
    std::vector<std::pair<double, int64> > S_ordered;
    CandidateNeigh candidate(graph, A, k);
    for (int64 w : candidate.N) S_ordered.emplace_back(std::make_pair(pi[w], w));
    //S_ordered is ordered by pageRank
    sort(S_ordered.begin(), S_ordered.end(), std::greater<>());
    for (auto & i : S_ordered) {
        int64 v = i.second;
        int64 u = candidate.source_participant(v);
        if (u == -1) continue;
        candidate.choose(u);
        seeds.emplace_back(v);
    }
    seedAvgDegree = candidate.avgDegree();
}

/*!
 * @brief Encapsulated operations for advanced version of IM solver : degree
 * @param graph : the graph
 * @param k : the number in the problem definition
 * @param A : the active participant set A
 * @param seeds : returns the seed set S = {S_1, S_2, ..., S_n}
 */
void advanced_degree_method(Graph &graph, int64 k, std::vector<int64> &A, std::vector<int64> &seeds, double &seedAvgDegree) {
    std::vector<std::pair<double, int64> > S_ordered;
    CandidateNeigh candidate(graph, A, k);
    for (int64 w : candidate.N) S_ordered.emplace_back(std::make_pair(graph.deg_out[w], w));
    //S_ordered is ordered by degree
    sort(S_ordered.begin(), S_ordered.end(), std::greater<>());
    for (auto & i : S_ordered) {
        int64 v = i.second;
        int64 u = candidate.source_participant(v);
        if (u == -1) continue;
        candidate.choose(u);
        seeds.emplace_back(v);
    }
    seedAvgDegree = candidate.avgDegree();
}

/*!
 * @brief Encapsulated operations for advanced version of IM solver : CELF
 * @param graph : the graph
 * @param k : the number in the problem definition
 * @param A : the active participant set A
 * @param seeds : returns the seed set S = {S_1, S_2, ..., S_n}
 */
void advanced_CELF_method(Graph &graph, int64 k, std::vector<int64> &A, std::vector<int64> &seeds, double &seedAvgDegree) {
    double cur = clock();
    int64 r = 0;
    CandidateNeigh candidate(graph, A, k);
    typedef std::pair<double, std::pair<int64, int64> > node0;
    std::priority_queue<node0> Q;
    seeds.resize(1); //Add a temporary space
    for (int64 u : candidate.N) {
        seeds[0] = u;
        if (!local_mg) {
            Q.push(make_pair(FI_simulation_new(graph, seeds, A), std::make_pair(u, 0)));
        } else {
            Q.push(make_pair(MG0[u], std::make_pair(u, 0)));
        }
    }
    double current_spread = 0;
    seeds.pop_back(); //Clear the temporary space
    if (verbose_flag) printf("\tInitialization time = %.5f\n", time_by(cur));
    while (!Q.empty()) {
        node0 Tp = Q.top();
        int64 v = Tp.second.first;
        int64 it_round = Tp.second.second;
        double mg = Tp.first;
        Q.pop();
        int64 u = candidate.source_participant(v);
        if (u == -1) continue;
        r++;
        if (it_round == seeds.size()) {
            candidate.choose(u);
            seeds.emplace_back(v);
            current_spread += mg;
            if (verbose_flag) {
                std::cout << "\tnode = " << v << "\tround = " << r << "\ttime = " << time_by(cur) << std::endl;
            }
        } else {
            seeds.emplace_back(v);
            Tp.first = FI_simulation_new(graph, seeds, A) - current_spread;
            seeds.pop_back();
            Tp.second.second = seeds.size();
            Q.push(Tp);
        }
    }
    seedAvgDegree = candidate.avgDegree();
    if (verbose_flag) printf("CELF advanced done. total time = %.3f\n", time_by(cur));
}

#endif //EXP_IMS_H
