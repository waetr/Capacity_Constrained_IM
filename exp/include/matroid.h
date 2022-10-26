//
// Created by asuka on 2022/8/28.
//

#ifndef IMS_H_MATROID_H
#define IMS_H_MATROID_H

#include "simulation.h"

/*!
 * @brief MC-simulation-based Threshold algorithm.
 * @param graph : the graph
 * @param k : the budget
 * @param A : active participant vector
 * @param seeds : Passing parameters, returns the seed set
 * @param candidate : the structure representing the partition matroid constraint
 * @param epsilon : decrement threshold per step
 * @param seedAvgDegree
 */
void Thresholding_CELF(Graph &graph, int64 k, std::vector<int64> &A, std::vector<int64> &seeds, CandidateNeigh &candidate,
                       double epsilon, double &seedAvgDegree) {
    double cur = clock();
    auto *mg = new std::pair<double, int64>[graph.n]();
    double W = 0;
    seeds.resize(1);
    for (auto u : candidate.N) {
        seeds[0] = u;
        if (!local_mg) mg[u] = std::make_pair(FI_simulation_new(graph, seeds, A), 0);
        else mg[u] = std::make_pair(MG0[u], 0);
        W = std::max(W, mg[u].first);
    }
    seeds.clear();
    double currentSpread = 0;

    //main loop
    int r = 0;
    for (int64 T = 1 - log(1.0 * A.size() * k / epsilon) / log(1.0 - epsilon); T > 0; T--) {
        r++;
        if (verbose_flag) std::cout << "w : " << W << " round : " << r << std::endl;
        for (auto u : candidate.N) {
            if (mg[u].second == -1) continue; //which means u was selected as seed
            if (mg[u].first < W) continue; //if f(u|S)<w, then after recalculating it also satisfies.
            int64 v = candidate.source_participant(u);
            if (v == -1) {
                mg[u].second = -1;
                continue;
            }
            if (mg[u].second < seeds.size()) { //when f(u|S) is out-of-date, update it
                seeds.emplace_back(u);
                mg[u] = std::make_pair(FI_simulation_new(graph, seeds, A) - currentSpread, seeds.size());
                seeds.pop_back();
            }
            if (mg[u].first >= W) {
                if (verbose_flag) std::cout << "\tnode = " << u << "\ttime = " << time_by(cur) << std::endl;
                seeds.emplace_back(u);
                candidate.choose(v);
                currentSpread += mg[u].first;
                mg[u].second = -1; //marked that u was selected
            }
        }
        W *= 1.0 - epsilon;
    }
    if (verbose_flag) printf("Thresholding CELF1 done. total time = %.3f\n", time_by(cur));
    delete[] mg;
    seedAvgDegree = candidate.avgDegree();
}

/*!
 * @brief given RR sets, the [greedy] algorithm select the seed set based on these RR sets.
 * @param graph : the graph
 * @param A : active participant vector
 * @param k : the budget
 * @param S : Passing parameters, returns the seed set
 * @param candidate : the structure representing the partition matroid constraint
 * @param RRI : container of RR sets
 * @return the coverage of S
 */
int64 IMMSelection(Graph &graph, std::vector<int64> &A, int64 k, std::vector<bi_node> &bi_seeds, CandidateNeigh &candidate,
                    RRContainer &RRI) {
    ///temporary varible
    int64 total = 0, selected = 0, unselected = 0;
    double cur = clock(), time_rrset = 0, cur_rrset;

    std::vector<bool> RISetCovered(RRI.R.size(), false);
    for (int64 i : candidate.N) nodeRemain[i] = true;
    auto *coveredNum_tmp = new int64[graph.n];
    memcpy(coveredNum_tmp, RRI.coveredNum, graph.n * sizeof(int64));
    std::priority_queue<std::pair<int64, int64>> Q;
    for (int64 i : candidate.N) Q.push(std::make_pair(coveredNum_tmp[i], i));

    int64 current_influence = 0;
    total = candidate.N.size();
    while (!Q.empty()) {
        int64 value = Q.top().first;
        int64 maxInd = Q.top().second;
        Q.pop();

        int64 u0 = candidate.source_participant(maxInd);
        if (u0 == -1) {
            unselected++;
            continue;
        }

        if (value > coveredNum_tmp[maxInd]) {
            Q.push(std::make_pair(coveredNum_tmp[maxInd], maxInd));
            continue;
        }
        selected++;

        cur_rrset = clock();
        bi_seeds.emplace_back(maxInd, u0);
        current_influence += coveredNum_tmp[maxInd];
        candidate.choose(u0);
        nodeRemain[maxInd] = false;
        for (int64 RIIndex : RRI.covered[maxInd]) {
            if (RISetCovered[RIIndex]) continue;
            for (int64 u : RRI.R[RIIndex]) {
                if (nodeRemain[u]) coveredNum_tmp[u]--;
            }
            RISetCovered[RIIndex] = true;
        }
        time_rrset += time_by(cur_rrset);
    }
//    stdFileOut << "IMM:\n";
//    stdFileOut << "total time = " << time_by(cur) << " time select seed = " << time_rrset << std::endl;
//    stdFileOut << "total = " << total << " selected = " << selected << " unselected = " << unselected << std::endl;
    for (int64 i : candidate.N) nodeRemain[i] = false;
    delete[] coveredNum_tmp;
    return current_influence;
}

/*!
 * @brief given RR sets, the [threshold] algorithm select the seed set based on these RR sets.
 * @param graph : the graph
 * @param A : active participant vector
 * @param k : the budget
 * @param S : Passing parameters, returns the seed set
 * @param candidate : the structure representing the partition matroid constraint
 * @param RRI : container of RR sets
 * @param epsilon : decrement threshold per step
 * @return the coverage of S
 */
int64 ThresholdSelection(Graph &graph, std::vector<int64> &A, int64 k, std::vector<bi_node> &bi_seeds, CandidateNeigh &candidate,
                          double epsilon, RRContainer &RRI) {
    ///temporary varible
    int64 total = 0, selected = 0, unselected = 0;
    double cur = clock(), time_rrset = 0, cur_rrset;
    ///initialization
    std::vector<bool> RISetCovered(RRI.R.size(), false);
    for (int64 i : candidate.N) nodeRemain[i] = true;
    auto *coveredNum_tmp = new int64[graph.n];
    memcpy(coveredNum_tmp, RRI.coveredNum, graph.n * sizeof(int64));
    int64 maxi = 1 - log(1.0 * A.size() * k / epsilon) / log(1.0 - epsilon);
    auto *threshold = new double[maxi + 1];
    auto thresholdList = new std::vector<int64>[maxi + 1];
    std::vector<std::pair<int64, int64> > N_sorted;
    int64 W = 0;
    for (auto u : candidate.N) {
        W = std::max(W, coveredNum_tmp[u]);
        N_sorted.emplace_back(coveredNum_tmp[u], u);
    }
    std::sort(N_sorted.begin(), N_sorted.end(), std::greater<>());
    threshold[0] = W;
    for (int64 i = 1; i <= maxi; i++) threshold[i] = (1.0 - epsilon) * threshold[i - 1];
    for (int64 i = 0, j = 1; i < N_sorted.size(); i++) {
        while (j <= maxi && N_sorted[i].first <= threshold[j]) j++;
        if (j > maxi) break;
        thresholdList[j].emplace_back(N_sorted[i].second);
    }

    int64 current_influence = 0;
    ///main loop
    for (int i = 1; i <= maxi; i++) {
        for (auto u : thresholdList[i]) {
            int64 u0 = candidate.source_participant(u);
            if (u0 == -1) {
                unselected++;
                continue;
            }

            if (coveredNum_tmp[u] <= threshold[i]) {
                if (i + 1 <= maxi) thresholdList[i + 1].emplace_back(u);
                continue;
            }

            cur_rrset = clock();
            bi_seeds.emplace_back(u, u0);
            candidate.choose(u0);
            nodeRemain[u] = false;
            current_influence += coveredNum_tmp[u];

            selected++;
            for (int64 RIIndex : RRI.covered[u]) {
                if (RISetCovered[RIIndex]) continue;
                for (int64 v : RRI.R[RIIndex]) {
                    if (nodeRemain[v]) coveredNum_tmp[v]--;
                }
                RISetCovered[RIIndex] = true;
            }
            time_rrset += time_by(cur_rrset);
        }
    }
    total = candidate.N.size();

//    stdFileOut << "Threshold(xi = " << epsilon << "):\n";
//    stdFileOut << "total time = " << time_by(cur) << " time select seed = " << time_rrset << std::endl;
//    stdFileOut << "total = " << total << " selected = " << selected << " unselected = " << unselected << std::endl;

    for (int64 i : candidate.N) nodeRemain[i] = false;
    delete[] coveredNum_tmp;
    delete[] threshold;
    delete[] thresholdList;
    return current_influence;
}

size_t choose_from_distributionP(std::vector<double> &p) {
    std::uniform_real_distribution<double> dist(0, p[p.size() - 1]);
    double x = dist(mt19937engine);
    return std::lower_bound(p.begin(), p.end(), x) - p.begin();
}

double intPower(double a, int64 b) {
    double res = 1;
    for (; b; b >>= 1) {
        if (b & 1) res *= a;
        a *= a;
    }
    return res;
}

/*!
 * @brief given RR sets, the [prob] algorithm select the seed set based on these RR sets.
 * @param graph : the graph
 * @param A : active participant vector
 * @param k : the budget
 * @param S : Passing parameters, returns the seed set
 * @param alpha : a valid upperbound of the curvature of spread function
 * @param RRI : container of RR sets
 * @return the coverage of S
 */
int64 prob(Graph &graph, std::vector<int64> &A, int64 k, std::vector<bi_node> &bi_seeds, double alpha,
           RRContainer &RRI) {
    ///temporary varible
    double cur = clock(), time_rrset = 0, cur_rrset;
    int64 total = 0, selected = 0, unselected = 0;
    std::set<int64> N_unselected;
    ///initialization
    std::set<int64> N[A.size()];
    for (int i = 0; i < A.size(); i++) {
        for (auto e : graph.g[A[i]])
            if (std::find(A.begin(), A.end(), e.v) == A.end()) {
                N_unselected.insert(e.v); //Unimportant statements
                N[i].insert(e.v);
                nodeRemain[e.v] = true;
            }
    }
    total = N_unselected.size(); //Unimportant statements
    std::vector<bool> RISetCovered(RRI.numOfRRsets(), false);
    auto *coveredNum_tmp = new int64[graph.n];
    memcpy(coveredNum_tmp, RRI.coveredNum, graph.n * sizeof(int64));

    ///main loop~
    int64 current_influence = 0, num_of_filled_N = 0;
    std::vector<int64> N_numbers(A.size(), 0); ///N_numbers[i] == k + 1 means that N[i] is full
    while (num_of_filled_N < A.size()) {
        for (int i = 0; i < A.size(); i++) {
            if (N_numbers[i] != k + 1 && (N_numbers[i] == k || N[i].empty())) {
                N_numbers[i] = k + 1;
                num_of_filled_N++;
            }
            if (N_numbers[i] == k + 1) continue;
            std::vector<int64> V;
            std::vector<double> p(N[i].size());
            V.assign(N[i].begin(), N[i].end()); //convert N(set) into V(vector)
            int64 alpha0 = ceil((1.0 + V.size()) / alpha) - 1;
            double p0 = 1;
            for (long j : V) p0 += 1.0 * coveredNum_tmp[j] / V.size();
            for (size_t j = 0; j < V.size(); j++) { //Dividing by p0 is to prevent the result from overflowing
                p[j] = ((j == 0) ? 0 : p[j - 1]) + intPower(1.0 * coveredNum_tmp[V[j]] / p0, alpha0);
            }
            auto v = V[choose_from_distributionP(p)];
            ///choose v
            cur_rrset = clock(); //Unimportant statements
            selected++;
            N_numbers[i]++;
            bi_seeds.emplace_back(v, A[i]);
            current_influence += coveredNum_tmp[v];
            nodeRemain[v] = false;
            for (int64 RIIndex : RRI.covered[v]) {
                if (RISetCovered[RIIndex]) continue;
                for (int64 u : RRI.R[RIIndex]) {
                    if (nodeRemain[u]) coveredNum_tmp[u]--;
                }
                RISetCovered[RIIndex] = true;
            }
            time_rrset += time_by(cur_rrset); //Unimportant statements
            ///delete v from all N[i]
            for (int j = 0; j < A.size(); j++) {
                N[j].erase(v);
            }
        }
    }
    //Unimportant statements
    N_unselected.clear();
    for (int i = 0; i < A.size(); i++) {
        for (auto u : N[i])
            N_unselected.insert(u);
    }
    unselected = N_unselected.size();

    for (long i : A)
        for (auto e : graph.g[i]) nodeRemain[e.v] = false;
    delete[] coveredNum_tmp;

//    stdFileOut << "total time = " << time_by(cur) << " time select seed = " << time_rrset << std::endl;
//    stdFileOut << "total = " << total << " selected = " << selected << " unselected = " << unselected << std::endl;
    return current_influence;
}

/*!
 * @brief given RR sets, the [d-prob] algorithm select the seed set based on these RR sets.
 * @param graph : the graph
 * @param A : active participant vector
 * @param k : the budget
 * @param S : Passing parameters, returns the seed set
 * @param RRI : container of RR sets
 * @return the coverage of S
 */
int64 prob_determined(Graph &graph, std::vector<int64> &A, int64 k, std::vector<bi_node> &bi_seeds, RRContainer &RRI) {
    ///temporary varible
    double cur = clock(), time_rrset = 0, cur_rrset;
    int64 total = 0, selected = 0, unselected = 0;
    std::set<int64> N_unselected;
    ///initialization
    std::set<int64> N[A.size()];
    for (int i = 0; i < A.size(); i++) {
        for (auto e : graph.g[A[i]])
            if (std::find(A.begin(), A.end(), e.v) == A.end()) {
                N_unselected.insert(e.v);
                N[i].insert(e.v);
                nodeRemain[e.v] = true;
            }
    }
    total = N_unselected.size();
    std::vector<bool> RISetCovered(RRI.R.size(), false);
    auto *coveredNum_tmp = new int64[graph.n];
    memcpy(coveredNum_tmp, RRI.coveredNum, graph.n * sizeof(int64));

    int64 current_influence = 0, N_empty = 0;
    std::vector<int64> Ni_empty(A.size(), 0);
    while (N_empty < A.size()) {
        for (int i = 0; i < A.size(); i++) { ///N_numbers[i] == k + 1 means that N[i] is full
            if (Ni_empty[i] != k + 1 && (Ni_empty[i] == k || N[i].empty())) {
                Ni_empty[i] = k + 1;
                N_empty++;
            }
            if (Ni_empty[i] == k + 1) continue;
            int64 v = -1;
            for(auto e : N[i])
                if(v == -1 || coveredNum_tmp[e] > coveredNum_tmp[v]) v = e;
            ///choose v
            cur_rrset = clock();
            selected++;
            Ni_empty[i]++;
            bi_seeds.emplace_back(v, A[i]);
            current_influence += coveredNum_tmp[v];
            nodeRemain[v] = false;
            for (int64 RIIndex : RRI.covered[v]) {
                if (RISetCovered[RIIndex]) continue;
                for (int64 u : RRI.R[RIIndex]) {
                    if (nodeRemain[u]) coveredNum_tmp[u]--;
                }
                RISetCovered[RIIndex] = true;
            }
            time_rrset += time_by(cur_rrset);
            ///delete v from all N[i]
            for (int j = 0; j < A.size(); j++) {
                N[j].erase(v);
            }
        }
    }

    N_unselected.clear();
    for (int i = 0; i < A.size(); i++) {
        for (auto u : N[i])
            N_unselected.insert(u);
    }

    unselected = N_unselected.size();

    for (long i : A)
        for (auto e : graph.g[i]) nodeRemain[e.v] = false;
    delete[] coveredNum_tmp;

//    stdFileOut << "total time = " << time_by(cur) << " time select seed = " << time_rrset << std::endl;
//    stdFileOut << "total = " << total << " selected = " << selected << " unselected = " << unselected << std::endl;
    return current_influence;
}

/*!
 * @brief a valid upperbound of the curvature of spread function
 * @param graph : the graph
 * @param A : ap set
 * @param RRI : RR sets
 * @return : a decimal, representing the upperbound
 */
double Alpha_upperbound(Graph &graph, std::vector<int64> &A, RRContainer &RRI) {
    std::set<int64> N;
    for (auto u : A)
        for (auto e : graph.g[u])
            if (std::find(A.begin(), A.end(), e.v) == A.end()) N.insert(e.v);
    std::vector<int64> RISetCovered(RRI.numOfRRsets(), 0), lastCovered(RRI.numOfRRsets(), 0);
    std::vector<int64> mg(graph.n, 0);
    for (auto u : N) {
        for (int64 RIIndex : RRI.covered[u]) {
            RISetCovered[RIIndex]++;
            lastCovered[RIIndex] = u;
        }
    }
    for (int64 i = 0; i < RRI.numOfRRsets(); i++) {
        if (RISetCovered[i] == 1) mg[lastCovered[i]]++;
    }
    double ans = 1.0;
    for (auto u : N) {
        if (RRI.coveredNum[u] == 0) continue;
        ans = std::min(ans, 1.0 * mg[u] / RRI.coveredNum[u]);
    }
    return 1.0 - ans;
}

#if false
/*!
 * @brief Explicit representation of the partition matroid.
 * The representation of the symbols in partition matroid is based on the paper "Comparing Apples and Oranges".
 * @param h : the number of partitions
 * @param k : the maximum number of elements that can be taken from a partition(i.e. numbers of nodes taken by an AP)
 * @param rank : rank of the partition matroid(=h*k)
 * @param N : N[i] stores the elements in i-th partition [i = 1,...,h]
 */
class partitionMatroid {
public:
    int32 h, k, rank;
    std::vector<node> *N;

    partitionMatroid(Graph &graph, std::vector<node> &A, int32 k0) {
        h = A.size();
        k = k0;
        rank = h * k;
        N = new std::vector<node>[A.size() + 1]();
        std::vector<node> f[graph.n]; //f[v] means in-coming active participant of v
        std::set<node> N_unique, A_unique;
        for (auto u : A) A_unique.insert(u);
        //update f
        for (int i = 0; i < A.size(); i++) {
            node u = A[i];
            for (auto &edge : graph.g[u]) {
                if (A_unique.find(edge.v) == A_unique.end()) {
                    f[edge.v].emplace_back(i + 1);
                    N_unique.insert(edge.v);
                }
            }
        }
        //update elements
        for (auto v : N_unique) {
            std::uniform_int_distribution<size_t> intDistribution(0, f[v].size() - 1);
            size_t x = intDistribution(random_engine);

            N[f[v][x]].emplace_back(v);
        }
    }

    ~partitionMatroid() {
        delete[] N;
    }
};

void Thresholding_CELF(Graph &graph, int32 k, std::vector<node> &A, std::vector<node> &seeds) {
    double epsilon = 0.05;
    double cur = clock();

    partitionMatroid matroid(graph, A, k);
    auto *mg = new std::pair<double, int32>[graph.n]();
    auto N_selected = new size_t[A.size() + 1]();
    double W = 0;
    seeds.resize(1);
    for (int32 i = 1; i <= matroid.h; i++) {
        for (auto u : matroid.N[i]) {
            seeds[0] = u;
            if (!local_mg) mg[u] = std::make_pair(FI_simulation_new(graph, seeds, A), 0);
            else mg[u] = std::make_pair(MG0[u], 0);
            W = std::max(W, mg[u].first);
        }
    }
    seeds.clear();
    double currentSpread = 0;

    int r = 0;
    //main loop
    for (int64 T = 1 - log(1.0 * matroid.rank / epsilon) / log(1.0 - epsilon); T > 0; T--) {
        r++;
        if (verbose_flag) std::cout << "w : " << W << " round : " << r << std::endl;
        for (int32 i = 1; i <= matroid.h; i++) {
            if (N_selected[i] == matroid.k) continue;
            for (auto u : matroid.N[i]) {
                if (mg[u].second == -1) continue; //which means u was selected as seed
                if (mg[u].first < W) continue; //if f(u|S)<w, then after recalculating it also satisfies.
                if (mg[u].second < seeds.size()) { //when f(u|S) is out-of-date, update it
                    seeds.emplace_back(u);
                    mg[u] = std::make_pair(FI_simulation_new(graph, seeds, A) - currentSpread, seeds.size());
                    seeds.pop_back();
                }
                if (mg[u].first >= W) {
                    if (verbose_flag) std::cout << "\tnode = " << u << "\ttime = " << time_by(cur) << std::endl;
                    seeds.emplace_back(u);
                    N_selected[i]++;
                    currentSpread += mg[u].first;
                    mg[u].second = -1; //marked that u was selected
                    if (N_selected[i] == matroid.k) break;
                }
            }
            W *= 1.0 - epsilon;
        }
    }
    if (verbose_flag) printf("Thresholding CELF done. total time = %.3f\n", time_by(cur));
    delete[] mg;
    delete[] N_selected;
}

void Thresholding_CELF_unused(Graph &graph, int32 k, std::vector<node> &A, std::vector<node> &seeds) {
    double epsilon = 0.05;
    double cur = clock();

    CandidateNeigh candidate(graph, A, k);
    typedef std::pair<double, std::pair<node, int64> > node0;
    std::priority_queue<node0> Q;
    double W = 0;
    seeds.resize(1);
    for (auto u : candidate.N) {
        seeds[0] = u;
        double mg;
        if (!local_mg) mg = FI_simulation(graph, seeds);
        else mg = MG0[u];
        W = std::max(W, mg);
        Q.push(make_pair(mg, std::make_pair(u, 0)));
    }
    seeds.clear();
    double currentSpread = 0;

    //main loop
    for (double w = W; w > epsilon * W / A.size() / k; w *= 1.0 - epsilon) {
        if (verbose_flag) std::cout << "w : " << w << std::endl;
        if (Q.empty()) break;
        while (!Q.empty()) {
            node0 Tp = Q.top();
            node v = Tp.second.first;
            int64 it_round = Tp.second.second;
            double mg = Tp.first;
            if (mg < w) break; //if f(v|S)<w, then after recalculating it also satisfies.
            Q.pop();
            node u = candidate.source_participant(v);
            if (u == -1) continue;
            if (it_round < seeds.size()) { //when f(v|S) is out-of-date, update it and reinsert into Q
                seeds.emplace_back(v);
                Tp.first = FI_simulation(graph, seeds) - currentSpread;
                seeds.pop_back();
                Tp.second.second = seeds.size();
                Q.push(Tp);
            } else if (mg >= w) {
                if (verbose_flag) std::cout << "\tnode = " << u << "\ttime = " << time_by(cur) << std::endl;
                candidate.choose(u);
                seeds.emplace_back(v);
                currentSpread += mg;
            }
        }
    }
    if (verbose_flag) printf("Thresholding CELF2 done. total time = %.3f\n", time_by(cur));
}

double SuperThresholdSelection(Graph &graph, std::vector<node> &A, int32 k, std::vector<node> &S, CandidateNeigh &candidate,
                          double epsilon, RRContainer &RRI) {
    ///temporary varible
    int64 total = 0, selected = 0, unselected = 0;
    double cur = clock(), time_rrset = 0, cur_rrset;
    ///initialization
    S.clear();
    std::vector<bool> RISetCovered(RRI.R.size(), false);
    auto *coveredNum_tmp = new int64[graph.n];
    memcpy(coveredNum_tmp, RRI.coveredNum, graph.n * sizeof(int64));
    int64 maxi = 1 - log(1.0 * A.size() * k / epsilon) / log(1.0 - epsilon);
    for (node i : candidate.N) nodeRemain[i] = true;
    std::set<node> NN(candidate.N.begin(), candidate.N.end());


    while(1) {
        ///initialize
        CandidateNeigh candidate0;
        std::vector<node> M;
        std::vector<std::pair<int64, node> > N_sorted;
        int64 W = 0;
        auto *threshold = new double[maxi + 1];
        auto thresholdList = new std::vector<node>[maxi + 1];
        auto *coveredNum_tmp_tmp = new int64[graph.n];
        auto *nodeRemain_ = new bool[graph.n];
        for (node i : NN) nodeRemain[i] = true;
        std::vector<bool> RISetCovered_(RRI.R.size(), false);
        memcpy(coveredNum_tmp_tmp, coveredNum_tmp, graph.n * sizeof(int64));
        for(int i = 0; i <= maxi; i++) thresholdList[i].clear();
        N_sorted.clear();
        candidate0.assign(graph.n, candidate);
        ///main part
        for (auto u : NN) {
            W = std::max(W, coveredNum_tmp_tmp[u]);
            N_sorted.emplace_back(coveredNum_tmp_tmp[u], u);
        }

        std::sort(N_sorted.begin(), N_sorted.end(), std::greater<>());

        threshold[0] = W;
        for (int64 i = 1; i <= maxi; i++) threshold[i] = (1.0 - epsilon) * threshold[i - 1];

        for (int64 i = 0, j = 1; i < N_sorted.size(); i++) {
            while (j <= maxi && N_sorted[i].first <= threshold[j]) j++;
            if(j > maxi) break;
            thresholdList[j].emplace_back(N_sorted[i].second);
        }

        ///main loop

        for (int i = 1; i <= maxi; i++) {
            for (auto u : thresholdList[i]) {
                node u0 = candidate0.source_participant(u);
                if (u0 == -1) {
                    unselected++;
                    continue;
                }

                if (coveredNum_tmp_tmp[u] <= threshold[i]) {
                    if (i + 1 <= maxi) thresholdList[i + 1].emplace_back(u);
                    continue;
                }

                cur_rrset = clock();
                M.emplace_back(u);
                candidate0.choose(u0);
                nodeRemain_[u] = false;

                selected++;
                for (int64 RIIndex : RRI.covered[u]) {
                    if (RISetCovered_[RIIndex]) continue;
                    for (node v : RRI.R[RIIndex]) {
                        if (nodeRemain_[v]) coveredNum_tmp_tmp[v]--;
                    }
                    RISetCovered_[RIIndex] = true;
                }
                time_rrset += time_by(cur_rrset);
            }
        }
        total += NN.size();

        delete[] threshold;
        delete[] thresholdList;
        delete[] coveredNum_tmp_tmp;

        ///really choose one into S
        if(M.empty()) break;
        auto x = M[mt19937engine() % M.size()];
        S.emplace_back(x);
        int y = candidate.source_participant(x);
        if(y == -1) continue;
        candidate.choose(y);
        for (int64 RIIndex : RRI.covered[x]) {
            if (RISetCovered[RIIndex]) continue;
            for (node v : RRI.R[RIIndex]) {
                if (nodeRemain[v]) coveredNum_tmp[v]--;
            }
            RISetCovered[RIIndex] = true;
        }
        NN.erase(x);
    }

    for (node i : candidate.N) nodeRemain[i] = false;
    stdFileOut << "Threshold(xi = " << epsilon << "):\n";
    stdFileOut << "total time = " << time_by(cur) << " time select seed = " << time_rrset << std::endl;
    stdFileOut << "total = " << total << " selected = " << selected << " unselected = " << unselected << std::endl;
    stdFileOut << "seed of Threshold = " << 1.0 * RRI.self_inf_cal(graph, S) / RRI.numOfRRsets() * graph.n;
    stdFileOut << "\tsize = " << S.size() << std::endl;

    delete[] coveredNum_tmp;

    return 0;
}

int64 OPIMSelection(Graph &graph, std::vector<node> &A, int32 k, std::vector<node> &S, CandidateNeigh &candidate,
                    RRContainer &RRI) {
    ///temporary varible
    int64 total = 0, selected = 0, unselected = 0;
    double cur = clock(), time_rrset = 0, cur_rrset;

    int32 kA = 0;
    for (node u : A) kA += std::min(k, (int32) graph.g[u].size());
    kA = std::min(kA, (int32) graph.n);
    S.clear();
    std::vector<bool> RISetCovered(RRI.R.size(), false);
    for (node i : candidate.N) nodeRemain[i] = true;
    auto *coveredNum_tmp = new int64[graph.n];
    memcpy(coveredNum_tmp, RRI.coveredNum, graph.n * sizeof(int64));
    std::vector<node> coveredNum_maxK;
    int64 coverage = 0, A_u = INT64_MAX;

    total = candidate.N.size();
    for (int32 i = 0;; i++) {
        if (i > 0) {
            node maxInd = -1;
            for (node v : candidate.N) {
                if (!nodeRemain[v]) continue;
                node u = candidate.source_participant(v);
                if (u == -1) {
                    nodeRemain[v] = false;
                    unselected++;
                    continue;
                }
                if (maxInd == -1 || coveredNum_tmp[v] > coveredNum_tmp[maxInd]) maxInd = v;
            }
            if (maxInd == -1) break;
            node u0 = candidate.source_participant(maxInd);

            selected++;
            cur_rrset = clock();
            coverage += coveredNum_tmp[maxInd];
            S.emplace_back(maxInd);
            candidate.choose(u0);
            nodeRemain[maxInd] = false;
            for (int64 RIIndex : RRI.covered[maxInd]) {
                if (RISetCovered[RIIndex]) continue;
                for (node u : RRI.R[RIIndex]) {
                    if (nodeRemain[u]) coveredNum_tmp[u]--;
                }
                RISetCovered[RIIndex] = true;
            }
            time_rrset += time_by(cur_rrset);
        }
        //find the max-k nodes with maximum marginal coverage
        int64 A_u_i = coverage;
        coveredNum_maxK.clear();
        for (node u : candidate.N) {
            if (nodeRemain[u])
                coveredNum_maxK.emplace_back(coveredNum_tmp[u]);
        }
        if (coveredNum_maxK.size() >= kA) {
            nth_element(coveredNum_maxK.begin(), coveredNum_maxK.end() - kA, coveredNum_maxK.end());
        }
        size_t begin0 = coveredNum_maxK.size() <= kA ? 0 : coveredNum_maxK.size() - kA;
        for (node j = begin0; j < coveredNum_maxK.size(); j++) {
            A_u_i += coveredNum_maxK[j];
        }
        A_u = std::min(A_u, A_u_i);
    }

    stdFileOut << "OPIM:\n";
    stdFileOut << "total time = " << time_by(cur) << " time select seed = " << time_rrset << std::endl;
    stdFileOut << "total = " << total << " selected = " << selected << " unselected = " << unselected << std::endl;
    stdFileOut << "seed of OPIM = " << 1.0 * RRI.self_inf_cal(graph, S) / RRI.numOfRRsets() * graph.n;
    stdFileOut << "\tsize = " << S.size() << std::endl;

    delete[] coveredNum_tmp;
    for (node i : candidate.N) nodeRemain[i] = false;
    return A_u;
}
#endif

#endif //IMS_H_MATROID_H
