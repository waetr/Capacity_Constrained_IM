//
// Created by asuka on 2022/11/17.
//

#ifndef EXP_OPIM_NEW_H
#define EXP_OPIM_NEW_H

/*!
 * @brief given RR sets, the [d-prob] algorithm select the seed set based on these RR sets.
 * @param graph : the graph
 * @param A : active participant vector
 * @param k : the budget
 * @param S : Passing parameters, returns the seed set
 * @param RRI : container of RR sets
 * @return the coverage of S
 */
int64
RR_OPIM_Selection(Graph &graph, std::vector<int64> &A, int64 k, std::vector<bi_node> &bi_seeds, RRContainer &RRI,
                  bool is_tightened) {
    ///temporary varible
    std::vector<bool> nodeRemain(graph.n, false);
    ///initialization
    std::set<int64> A_reorder(A.begin(), A.end());
    for (int i = 0; i < A.size(); i++) {
        for (auto e : graph.g[A[i]])
            if (A_reorder.find(e.v) == A_reorder.end()) {
                nodeRemain[e.v] = true;
            }
    }
    std::vector<bool> RISetCovered(RRI.R.size(), false);
    auto *coveredNum_tmp = new int64[graph.n];
    memcpy(coveredNum_tmp, RRI.coveredNum, graph.n * sizeof(int64));

    int64 current_influence = 0, N_empty = 0;

    int64 xx = INT64_MAX;

    std::vector<int64> Ni_empty(A.size(), 0);
    while (N_empty < A.size()) {
        if (is_tightened) {
            int64 tight_mg = 0;
            std::vector<int64> tight_list;
            for (int i = 0; i < A.size(); i++) {
                tight_list.clear();
                for (auto e : graph.g[A[i]]) {
                    if (coveredNum_tmp[e.v] > 0)
                        tight_list.emplace_back(coveredNum_tmp[e.v]);
                }
                int64 k_max = std::min((int64) tight_list.size(), k);
                std::nth_element(tight_list.begin(), tight_list.begin() + k_max - 1, tight_list.end(),
                                 std::greater<>());
                for (int j = 0; j < k_max; j++) {
                    tight_mg += tight_list[j];
                }
            }
            xx = std::min(xx, current_influence + tight_mg);
        }
        for (int i = 0; i < A.size(); i++) { ///N_numbers[i] == k + 1 means that N[i] is full
            if (Ni_empty[i] != k + 1 && Ni_empty[i] == k) {
                Ni_empty[i] = k + 1;
                N_empty++;
            }
            if (Ni_empty[i] == k + 1) continue;

            int64 v = -1;
            for (auto e : graph.g[A[i]])
                if (nodeRemain[e.v] && (v == -1 || coveredNum_tmp[e.v] > coveredNum_tmp[v])) v = e.v;

            if (Ni_empty[i] != k + 1 && v == -1) {
                Ni_empty[i] = k + 1;
                N_empty++;
            }
            if (Ni_empty[i] == k + 1) continue;
            ///choose v
            Ni_empty[i]++;
            bi_seeds.emplace_back(v, A[i]);
            current_influence += coveredNum_tmp[v];
            nodeRemain[v] = false;
            for (int64 RIIndex : RRI.covered[v]) {
                if (RISetCovered[RIIndex]) continue;
                for (int64 u : RRI.R[RIIndex]) {
                    coveredNum_tmp[u]--;
                }
                RISetCovered[RIIndex] = true;
            }
        }
    }
    delete[] coveredNum_tmp;
    return is_tightened ? xx : current_influence;
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
void
RR_OPIM_Main(Graph &graph, std::vector<int64> &A, int64 k, double eps, double delta, std::vector<bi_node> &bi_seeds,
             bool is_tightened = false) {
    double timer = 0, timer0 = 0, cur;
    double sum_log = 0;
    for (auto ap:A) sum_log += logcnk(graph.deg_out[ap], k);
    double C_0 = 2.0 * sqr(
            0.5 * sqrt(log(6.0 / delta)) + sqrt(0.5 * (sum_log + log(6.0 / delta))));
    RRContainer R1(graph, A, true), R2(graph, A, true);
    cur = clock();
    R1.resize(graph, (size_t) C_0);
    R2.resize(graph, (size_t) C_0);
    timer += time_by(cur);
    auto i_max = (int64) (log2(1.0 * graph.n / eps / eps / k / A.size()) + 1);
    double d0 = log(3.0 * i_max / delta);
    for (int64 i = 1; i <= i_max; i++) {
        bi_seeds.clear();
        cur = clock();
        auto upperC = (is_tightened ? 1.0 : 2.0) * RR_OPIM_Selection(graph, A, k, bi_seeds, R1, is_tightened);
        auto lowerC = (double) R2.self_inf_cal(graph, bi_seeds);
        //printf("theta0: %zu upper: %.3f lower: %.3f\n", R1.numOfRRsets(), upperC, lowerC);
        timer0 += time_by(cur);
        double lower = sqr(sqrt(lowerC + 2.0 * d0 / 9.0) - sqrt(d0 / 2.0)) - d0 / 18.0;
        double upper = sqr(sqrt(upperC + d0 / 2.0) + sqrt(d0 / 2.0));
        double a0 = lower / upper;
        if (a0 > 0.5 - eps || i == i_max) break;
        cur = clock();
        R1.resize(graph, R1.R.size() * 2ll);
        R2.resize(graph, R2.R.size() * 2ll);
        timer += time_by(cur);
    }
    //printf("number of sets: %ld time = %.3f time0 = %.3f\n", R1.numOfRRsets() * 2, timer, timer0);
}

int64
Greedy_OPIM_Selection(Graph &graph, std::vector<int64> &A, int64 k, std::vector<bi_node> &bi_seeds, RRContainer &RRI,
                      bool is_tightened) {
    ///temporary varible
    CandidateNeigh candidate(graph, A, k);
    std::vector<bool> RISetCovered(RRI.R.size(), false);
    std::vector<bool> nodeRemain(graph.n, false);
    for (int64 i : candidate.N) nodeRemain[i] = true;
    auto *coveredNum_tmp = new int64[graph.n];
    memcpy(coveredNum_tmp, RRI.coveredNum, graph.n * sizeof(int64));
    std::priority_queue<std::pair<int64, int64>> Q;
    for (int64 i : candidate.N) Q.push(std::make_pair(coveredNum_tmp[i], i));

    int64 xx = INT64_MAX;

    int64 current_influence = 0;
    while (!Q.empty()) {
        int64 value = Q.top().first;
        int64 maxInd = Q.top().second;
        Q.pop();

        int64 u0 = candidate.source_participant(maxInd);
        if (u0 == -1) {
            nodeRemain[maxInd] = false;
            continue;
        }

        if (value > coveredNum_tmp[maxInd]) {
            Q.push(std::make_pair(coveredNum_tmp[maxInd], maxInd));
            continue;
        }
        if (is_tightened) {
            int64 tight_mg = 0;
            std::vector<int64> tight_list;
            for (auto item : candidate.N) {
                if (coveredNum_tmp[item] > 0)
                    tight_list.emplace_back(coveredNum_tmp[item]);
            }
            int64 k_max = std::min(tight_list.size(), k * A.size());
            std::nth_element(tight_list.begin(), tight_list.begin() + k_max - 1, tight_list.end(), std::greater<>());
            //std::sort(tight_list.begin(), tight_list.end(), std::greater<>());
            for (int j = 0; j < k_max; j++) {
                tight_mg += tight_list[j];
            }
            xx = std::min(xx, tight_mg + current_influence);
        }
        //select
        bi_seeds.emplace_back(maxInd, u0);
        current_influence += coveredNum_tmp[maxInd];
        candidate.choose(u0);
        nodeRemain[maxInd] = false;
        for (int64 RIIndex : RRI.covered[maxInd]) {
            if (RISetCovered[RIIndex]) continue;
            for (int64 u : RRI.R[RIIndex]) {
                coveredNum_tmp[u]--;
            }
            RISetCovered[RIIndex] = true;
        }
    }
    delete[] coveredNum_tmp;
    return is_tightened ? xx : current_influence;
}

void
Greedy_OPIM_Main(Graph &graph, std::vector<int64> &A, int64 k, double eps, double delta, std::vector<bi_node> &bi_seeds,
                 bool is_tightened = false) {
    double timer = 0, timer0 = 0, cur;
    double sum_log = 0;
    for (auto ap:A) sum_log += logcnk(graph.deg_out[ap], k);
    double C_0 = 2.0 * sqr(
            0.5 * sqrt(log(6.0 / delta)) + sqrt(0.5 * (sum_log + log(6.0 / delta))));
    RRContainer R1(graph, A, true), R2(graph, A, true);
    cur = clock();
    R1.resize(graph, (size_t) C_0);
    R2.resize(graph, (size_t) C_0);
    timer += time_by(cur);
    auto i_max = (int64) (log2(1.0 * graph.n / eps / eps / k / A.size()) + 1);
    double d0 = log(3.0 * i_max / delta);
    for (int64 i = 1; i <= i_max; i++) {
        bi_seeds.clear();
        cur = clock();
        auto upperC = (is_tightened ? 1.0 : 2.0) * Greedy_OPIM_Selection(graph, A, k, bi_seeds, R1, is_tightened);
        auto lowerC = (double) R2.self_inf_cal(graph, bi_seeds);
        timer0 += time_by(cur);
        double lower = sqr(sqrt(lowerC + 2.0 * d0 / 9.0) - sqrt(d0 / 2.0)) - d0 / 18.0;
        double upper = sqr(sqrt(upperC + d0 / 2.0) + sqrt(d0 / 2.0));
        double a0 = lower / upper;
        //printf("upper: %.3f lower: %.3f\n", upperC, lowerC);
        if (a0 > 0.5 - eps || i == i_max) break;
        cur = clock();
        R1.resize(graph, R1.R.size() * 2ll);
        R2.resize(graph, R2.R.size() * 2ll);
        timer += time_by(cur);
    }

    //printf("number of sets: %ld time = %.3f time0 = %.3f\n", R1.numOfRRsets() * 2, timer, timer0);
}

/*!
 * @brief Encapsulated method using local OPIM
 * @param graph : the graph
 * @param k : the number in the problem definition
 * @param A : the active participant set A
 * @param seeds : returns the seed set S (each element is a pair <node, AP>)
 */
double method_local_OPIM(Graph &graph, int64 k, std::vector<int64> &A, std::vector<bi_node> &seeds) {
    double cur = clock(), eps = 0.1, delta = 1.0 / graph.n, approx = 1.0 - 1.0 / exp(1);
    std::set<int64> seeds_unique;
    for(auto ap : A) {
        std::vector<int64> ap_vec = {ap};
        std::vector<bi_node> bi_seeds;
        double C_0 = 2.0 * sqr(
                approx * sqrt(log(6.0 / delta)) + sqrt(approx * (logcnk(graph.deg_out[ap], k) + log(6.0 / delta))));
        RRContainer R1(graph, A, true), R2(graph, A, true);
        R1.resize(graph, (size_t) C_0);
        R2.resize(graph, (size_t) C_0);
        auto i_max = (int64) (log2(1.0 * graph.n / eps / eps / k) + 1);
        double d0 = log(3.0 * i_max / delta);
        for (int64 i = 1; i <= i_max; i++) {
            bi_seeds.clear();
            auto upperC = Greedy_OPIM_Selection(graph, ap_vec, k, bi_seeds, R1, true);
            auto lowerC = (double) R2.self_inf_cal(graph, bi_seeds);
            double lower = sqr(sqrt(lowerC + 2.0 * d0 / 9.0) - sqrt(d0 / 2.0)) - d0 / 18.0;
            double upper = sqr(sqrt(upperC + d0 / 2.0) + sqrt(d0 / 2.0));
            double a0 = lower / upper;
            //printf("%.3f %ld\n", a0, R1.numOfRRsets());
            if (a0 > approx - eps || i == i_max) break;
            R1.resize(graph, R1.R.size() * 2ll);
            R2.resize(graph, R2.R.size() * 2ll);
        }
        for (auto u : bi_seeds) {
            if (seeds_unique.find(u.first) == seeds_unique.end()) {
                seeds_unique.insert(u.first);
                seeds.emplace_back(u);
            }
        }
    }
    return time_by(cur);
}

#endif //EXP_OPIM_NEW_H