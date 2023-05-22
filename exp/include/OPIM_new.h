#ifndef EXP_OPIM_NEW_H
#define EXP_OPIM_NEW_H

static int64 *coveredNum_tmp;
static bool *nodeRemain;

/*!
 * @brief given RR sets, using RR-method to select the seed set based on these RR sets.
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
    assert(bi_seeds.empty());
    ///temporary varible
    ///initialization
    std::set<int64> A_reorder(A.begin(), A.end());
    for (int i = 0; i < A.size(); i++) {
        for (auto e : graph.g[A[i]])
            if (A_reorder.find(e.v) == A_reorder.end()) {
                nodeRemain[e.v] = true;
            }
    }
    A_reorder.clear();
    std::vector<bool> RISetCovered(RRI.R.size(), false);
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
                for (auto e : graph.g[A[i]])
                    if (coveredNum_tmp[e.v] > 0) tight_list.emplace_back(coveredNum_tmp[e.v]);
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
    for (int i = 0; i < A.size(); i++) {
        for (auto e : graph.g[A[i]])
            nodeRemain[e.v] = false;
    }
    return is_tightened ? xx : current_influence;
}


int64
MG_OPIM_Selection(Graph &graph, std::vector<int64> &A, int64 k, std::vector<bi_node> &bi_seeds, RRContainer &RRI,
                  bool is_tightened) {
    assert(bi_seeds.empty());
    ///temporary varible
    CandidateNeigh candidate(graph, A, k);
    std::vector<bool> RISetCovered(RRI.R.size(), false);
    for (int64 i : candidate.N) nodeRemain[i] = true;
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
    for (int64 i : candidate.N) nodeRemain[i] = false;
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
double
method_FOPIM(Graph &graph, int64 k, std::vector<int64> &A, std::vector<bi_node> &bi_seeds, double eps,
             std::string type) {
    assert(bi_seeds.empty());
    const double delta = 1.0 / graph.n, xi = 0.05;
    const double approx = 0.5;
    std::vector<bi_node> pre_seeds;
    method_random(graph, k, A, pre_seeds);
    int64 opt_lower_bound = pre_seeds.size();
    RRContainer R1(graph, A, true), R2(graph, A, true);
    coveredNum_tmp = new int64[graph.n];
    nodeRemain = new bool[graph.n];
    auto start_time = std::chrono::high_resolution_clock::now();
    double time1 = 0, time2 = 0, cur;
    double sum_log = 0;
    for (auto ap:A) sum_log += logcnk(graph.deg_out[ap], k);
    double C_max = 2.0 * graph.n * sqr(
            approx * sqrt(log(6.0 / delta)) + sqrt(approx * (sum_log + log(6.0 / delta)))) / eps / eps /
                   opt_lower_bound;
    double C_0 = 2.0 * sqr(approx * sqrt(log(6.0 / delta)) + sqrt(approx * (sum_log + log(6.0 / delta)))) /
                 sqrt(opt_lower_bound);
    cur = clock();
    R1.resize(graph, (size_t) C_0);
    R2.resize(graph, (size_t) C_0);
    time1 += time_by(cur);
    auto i_max = (int64) (log2(C_max / C_0) + 1);
    //printf("lower-bound of optimal: %ld i_max: %ld\n", opt_lower_bound, i_max);
    double d0 = log(3.0 * i_max / delta);
    for (int64 i = 1; i <= i_max; i++) {
        bi_seeds.clear();
        cur = clock();
        double upperC;
        if (type == "RR+") upperC = RR_OPIM_Selection(graph, A, k, bi_seeds, R1, true);
        else if (type == "RR") upperC = RR_OPIM_Selection(graph, A, k, bi_seeds, R1, false) / approx;
        else if (type == "MG") upperC = MG_OPIM_Selection(graph, A, k, bi_seeds, R1, false) / approx;
        else {
            std::cerr << "error: invalid FOPIM type!\n";
            exit(1);
        }
        auto lowerC = (double) R2.self_inf_cal(graph, bi_seeds);
        time2 += time_by(cur);
        double lower = sqr(sqrt(lowerC + 2.0 * d0 / 9.0) - sqrt(d0 / 2.0)) - d0 / 18.0;
        double upper = sqr(sqrt(upperC + d0 / 2.0) + sqrt(d0 / 2.0));
        double a0 = lower / upper;
        //printf("a0:%.3f theta0:%zu lowerC: %.3f upperC: %.3f\n", a0, R1.R.size(), lowerC, upperC);
        if (a0 > 0.5 - eps || i == i_max) break;
        cur = clock();
        R1.resize(graph, R1.R.size() * 2ll);
        R2.resize(graph, R2.R.size() * 2ll);
        time1 += time_by(cur);
    }
    printf("time1: %.3f time2: %.3f size: %zu, size1: %zu\n", time1, time2, R1.numOfRRsets(), R1.sizeOfRRsets());
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    delete[] coveredNum_tmp;
    delete[] nodeRemain;
    return elapsed.count();
}


/*!
 * @brief Encapsulated method using local OPIM
 * @param graph : the graph
 * @param k : the number in the problem definition
 * @param A : the active participant set A
 * @param seeds : returns the seed set S (each element is a pair <node, AP>)
 */
double method_local_OPIM(Graph &graph, int64 k, std::vector<int64> &A, double eps, std::vector<bi_node> &seeds) {
    assert(seeds.empty());
    coveredNum_tmp = new int64[graph.n];
    nodeRemain = new bool[graph.n];
    auto start_time = std::chrono::high_resolution_clock::now();
    double delta = 1.0 / graph.n, approx = 1.0 - 1.0 / exp(1);
    std::set<int64> seeds_unique;
    for (auto ap : A) {
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
            auto upperC = MG_OPIM_Selection(graph, ap_vec, k, bi_seeds, R1, true);
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
            if (std::find(A.begin(), A.end(), u.first) == A.end() &&
                seeds_unique.find(u.first) == seeds_unique.end()) {
                seeds_unique.insert(u.first);
                seeds.emplace_back(u);
            }
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    delete[] coveredNum_tmp;
    delete[] nodeRemain;
    return elapsed.count();
}

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
    for (int64 i : candidate) nodeRemain[i] = true;
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
    for (int64 i : candidate) nodeRemain[i] = false;
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
                           (iota * log(graph.n) + logcnk(candidate.size(), k) + log(log2(graph.n))) /
                           sqr(epsilon_prime) *
                           pow(2.0, i));
        RRI.resize(graph, ci);

        double ept = IMMNodeSelection(graph, candidate, k, S_tmp, RRI);
        if (ept > 1.0 / pow(2.0, i)) {
            LB = ept * graph.n / (1.0 + epsilon_prime);
            break;
        }
    }
    double e = exp(1);
    double alpha = sqrt(iota * log(graph.n) + log(2));
    double beta = sqrt((1.0 - 1.0 / e) * (logcnk(candidate.size(), k) + iota * log(graph.n) + log(2)));
    auto C = (int64) (2.0 * graph.n * sqr((1.0 - 1.0 / e) * alpha + beta) / LB / sqr(eps));
    RRI.resize(graph, C);
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
void IMM_full(Graph &G, std::vector<int64> &candidate, int64 k, double eps, double iota, std::vector<int64> &S,
              RRContainer &RRI) {
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
double method_local_IMM(Graph &graph, int64 k, std::vector<int64> &A, double eps, std::vector<bi_node> &seeds) {
    assert(seeds.empty());
    coveredNum_tmp = new int64[graph.n];
    nodeRemain = new bool[graph.n];
    auto start_time = std::chrono::high_resolution_clock::now();
    std::set<int64> seeds_reorder;
    for (int64 u : A) {
        std::vector<int64> neighbours, one_seed;
        for (auto &edge : graph.g[u]) {
            if (std::find(A.begin(), A.end(), edge.v) == A.end())
                neighbours.emplace_back(edge.v);
        }
        std::vector<int64> ap = {u};
        RRContainer RRI(graph, A, true);
        IMM_full(graph, neighbours, k, eps, 1, one_seed, RRI);
        for (int64 w : one_seed) {
            assert(std::find(A.begin(), A.end(), w) == A.end());
            if (seeds_reorder.find(w) == seeds_reorder.end()) {
                seeds_reorder.insert(w);
                seeds.emplace_back(w, u);
            }
        }
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    delete[] coveredNum_tmp;
    delete[] nodeRemain;
    return elapsed.count();
}

#endif //EXP_OPIM_NEW_H
