//
// Created by asuka on 2022/7/15.
//

#ifndef EXP_MODELS_H
#define EXP_MODELS_H

#include <vector>
#include <random>
#include <iostream>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <cstring>
#include <sstream>
#include <set>
#include <queue>
#include <stack>
#include <bitset>

#define MAX_NODE_SIZE 50000000
#define graph_type int8_t
#define DIRECTED_G 0
#define UNDIRECTED_G 1

#define model_type int8_t
#define NONE 0
#define IC 1
#define LT 2
#define IC_M 3

#define IM_solver int8_t
#define ENUMERATION 0
#define DEGREE 1
#define PAGERANK 2
#define CELF_NORMAL 3
#define DEGREE_ADVANCED 4
#define PAGERANK_ADVANCED 5
#define CELF_ADVANCED 6
#define IMM_NORMAL 7
#define IMM_ADVANCED 8
#define OPIM_NORMAL 9
#define OPIM_ADVANCED 10
#define CELF_THRESHOLD 11
#define CELF_THRESHOLD1 12
#define CELF_THRESHOLD2 13

typedef int64_t int64;
typedef std::pair<int64, int64 > bi_node; //first: seed second: source AP

const std::string solver_name[] = {"ENUMERATION", "DEGREE", "PAGERANK", "CELF", "DEGREE_ADVANCED", "PAGERANK_ADVANCED",
                                   "CELF_ADVANCED", "IMM_NORMAL", "IMM_ADVANCED", "OPIM_NORMAL", "OPIM_ADVANCED",
                                   "CELF_THRESHOLD", "CELF_THRESHOLD1", "CELF_THRESHOLD2"};

std::random_device rd__;
std::minstd_rand random_engine(rd__());
std::mt19937 mt19937engine(rd__());
std::uniform_real_distribution<double> distrib(0.0, 1.0);

std::ofstream stdFileOut;
int8_t verbose_flag;
int64_t MC_iteration_rounds = 10000;

double logcnk(int n, int k) {
    if(k >= n) return 0;
    k = k < n - k ? k : n - k;
    double res = 0;
    for (auto i = 1; i <= k; i++) res += log(double(n - k + i) / i);
    return res;
}

template<class T>
T sqr(T x) {
    return x * x;
}

/*!
 * @brief Random number generator that generates real number between [0,1)
 * @return A random real number between [0,1)
 */
inline double random_real() {
    return distrib(random_engine);
}

/*!
 * @brief calculate the interval from start time.
 * @param start : start timestamp
 * @return the length of the interval
 */
double time_by(double start) {
    return (clock() - start) / CLOCKS_PER_SEC;
}

/*!
 * @brief Print all elements in a vector.
 * @param S : the set
 */
void print_set(std::vector<int64> &S, const std::string &Prefix = "") {
    std::vector<int64> S_ordered = S;
    std::sort(S_ordered.begin(), S_ordered.end());
    std::cout << Prefix;
    std::cout << "{";
    for (int64 i = 0; i < S_ordered.size(); i++) {
        std::cout << S_ordered[i];
        if (i != S_ordered.size() - 1) std::cout << ",";
    }
    std::cout << "}[" << S_ordered.size() << "]";
}

/*!
 * @brief Print the set if the size of the set is not exceeding 100.
 * @param S : the set
 */
void print_small_set(std::vector<int64> &S, const std::string &Prefix = "") {
    if (S.size() <= 100) print_set(S, Prefix);
    else {
        std::cout << Prefix;
        std::cout << "{...}[" << S.size() << "]";
    }
}

/*!
 * @brief Print all elements to the file in a vector.
 * @param S : the set
 */
void print_set_f(std::vector<int64> &S, const std::string &Prefix = "") {
    if (S.size() > 100) return;
    std::vector<int64> S_ordered = S;
    std::sort(S_ordered.begin(), S_ordered.end());
    stdFileOut << Prefix;
    stdFileOut << "{";
    for (int64 i = 0; i < S_ordered.size(); i++) {
        stdFileOut << S_ordered[i];
        if (i != S_ordered.size() - 1) stdFileOut << ",";
    }
    stdFileOut << "}[" << S_ordered.size() << "]";
}

inline int64 fast_read(std::ifstream &inFile) {
    int64 x = 0;
    auto ch = inFile.get();
    if (ch == EOF) return -1;
    while (ch < '0' || ch > '9') {
        if (ch == EOF) return -1;
        ch = inFile.get();
    }
    while (ch >= '0' && ch <= '9') {
        x = x * 10 + ch - '0';
        ch = inFile.get();
    }
    return x;
}

void ASSERT_SEED(std::vector<bi_node> &seeds) {
    std::set<int64> unique_;
    for(auto &e:seeds) {
        if(unique_.find(e.first) != unique_.end()) {
            printf("FATAL: SEED GET INTO RUIN!\n");
            return;
        }
        unique_.insert(e.first);
    }
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

#endif //EXP_MODELS_H
