//
// Created by asuka on 2022/7/15.
//

#ifndef EXP_MODELS_H
#define EXP_MODELS_H
#include <vector>
#include <random>
#include <iostream>
#include <algorithm>

#define MAX_NODE_SIZE 5000000
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

typedef int64_t node;
typedef int32_t int32;
typedef int64_t int64;

const std::string solver_name[] = {"ENUMERATION", "DEGREE", "PAGERANK", "CELF", "DEGREE_ADVANCED", "PAGERANK_ADVANCED",
                        "CELF_ADVANCED", "IMM_NORMAL", "IMM_ADVANCED", "OPIM_NORMAL", "OPIM_ADVANCED"};

std::random_device rd__;
std::minstd_rand random_engine(rd__());
std::uniform_real_distribution<double> distrib(0.0, 1.0);

std::ofstream out;
int8_t verbose_flag, local_mg;
int64_t MC_iteration_rounds = 10000;

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
void print_set(std::vector<node> &S, const std::string &Prefix = "") {
    std::vector<node> S_ordered = S;
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
void print_small_set(std::vector<node> &S, const std::string &Prefix = "") {
    if(S.size() <= 100) print_set(S, Prefix);
    else {
        std::cout << Prefix;
        std::cout << "{...}[" << S.size() << "]";
    }
}

/*!
 * @brief Print all elements to the file in a vector.
 * @param S : the set
 */
void print_set_f(std::vector<node> &S, const std::string &Prefix = "") {
    if(S.size() > 100) return;
    std::vector<node> S_ordered = S;
    std::sort(S_ordered.begin(), S_ordered.end());
    out << Prefix;
    out << "{";
    for (int64 i = 0; i < S_ordered.size(); i++) {
        out << S_ordered[i];
        if (i != S_ordered.size() - 1) out << ",";
    }
    out << "}[" << S_ordered.size() << "]";
}

#endif //EXP_MODELS_H
