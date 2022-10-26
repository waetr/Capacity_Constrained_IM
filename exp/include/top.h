//
// Created by lenovo on 2022/7/19.
//

#ifndef EXP_TOP_H
#define EXP_TOP_H

#include "IMs.h"
#include "argparse.h"
#include "IMM.h"
#include "OPIM.h"
#include "matroid.h"

std::string graphFilePath;

/*!
 * @brief Initialize the file path of the graph, verbose flag, etc. through command line arguments.
 * @param argc
 * @param argv
 */
void init_commandLine(int argc, char const *argv[]) {
    auto args = util::argparser("The experiment of BIM.");
    args.set_program_name("exp")
            .add_help_option()
            .add_argument<std::string>("input", "initialize graph file")
            .add_option("-v", "--verbose", "output verbose message or not")
            .add_option<std::string>("-l", "--local", "use local value as single spread or not", "")
            .add_option<int64>("-r", "--rounds", "number of MC simulation iterations per time, default is 10000", 10000)
            .parse(argc, argv);
    graphFilePath = "../data/" + args.get_argument_string("input");
    if (args.has_option("--verbose")) {
        verbose_flag = 1;
        std::cout << "verbose flag set to 1\n";
    }
    if (!args.get_option_string("--local").empty()) {
        local_mg = 1;
        std::string MGPath = "../data/" + args.get_option_string("--local");
        std::cout << "local spread flag set to 1, file path = " << MGPath << std::endl;
        std::ifstream inFile(MGPath, std::ios::in);
        if (!inFile.is_open()) {
            std::cerr << "(get error) local file not found: " << args.get_option_string("--local") << std::endl;
            std::exit(-1);
        }
        int64 cnt = 0;
        while (inFile.good()) inFile >> MG0[cnt++];
        inFile.close();
    }
    MC_iteration_rounds = args.get_option_int64("--rounds");
    std::cout << "MC_iteration_rounds set to " << MC_iteration_rounds << std::endl;
}

/*!
 * @brief A BIM Solver collection.
 * @param graph : the graph
 * @param k : the size of seed set of each active participant
 * @param A : the active participant set
 * @param seeds : returns the union of the seed sets as a passing parameter
 * @param solver : identifier of the BIM solver
 * @return time for running BIM Solver
 */
double solvers(Graph &graph, int64 k, std::vector<int64> &A, std::vector<int64> &seeds, IM_solver solver, double &seedAvgDegree) {
    double cur = clock();
    seeds.clear();
    switch (solver) {
        case ENUMERATION:
            enumeration_method(graph, k, A, seeds);
            break;
        case DEGREE:
            degree_method(graph, k, A, seeds);
            break;
        case PAGERANK:
            pgrank_method(graph, k, A, seeds);
            break;
        case CELF_NORMAL:
            CELF_method(graph, k, A, seeds);
            break;
        case DEGREE_ADVANCED:
            advanced_degree_method(graph, k, A, seeds, seedAvgDegree);
            break;
        case PAGERANK_ADVANCED:
            advanced_pgrank_method(graph, k, A, seeds, seedAvgDegree);
            break;
        case CELF_ADVANCED:
            advanced_CELF_method(graph, k, A, seeds, seedAvgDegree);
            break;
        case IMM_NORMAL:
            IMM_method(graph, k, A, seeds);
            break;
        case IMM_ADVANCED:
            advanced_IMM_method(graph, k, A, seeds);
            break;
        case OPIM_NORMAL:
            OPIM_method(graph, k, A, seeds);
            break;
        case OPIM_ADVANCED:
            advanced_OPIM_method(graph, k, A, seeds);
            break;
        case CELF_THRESHOLD:
            //Thresholding_CELF(graph, k, A, seeds);
            break;
        case CELF_THRESHOLD1:
            //Thresholding_CELF1(graph, k, A, seeds, seedAvgDegree);
            break;
        case CELF_THRESHOLD2:
            //Thresholding_CELF2(graph, k, A, seeds, seedAvgDegree);
            break;
        default:
            break;
    }
    print_small_set(seeds, " Seed set using " + solver_name[solver] + ": "), puts("");

    if (verbose_flag) printf(" total time = %.3f\n", time_by(cur));
    return time_by(cur);
}

/*!
 * @brief A procedure that automates batch processing experiments
 * @param A_batch : the collection of ap sizes
 * @param k_batch : the collection of ks
 * @param solver_batch : the collection of BIM solvers
 * @param type : the type of the diffusion model (IC/IC_M)
 * @param rounds : take the number of replicates as the average
 */
void Run_simulation(std::vector<int64> &A_batch, std::vector<int64> &k_batch, std::vector<IM_solver> &solver_batch, model_type type,
                    int64 rounds = 3) {
    std::fstream file_eraser("../output/result.out", std::ios::out);
    file_eraser.close();
    stdFileOut.open("../output/result.out", std::ios::app);
    if(!stdFileOut.is_open()) {
        std::cerr << "(file error)Cannot load result.out! Maybe losing /output?\n";
        exit(-1);
    }
    //load graph from absolute path
    Graph G(graphFilePath, DIRECTED_G);

    //set diffusion model
    G.set_diffusion_model(type, 15);

    //Instantiate the active participant set A and seed set
    std::vector<int64> A, seeds;

    double result[20][500], timer[20][500], seedSize[20][500], seedsAvgDegree[20][500];

    for (int64 A_size : A_batch) {
        double apAvgDegree = 0, overlap_ratio = 0;
        memset(result, 0, sizeof(result));
        memset(timer, 0, sizeof(timer));
        memset(timer, 0, sizeof(seedSize));
        memset(timer, 0, sizeof(seedsAvgDegree));
        for (int64 r_ = 1; r_ <= rounds; r_++) {
            generate_seed(G, A, A_size);
            print_set(A, "active participant: "), puts("");
            print_set_f(A, "active participant: "), stdFileOut << '\n';
            overlap_ratio += estimate_neighbor_overlap(G, A) / rounds;
            for(auto u : A) apAvgDegree += (double) G.deg_out[u] / A.size() / rounds;
            for (int64 k : k_batch) {
                std::cout << "Working on A_size = " << A_size << ", round = " << r_ << ", k = " << k << std::endl;
                stdFileOut << "Working on A_size = " << A_size << ", round = " << r_ << ", k = " << k << std::endl;
                for (IM_solver solver_used : solver_batch) {
                    double seedAvgDegree = 0;
                    timer[solver_used][k] += solvers(G, k, A, seeds, solver_used, seedAvgDegree) / rounds;
                    result[solver_used][k] += FI_simulation(G, seeds) / rounds;
                    seedSize[solver_used][k] += (double) seeds.size() / rounds;
                    seedsAvgDegree[solver_used][k] += (double) seedAvgDegree / rounds;
                }
            }
        }
        std::cout << "*************\nSize " << A_size << " DONE!\n*************\n";

        stdFileOut << "Active participant size = " << A_size << std::endl;
        stdFileOut << "Mean overlap ratio = " << overlap_ratio << std::endl;
        stdFileOut << "Mean ap's degree = " << apAvgDegree << std::endl;
        for (int64 k : k_batch) {
            stdFileOut << "\tk = " << k << std::endl;
            for (IM_solver solver_used : solver_batch) {
                stdFileOut << "\t\tseed quality of " << solver_name[solver_used] << " = " << result[solver_used][k];
                stdFileOut << " time = " << timer[solver_used][k];
                stdFileOut << " size = " << seedSize[solver_used][k] << std::endl;
                stdFileOut << "\t\t\tseed's avg degree = " << seedsAvgDegree[solver_used][k] << std::endl;
            }
        }
    }
    stdFileOut.close();
}

#endif //EXP_TOP_H
