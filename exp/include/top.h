#ifndef EXP_TOP_H
#define EXP_TOP_H

#include "IMs.h"
#include "argparse.h"
#include "OPIM_new.h"

static std::string graphFilePath;
static std::string ApFilePath;

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
            .add_option<int64>("-r", "--rounds", "number of MC simulation iterations per time, default is 10000", 10000)
            .parse(argc, argv);
    graphFilePath = "../data/" + args.get_argument_string("input") + ".txt";
    ApFilePath = "../data/" + args.get_argument_string("input") + ".ap";
    MC_iteration_rounds = args.get_option_int64("--rounds");
    //std::cout << "MC_iteration_rounds set to " << MC_iteration_rounds << std::endl;
}

std::vector<std::vector<int64>> AP_from_file(const std::string &filename) {
    std::vector<std::vector<int64>> res;
    std::ifstream inFile(filename, std::ios::in);
    if (!inFile.is_open()) {
        std::cerr << "(get error) AP file not found: " << filename << std::endl;
        std::exit(-1);
    }
    std::string line, word;
    std::istringstream sin;
    std::vector<int64> one_AP;
    while (getline(inFile, line)) {
        if (line.empty()) continue;
        sin.clear();
        sin.str(line);
        one_AP.clear();
        while (std::getline(sin, word, ','))
            one_AP.emplace_back(std::stoi(word));
        res.emplace_back(one_AP);
    }
    inFile.close();
    return res;
}

#endif //EXP_TOP_H
