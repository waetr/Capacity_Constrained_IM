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

static std::string graphFilePath;

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
            .add_option<int64>("-r", "--rounds", "number of MC simulation iterations per time, default is 10000", 10000)
            .parse(argc, argv);
    graphFilePath = "../data/" + args.get_argument_string("input");
    if (args.has_option("--verbose")) {
        verbose_flag = 1;
        std::cout << "verbose flag set to 1\n";
    }
    MC_iteration_rounds = args.get_option_int64("--rounds");
    std::cout << "MC_iteration_rounds set to " << MC_iteration_rounds << std::endl;
}

#endif //EXP_TOP_H
