#include "top.h"
using namespace std;

int main(int argc, char const *argv[]) {
    init_commandLine(argc, argv);
    vector<node> A_batch = {5};
    vector<int32> k_batch = {5};
    vector<IM_solver> solver_batch = {CELF_NORMAL, CELF_ADVANCED, CELF_THRESHOLD1, CELF_THRESHOLD2};
    Run_simulation(A_batch, k_batch, solver_batch, IC, 3);
    return 0;
}