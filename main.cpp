#include "ProblemInstance.h"
// #include "GreedySolver.h"
#include "RandomSolver.h"
#include <iostream>

int main() {
    ///the pipeline with random
    ProblemInstance prob("./TSPA.csv");
    //inside the solver result printed
    RandomSolver solver(prob);
    std::vector<int> result = solver.solve();
    return 0;
}