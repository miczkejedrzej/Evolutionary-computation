#include "ProblemInstance.h"
#include "GreedySolver.h"
#include "RandomSolver.h"
#include <iostream>

int main() {
    ///the pipeline with random
    ProblemInstance prob("./TSPA.csv");
    //inside the solver result printed
    RandomSolver solver(prob);
    std::vector<int> result = solver.solve();

    GreedySolver greedySolver(prob, GreedyMode::NearestNeighbourEnd);
    std::vector<int> greedyResult = greedySolver.solve();
    std::cout << "Greedy Result (NearestNeighbourEnd): " << prob.FullDistanceAndCost(greedyResult) << std::endl;

    GreedySolver greedySolver2(prob, GreedyMode::NearestNeighbour);
    greedyResult = greedySolver2.solve();
    std::cout << "Greedy Result: (NearestNeighbour): " << prob.FullDistanceAndCost(greedyResult) << std::endl;

    return 0;
}