
#include "RandomSolver.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include <random>

std::vector<int> RandomSolver::slicing(const std::vector<int> arr, int X) {

    auto end = arr.begin() + X;
    auto start = arr.begin();
    std::vector<int> result(X + 1);
    std::copy(start, end, result.begin());
    return result;
}

std::vector<int> RandomSolver::solve() {
    //random generator function 
    std::mt19937 g(std::random_device{}());
    int numCitiesInCycle = static_cast<int>(problem.GetNumberCitiesInCycle());
    std::vector<int> indices = problem.GiveIndices();
    //randomly shuffling the array
    std::shuffle(indices.begin(), indices.end(), g);
    std::vector<int> result = slicing(indices,numCitiesInCycle - 1);
    int64_t total_cost = problem.FullDistanceAndCost(result);
    // std::cout << "Random cost is " << total_cost << std::endl;
    return result;
}


RandomSolver::RandomSolver(const ProblemInstance& prob)
    : Solver(prob) {}  