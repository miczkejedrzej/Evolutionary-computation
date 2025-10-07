#include <vector>
#include <iostream>
#include <random>
#include "GreedySolver.h"

// enum class GreedyMode {
//     NearestNeighbour,
//     NearestNeighbourEnd,
//     GreedyCycle
// };

GreedySolver::GreedySolver(const ProblemInstance& prob, GreedyMode mode)
    : Solver(prob), mode_(mode) {}

std::vector<int> GreedySolver::solve() {
    switch (mode_) {
        case GreedyMode::NearestNeighbour:
            return solveNearestNeighbour();
        case GreedyMode::NearestNeighbourEnd:
            return solveNearestNeighbourEnd();
        case GreedyMode::GreedyCycle:
            return solveGreedyCycle();
        default:
            return {};
    }
}

std::vector<int> GreedySolver::solveNearestNeighbour() {
    std::vector<int> result;
    std::cout << "Solving with Nearest Neighbour...\n";
    // TODO: implement

    // More complicated - keep track of distances and update them based on last node added

    return result;
}

std::vector<int> GreedySolver::solveNearestNeighbourEnd() {
    std::cout << "Solving with Nearest Neighbour from End...\n";
    
    std::mt19937 g(156064);
    std::vector<int> unvisited = problem.GiveIndices();
    std::vector<int> visited;
    visited.reserve(problem.GetNumberCitiesInCycle() * sizeof(int));

    int currIdx = g() % unvisited.size();
    visited.push_back(unvisited[currIdx]);
    unvisited.erase(unvisited.begin() + currIdx);

    for (int i = 1; i < problem.GetNumberCitiesInCycle(); ++i) {
        int nearestIdx = -1;
        int64_t nearestDist = INT64_MAX;
        for (int j = 0; j < unvisited.size(); ++j) {
            int64_t dist = problem.GetCostAndDistance(visited[i - 1], unvisited[j]);
            if (dist < nearestDist) {
                nearestDist = dist;
                nearestIdx = j;
            }
        }
        visited.push_back(unvisited[nearestIdx]);
        unvisited.erase(unvisited.begin() + nearestIdx);
    }

    return visited;
}

std::vector<int> GreedySolver::solveGreedyCycle() {
    std::vector<int> result;
    std::cout << "Solving with Greedy Cycle...\n";
    // TODO: implement
    return result;
}
