#include <vector>
#include <iostream>
#include <random>
#include <list>
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

// Something does not work here
std::vector<int> GreedySolver::solveNearestNeighbour() {
    std::cout << "Solving with Nearest Neighbour...\n";

    std::vector<int> visited;
    visited.reserve(problem.GetNumberCitiesInCycle());
    
    // Choose a random starting point
    std::mt19937 g(156064);
    std::vector<int> unvisited = problem.GiveIndices();

    int currIdx = g() % unvisited.size();
    visited.push_back(unvisited[currIdx]);
    unvisited.erase(unvisited.begin() + currIdx);

    // Take its distance to all other points as a distance vector
    std::vector<int64_t> distances(unvisited.size());
    int nearestIdx = -1;
    int64_t nearestDist = INT64_MAX;
    for (int i = 0; i < unvisited.size(); ++i) {
        distances[i] = problem.GetCostAndDistance(currIdx, unvisited[i]);
        if (distances[i] < nearestDist) {
            nearestDist = distances[i];
            nearestIdx = i;
        }
    }

    currIdx = unvisited[nearestIdx];
    visited.push_back(unvisited[nearestIdx]);
    unvisited.erase(unvisited.begin() + nearestIdx);
    distances.erase(distances.begin() + nearestIdx);

    for (int i = 0; i < unvisited.size(); ++i) {
        int64_t currDist = problem.GetCostAndDistance(currIdx, unvisited[i]);
        if (currDist < distances[i]) {
            distances[i] = currDist;
        }
    }

    // Add the nearest point to the visited list by iterating whrough all the list and checking differences between distances
    for (int i = 2; i < problem.GetNumberCitiesInCycle(); ++i) {
        int nearestIdx = -1;
        int64_t nearestDist = INT64_MAX;
        for (int j = 0; j < unvisited.size(); ++j) {
            if (distances[j] < nearestDist) {
                nearestDist = distances[j];
                nearestIdx = j;
            }
        }

        int shortestIdx = -1;
        int64_t smallestDiff = INT64_MAX;
        for (int k = 0; k < visited.size(); ++k) {
            int nextIdx = (k + 1) % visited.size();
            int64_t diff = problem.GetCostAndDistance(visited[k], unvisited[nearestIdx]) +
                           problem.GetCostAndDistance(unvisited[nearestIdx], visited[nextIdx]) -
                           problem.GetCostAndDistance(visited[k], visited[nextIdx]);
            
            if (diff < smallestDiff) {
                smallestDiff = diff;
                shortestIdx = nextIdx;
            }
        }

        visited.insert(visited.begin() + shortestIdx, unvisited[nearestIdx]);

        unvisited.erase(unvisited.begin() + nearestIdx);
        distances.erase(distances.begin() + nearestIdx);
        for (int i = 0; i < unvisited.size(); ++i) {
            int64_t newDistance = problem.GetCostAndDistance(visited[shortestIdx], unvisited[i]);
            if (newDistance < distances[i]) {
                distances[i] = newDistance;
            }            
        }
    }

    return visited;
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
