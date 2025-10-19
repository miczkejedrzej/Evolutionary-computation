#include <vector>
#include <unordered_set>
#include <algorithm>
#include <cassert>
#include <climits>
#include <iostream>
#include "../../src/ProblemInstance.h"
#include "../../src/Solver.h"
#include "GreedySolver.h"


// struct to keep n=capacity best positions with cost for insertion of given unvisited city
struct Node_List {
    std::vector<int> position;
    std::vector<int> cityIndex;
    std::vector<int> cost;
    size_t capacity;

    explicit Node_List(size_t n = 2) : capacity(n) {} // always 2 for 2-regret

    void insert(int pos, int city, int c) {
        if (cost.size() < capacity) {
            auto it = std::lower_bound(cost.begin(), cost.end(), c);
            size_t index = it - cost.begin();
            cost.insert(it, c);
            position.insert(position.begin() + index, pos);
            cityIndex.insert(cityIndex.begin() + index, city);
        } 
        else if (!cost.empty() && c < cost.back()) {
            auto it = std::lower_bound(cost.begin(), cost.end(), c);
            size_t index = it - cost.begin();
            cost.insert(it, c);
            position.insert(position.begin() + index, pos);
            cityIndex.insert(cityIndex.begin() + index, city);
            // keep only top 2
            cost.pop_back();
            position.pop_back();
            cityIndex.pop_back();
        }
    }

    void reset() {
        position.clear();
        cityIndex.clear();
        cost.clear();
    }

    std::vector<int> get_cost() const { return cost; }
    int get_best_pos() const { return position.empty() ? -1 : position.front(); }
};


GreedySolver::GreedySolver(const ProblemInstance& prob, int start_idx, GreedyMode mode, Heuristic heuristic,
float weightObjective)
    : Solver(prob), mode_(mode), heuristic(heuristic), starting_index(start_idx) {}

int GreedySolver::getStartingIndex() { return starting_index; }

void GreedySolver::AssertHamiltonian(const std::vector<int>& visited, int citiesNumber) {
    std::unordered_set<int> uniqueCities(visited.begin(), visited.end());
    assert((int)visited.size() == citiesNumber);
    assert((int)uniqueCities.size() == citiesNumber);
}

void GreedySolver::setWeight(float weight){
    GreedySolver::weight = weight;
}

float GreedySolver::getWeight(){
    return weight;
}
int GreedySolver::ObjectiveEvaluation(const Node_List& node_list, Heuristic heuristic) {
    const auto& costs = node_list.get_cost();
    if (costs.empty()) return INT_MAX;

    int regret2 = (costs[1]-costs[0]);

    if (heuristic == Heuristic::Regret)
        return regret2; // just 2-regret
    else // Hybrid
    {
        float weightObjective = getWeight();
        return static_cast<int>((1.0 - weightObjective) * (regret2) + weightObjective * (-costs[0]));
    }
}




std::vector<int> GreedySolver::solveElementNearestNeighbour(
    int targetSize, const std::vector<int>& visited, const std::vector<int>& unvisited, Heuristic heuristic)
{
    Node_List nodeList(2);
    int bestCityIdx = -1;
    int bestInsertPos = -1;
    int bestScore = INT_MIN; // we maximize

    for (size_t cityIdx = 0; cityIdx < unvisited.size(); ++cityIdx) {
        int city = unvisited[cityIdx];
        // to store best positions to insert  for given unvisited city
        nodeList.reset();

        // compute delta costs for all positions
        for (size_t pos = 0; pos <= visited.size(); ++pos) {
            int deltaCost;
            if (pos == 0)
                deltaCost = problem.GetDistance( visited.front(),city) + problem.GetCost(city);
            else if (pos == visited.size())
                deltaCost = problem.GetDistance(visited.back(), city) +problem.GetCost(city) ;
            else {
                int oldEdge = problem.GetDistance(visited[pos - 1], visited[pos]) +
                problem.GetCost(visited[pos]) ;
                int newEdges = problem.GetDistance(visited[pos - 1], city) + 
                                problem.GetCost(city) +
                               problem.GetDistance(city, visited[pos]) + problem.GetCost(visited[pos]);
                deltaCost = newEdges - oldEdge;
            }
            nodeList.insert(pos, city, deltaCost);
        }

        int score = ObjectiveEvaluation(nodeList, heuristic);
        if (score > bestScore) {
            bestScore = score;
            bestCityIdx = cityIdx;
            bestInsertPos = nodeList.get_best_pos();
        }
    }

    return {bestInsertPos, unvisited[bestCityIdx], bestScore};
}


std::vector<int> GreedySolver::solveElementGreedyCycle(
    int targetSize, const std::vector<int>& visited, const std::vector<int>& unvisited, Heuristic heuristic)
{
    Node_List nodeList(2);
    int bestCityIdx = -1;
    int bestInsertPos = -1;
    int bestScore = INT_MIN; // we maximize cost

    for (size_t cityIdx = 0; cityIdx < unvisited.size(); ++cityIdx) {
        int city = unvisited[cityIdx];
        nodeList.reset();
        size_t nVisited = visited.size();

        for (size_t pos = 0; pos <= nVisited; ++pos) {
            int deltaCost;
            if (pos == 0)
                deltaCost = problem.GetDistance(visited.front(), city) +
                            problem.GetCost(city) +
                            problem.GetDistance(city, visited.back()) -
                            problem.GetDistance(visited.front(), visited.back());
            else if (pos == nVisited)
                deltaCost = problem.GetDistance(visited.back(), city) +
                            problem.GetCost(city) +
                            problem.GetDistance(city, visited.front()) -
                            problem.GetDistance(visited.back(), visited.front());
            else {
                int oldEdge = problem.GetDistance(visited[pos - 1], visited[pos])+
                problem.GetCost(visited[pos]);
                int newEdges = problem.GetDistance(visited[pos - 1], city) +
                                problem.GetCost(city) +
                               problem.GetDistance(city, visited[pos]) +
                               problem.GetCost(visited[pos]);
                deltaCost = newEdges - oldEdge;
            }
            nodeList.insert(pos, city, deltaCost);
        }

        int score = ObjectiveEvaluation(nodeList, heuristic);
        // we maximize regret or weighted regret,-cost 
        if (score > bestScore) {
            bestScore = score;
            bestCityIdx = cityIdx;
            bestInsertPos = nodeList.get_best_pos();
        }
    }

    return {bestInsertPos, unvisited[bestCityIdx], bestScore};
}


std::vector<int> GreedySolver::solve() {
    int targetSize = problem.GetNumberCitiesInCycle();
    std::vector<int> visited;
    visited.reserve(targetSize);
    std::vector<int> unvisited = problem.GiveIndices();

    int start_idx = getStartingIndex();
    visited.push_back(unvisited[start_idx]);
    unvisited.erase(unvisited.begin() + start_idx);

    while (visited.size() < targetSize) {
        std::vector<int> step;

        if (mode_ == GreedyMode::NearestNeighbour){
            step = solveElementNearestNeighbour(targetSize, visited, unvisited, heuristic);
            }
        else{
            step = solveElementGreedyCycle(targetSize, visited, unvisited, heuristic);}

        int insertPos = step[0];
        int city = step[1];

        visited.insert(visited.begin() + insertPos, city);
        unvisited.erase(std::remove(unvisited.begin(), unvisited.end(), city), unvisited.end());
    }

    AssertHamiltonian(visited, targetSize);
    return visited;
}
