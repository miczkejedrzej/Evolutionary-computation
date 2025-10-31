#include "LocalSearchSolver.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <valarray>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <random>
#include <cassert>
#include <unordered_set>
#include <chrono>
LocalSearchSolver::LocalSearchSolver(const ProblemInstance& prob,
                                     int randomSeed)
    : Solver(prob),
      rng_(randomSeed)
{}

// ---------------- STRUCT FOR CANDIDATE EDGES ----------------

struct Node_List {
    std::vector<int> startIndex;
    std::vector<int> endIndex;
    std::vector<int> cost;
    size_t capacity;

    explicit Node_List(size_t n = 10) : capacity(n) {}

    void insert(int start, int end, int c) {
        if (cost.size() < capacity) {
            auto it = std::lower_bound(cost.begin(), cost.end(), c);
            size_t index = it - cost.begin();
            cost.insert(it, c);
            startIndex.insert(startIndex.begin() + index, start);
            endIndex.insert(endIndex.begin() + index, end);
        } 
        else if (!cost.empty() && c < cost.back()) {
            auto it = std::lower_bound(cost.begin(), cost.end(), c);
            size_t index = it - cost.begin();
            cost.insert(it, c);
            startIndex.insert(startIndex.begin() + index, start);
            endIndex.insert(endIndex.begin() + index, end);
            
            cost.pop_back();
            startIndex.pop_back();
            endIndex.pop_back();
        }
    }

    void reset() {
        startIndex.clear();
        endIndex.clear();
        cost.clear();
    }

    // std::vector<int> get_cost() const { return cost; }
    // int get_best_pos() const { return startIndex.empty() ? -1 : startIndex.front(); }
    std::vector<int> get_best_indices() const {
        std::vector<int> indices;
        for (size_t i = 0; i < startIndex.size(); ++i) {
            indices.emplace_back(endIndex[i]);
        }
        return indices;
    }
};

void LocalSearchSolver::fillCandidateEdges() {
    candidateEdges.resize(problem.getNumCities());
    Node_List nodeList(10);

    for (int i = 0; i < static_cast<int>(problem.getNumCities()); ++i) {
        nodeList.reset();
        candidateEdges[i].clear();
        for (int j = 0; j < static_cast<int>(problem.getNumCities()); ++j) {
            if (i == j) continue;
            nodeList.insert(i, j, problem.GetCostAndDistance(i, j));
        }

        candidateEdges[i] = nodeList.get_best_indices();
    }
}

// ---------------- DELTA CALCULATIONS ----------------

int LocalSearchSolver::calculateDeltaInter(const std::vector<int>& solution, int selIdx, int unselIdx, const std::vector<int>& unselected) {
    int oldNode = solution[selIdx];
    int newNode = unselected[unselIdx];

    int prev = (selIdx == 0) ? solution.back() : solution[selIdx - 1];
    int next = (selIdx == solution.size() - 1) ? solution.front() : solution[selIdx + 1];

    int newEdge = problem.GetDistance(prev, newNode) + problem.GetDistance(newNode, next) + problem.GetCost(newNode);
    int oldEdge = problem.GetDistance(prev, oldNode) + problem.GetDistance(oldNode, next) + problem.GetCost(oldNode);

    int delta  = newEdge - oldEdge;      
    return delta;
}

int LocalSearchSolver::calculateDeltaIntraTwoEdge(const std::vector<int>& solution, int i, int j) {
    if (j <= i+1 || ((i==0) && (j == solution.size()-1))) return INT_MAX; // skip adjacent edges

    // it is pointles to consider adjacent edges, since no change would be made 
    // similarly to intra node exchange, our case is symmetric so we do not need to consider cases of 
    //  i >  j 

    int next_i = (i==solution.size()-1) ? solution.front() : solution[i+1];
    int next_j = (j==solution.size()-1) ? solution.front() : solution[j+1];
    int oldEdges = problem.GetDistance(solution[i],next_i) +  problem.GetDistance(solution[j],next_j);
    int newEdges = problem.GetDistance(solution[i],solution[j]) + problem.GetDistance(next_i,next_j);

    return newEdges - oldEdges;
}


// ---------------- Helper Function ----------------

std::vector<int> LocalSearchSolver::slicing(const std::vector<int>& arr, int X) {
    // Ensure X is within bounds
    assert(X >= 0 && X < static_cast<int>(arr.size()));

    // Copy first X+1 elements
    std::vector<int> result(arr.begin(), arr.begin() + X + 1);
    return result;
}

void LocalSearchSolver::AssertHamiltonian(const std::vector<int>& visited, int citiesNumber) {
    std::unordered_set<int> uniqueCities(visited.begin(), visited.end());
    assert((int)visited.size() == citiesNumber);
    assert((int)uniqueCities.size() == citiesNumber);
}

// ---------------- MOVE SELECTION ----------------

MoveDelta LocalSearchSolver::findBestMove(const std::vector<int>& solution, const std::vector<int>& unselected) {
    // auto start = std::chrono::high_resolution_clock::now();
    
    MoveDelta bestMove{-1, -1, INT_MAX, MoveType::InterNode};

    // Find all of the candidate edges
    // Use the Node_List from assignment 2 to get them

    std::vector<int> allIndices = problem.GiveIndices();

    for (int startNodeIdx = 0; startNodeIdx < solution.size(); ++startNodeIdx) {
        for (int endNode : candidateEdges[solution[startNodeIdx]]) {
            auto it = std::find(solution.begin(), solution.end(), endNode);
            if (it != solution.end()) {
                int idx_i = startNodeIdx;
                int idx_j = std::distance(solution.begin(), it);
                int delta = calculateDeltaIntraTwoEdge(solution, idx_i, idx_j);
                if (delta < bestMove.delta) {
                    bestMove = {idx_i, idx_j, delta, MoveType::IntraEdge};
                }
            } else {
                int idx_i = startNodeIdx;
                int idx_j = std::distance(unselected.begin(), std::find(unselected.begin(), unselected.end(), endNode));
                int delta = calculateDeltaInter(solution, idx_i, idx_j, unselected);
                if (delta < bestMove.delta) {
                    bestMove = {idx_i, idx_j, delta, MoveType::InterNode};
                }
            }
        }
    }

    // Use only candidate edges as parameters
    // Use different functions depending on if the other node is also in the solution or not

    // auto end = std::chrono::high_resolution_clock::now();
    // int duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    // std::cout << "Single move selection time: " << duration << " ns" << std::endl;

    return bestMove;
}

// ---------------- STARTING SOLUTION ----------------

std::vector<int> LocalSearchSolver::initializeSolution() {
    std::vector<int> sol = problem.GiveIndices();
    std::shuffle(sol.begin(), sol.end(), rng_);
    int numCitiesInCycle = static_cast<int>(problem.GetNumberCitiesInCycle());
    std::vector<int> result = slicing(sol,numCitiesInCycle - 1);
    return result;
}

// ---------------- MAIN SOLVE LOOP ----------------

std::vector<int> LocalSearchSolver::solve() {
    fillCandidateEdges();
    std::vector<int> solution = initializeSolution();
    std::vector<int> unselected = problem.GiveIndices();

    // Remove selected from unselected
    for (int idx : solution)
        unselected.erase(std::remove(unselected.begin(), unselected.end(), idx), unselected.end());
    while (true) {
        MoveDelta move = findBestMove(solution, unselected);

        if (move.delta >= 0) break;

        switch (move.type) {
            case MoveType::InterNode:
                std::swap(solution[move.i1], unselected[move.i2]);
                break;
            case MoveType::IntraEdge:
                if (move.i1 > move.i2) std::swap(move.i1, move.i2);
                    std::reverse(solution.begin() +move.i1 +1, solution.begin()+move.i2+ 1);
                
                break;
        }
    }
    AssertHamiltonian(solution,problem.GetNumberCitiesInCycle());
    return solution;
}
