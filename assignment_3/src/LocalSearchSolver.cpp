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
LocalSearchSolver::LocalSearchSolver(const ProblemInstance& prob,
                                     LocalSearchType searchType,
                                     MoveType intraType,
                                     StartSolutionType startType,
                                     int randomSeed,
                                    int startingIndex)
    : Solver(prob),
      searchType_(searchType),
      intraType_(intraType),
      startType_(startType),
      rng_(randomSeed),
      startingIndex(startingIndex)
{}

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

int LocalSearchSolver::calculateDeltaIntraTwoNode(const std::vector<int>& solution, int i, int j) {
    // this condition made since teh case is symmetric , ie exchanging positiong between nodes 
    // i and j is tantamount to exchanging to position between j and i  
    if (j <= i) return INT_MAX;

    int n = static_cast<int>(solution.size());
    if (i == 0 && j == n-1) {
    // Only three edges: (last, first), (first, next_i), (prev_j, last)
    int prev_j = solution[j-1];
    int first = solution[i];
    int last = solution[j];
    int next_i = solution[i+1];

    int oldEdges = problem.GetDistance(first, next_i) + problem.GetDistance(prev_j, last);
    int newEdges = problem.GetDistance(prev_j, first)  + problem.GetDistance(last, next_i);

    return newEdges - oldEdges;
}

    // ----- Adjacent nodes -----
    if (j == i + 1 || (i==0 && j == n-1)) {
        int prev = (i == 0) ? solution.back() : solution[i - 1];
        int a = solution[i];
        int b = solution[j];
        int next = (j == n - 1) ? solution.front() : solution[j + 1];

        int oldEdges = problem.GetDistance(prev, a)  + problem.GetDistance(b, next);
        int newEdges = problem.GetDistance(prev, b)  + problem.GetDistance(a, next);

        return newEdges - oldEdges;
    }


    int prev_i = (i == 0) ? solution.back() : solution[i - 1];
    int next_i = (i == solution.size() - 1) ? solution.front() : solution[i + 1];

    int prev_j = (j == 0) ? solution.back() : solution[j - 1];
    int next_j = (j == solution.size() - 1) ? solution.front() : solution[j + 1];

    int oldEdges = problem.GetDistance(prev_i, solution[i]) + problem.GetDistance(solution[i], next_i)
                 + problem.GetDistance(prev_j, solution[j]) + problem.GetDistance(solution[j], next_j);

    int newEdges = problem.GetDistance(prev_i, solution[j]) + problem.GetDistance(solution[j], next_i)
                 + problem.GetDistance(prev_j, solution[i]) + problem.GetDistance(solution[i], next_j);

    int delta  = newEdges - oldEdges;   
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
    MoveDelta bestMove{-1, -1, INT_MAX, MoveType::InterRoute};

    // Inter-route moves
    for (size_t selIdx = 0; selIdx < solution.size(); ++selIdx) {
        for (size_t unselIdx = 0; unselIdx < unselected.size(); ++unselIdx) {
            int delta = calculateDeltaInter(solution, selIdx, unselIdx, unselected);
            if (delta < bestMove.delta) {
                bestMove = {static_cast<int>(selIdx), static_cast<int>(unselIdx), delta, MoveType::InterRoute};
            }
        }
    }

    // Intra-route moves
    for (size_t i = 0; i < solution.size(); ++i) {
        for (size_t j = i + 1; j < solution.size(); ++j) {
            int delta = (intraType_ == MoveType::IntraNodeSwap)
                        ? calculateDeltaIntraTwoNode(solution, i, j)
                        : calculateDeltaIntraTwoEdge(solution, i, j);
            if (delta < bestMove.delta) {
                bestMove = {static_cast<int>(i), static_cast<int>(j), delta, intraType_};
            }
        }
    }

    return bestMove;
}

MoveDelta LocalSearchSolver::findRandomGreedyMove(const std::vector<int>& solution, const std::vector<int>& unselected) {
    struct Candidate { bool isInter; int i; int j; };
    std::vector<Candidate> candidates;

    for (size_t i = 0; i < solution.size(); ++i)
        for (size_t j = 0; j < unselected.size(); ++j)
            candidates.push_back({true, static_cast<int>(i), static_cast<int>(j)});

    for (size_t i = 0; i < solution.size(); ++i)
        for (size_t j = i + 1; j < solution.size(); ++j)
            candidates.push_back({false, static_cast<int>(i), static_cast<int>(j)});

    std::shuffle(candidates.begin(), candidates.end(), rng_);

    for (const auto& cand : candidates) {
        int delta = cand.isInter
                    ? calculateDeltaInter(solution, cand.i, cand.j, unselected)
                    : (intraType_ == MoveType::IntraNodeSwap
                        ? calculateDeltaIntraTwoNode(solution, cand.i, cand.j)
                        : calculateDeltaIntraTwoEdge(solution, cand.i, cand.j));

        if (delta < 0)
            return {cand.i, cand.j, delta, cand.isInter ? MoveType::InterRoute : intraType_};
    }

    return {-1, -1, INT_MAX, MoveType::InterRoute};
}

// ---------------- STARTING SOLUTION ----------------

std::vector<int> LocalSearchSolver::initializeSolution() {
    if (startType_ == StartSolutionType::Random) {
        std::vector<int> sol = problem.GiveIndices();
        std::shuffle(sol.begin(), sol.end(), rng_);
        int numCitiesInCycle = static_cast<int>(problem.GetNumberCitiesInCycle());
        std::vector<int> result = slicing(sol,numCitiesInCycle - 1);
        return result;
    }
    
    int starting_index = giveStartingIndex(); // or any index you want to start from
    GreedySolver greedySolver(problem, starting_index, GreedyMode::NearestNeighbour, Heuristic::HybridRegretObjective);
    greedySolver.setWeight(0.5f); // Set weight for hybrid heuristic
        // Solve using greedy heuristic
    std::vector<int> greedy_solution = greedySolver.solve();
    return greedy_solution;
}

int LocalSearchSolver::giveStartingIndex(){
    return startingIndex;
}

// ---------------- MAIN SOLVE LOOP ----------------

std::vector<int> LocalSearchSolver::solve() {
    std::vector<int> solution = initializeSolution();
    std::vector<int> unselected = problem.GiveIndices();

    // Remove selected from unselected
    for (int idx : solution)
        unselected.erase(std::remove(unselected.begin(), unselected.end(), idx), unselected.end());
    while (true) {
        MoveDelta move = (searchType_ == LocalSearchType::Steepest)
                         ? findBestMove(solution, unselected)
                         : findRandomGreedyMove(solution, unselected);

        if (move.delta >= 0) break;

        switch (move.type) {
            case MoveType::InterRoute:
                std::swap(solution[move.i1], unselected[move.i2]);
                break;
            case MoveType::IntraNodeSwap:
                std::swap(solution[move.i1], solution[move.i2]);
                break;
            case MoveType::IntraEdgeSwap:
                // if (move.i1 > move.i2) std::swap(move.i1, move.i2);
                // std::reverse(solution.begin() + move.i1, solution.begin() + move.i2 + 1);
                if (move.i1 > move.i2) std::swap(move.i1, move.i2);
                    std::reverse(solution.begin() +move.i1 +1, solution.begin()+move.i2+ 1);
                
                break;
        }
    }
    AssertHamiltonian(solution,problem.GetNumberCitiesInCycle());
    return solution;
}
