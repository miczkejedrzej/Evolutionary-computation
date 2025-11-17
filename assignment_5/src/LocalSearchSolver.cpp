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
#include <unordered_map>
LocalSearchSolver::LocalSearchSolver(const ProblemInstance& prob,
                                     int randomSeed)
    : Solver(prob),
      rng_(randomSeed)
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
    // ---- Choose steepest descent move ----
    MoveDelta best{-1, -1, INT_MAX, MoveType::InterNode};

    std::unordered_map<int, int> nextMoves;
    nextMoves.reserve(solution.size());
    nextMoves[solution.back()] = solution.front();
    for (int i = 0; i < solution.size() - 1; i++) {
        nextMoves[solution[i]] = solution[i + 1];
    }

    for (auto &m : moveMemory_) {
        // Check if the intra-edge swap move is applicable
        if (m.first.type == MoveType::IntraEdge) {
            if (m.first.i1next != nextMoves[m.first.i1] || m.first.i2next != nextMoves[m.first.i2])
                continue;
        }

        if (m.second < best.delta) {
            best.i1 = m.first.i1;
            best.i2 = m.first.i2;
            best.delta = m.second;
            best.type = m.first.type;
        }
    }

    if (best.delta < 0) {
        auto it1 = std::find(solution.begin(), solution.end(), best.i1);
        assert(it1 != solution.end());
        best.i1 = std::distance(solution.begin(), it1);

        if (best.type == MoveType::InterNode) {
            auto it2 = std::find(unselected.begin(), unselected.end(), best.i2);    
            assert(it2 != unselected.end());
            best.i2 = std::distance(unselected.begin(), it2);
        }
        else {
            auto it2 = std::find(solution.begin(), solution.end(), best.i2);    
            assert(it2 != solution.end());
            best.i2 = std::distance(solution.begin(), it2);
        }
    }

    return best;
}

// ---------------- STARTING SOLUTION ----------------

bool adjacent(int a, int b, int n) {
    return ( (a+1)%n == b ) || ( (b+1)%n == a );
}

std::vector<int> LocalSearchSolver::initializeSolution() {
    std::vector<int> sol = problem.GiveIndices();
    std::shuffle(sol.begin(), sol.end(), rng_);
    int numCitiesInCycle = static_cast<int>(problem.GetNumberCitiesInCycle());
    std::vector<int> result = slicing(sol,numCitiesInCycle - 1);
    return result;
}

// ---------------- MAIN SOLVE LOOP ----------------

std::vector<int> LocalSearchSolver::solve() {
    std::vector<int> solution = initializeSolution();
    std::vector<int> unselected = problem.GiveIndices();

    int moveNum = 0;
    int lastCost = problem.FullDistanceAndCost(solution);
    std::string lastMoveType = "None";

    // Remove selected from unselected
    for (int idx : solution)
        unselected.erase(std::remove(unselected.begin(), unselected.end(), idx), unselected.end());

    for (int i = 0; i < (int)solution.size(); i++) {
        for (int j = 0; j < (int)solution.size(); j++) {
            if (i == j) continue;
            int delta = calculateDeltaIntraTwoEdge(solution, i, j);
            int a = solution[i];
            int b = solution[(i+1)%solution.size()];
            int c = solution[j];
            int d = solution[(j+1)%solution.size()];

            moveMemory_[StoredMove{MoveType::IntraEdge, a, c, b, d}] = delta;
            // moveMemory_.push_back(StoredMove{MoveType::IntraEdge, a, c, delta, b, d});
        }

        for (int j = 0; j < (int)unselected.size(); j++) {
            int delta = calculateDeltaInter(solution, i, j, unselected);
            int oldNode = solution[i];
            int newNode = unselected[j];

            moveMemory_[StoredMove{MoveType::InterNode, oldNode, newNode, -1, -1}] = delta;
            // moveMemory_.push_back(StoredMove{MoveType::InterNode, oldNode, newNode, delta});
        }
    }
    
    while (true) {
        moveNum++;

        MoveDelta move = findBestMove(solution, unselected);
        int i1prev = solution[move.i1 - 1 < 0 ? solution.size() - 1 : move.i1 - 1];
        int i1next = solution[(move.i1 + 1) % solution.size()];
        int i1 = solution[move.i1];
        int i2;
        int i2prev, i2next;
        if (move.type == MoveType::InterNode) {
            i2 = unselected[move.i2];
        } else {
            i2 = solution[move.i2];
            i2prev = solution[move.i2 - 1 < 0 ? solution.size() - 1 : move.i2 - 1];
            i2next = solution[(move.i2 + 1) % solution.size()];
        }

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

        // ---- After finding the best move, modify all the changed edges deltas in the memory ----
        std::unordered_map<int, int> solutionIndexMap;
        solutionIndexMap.reserve(solution.size());
        for (int idx = 0; idx < (int)solution.size(); idx++)
            solutionIndexMap[solution[idx]] = idx;
        
        std::unordered_map<int, int> unselectedIndexMap;
        unselectedIndexMap.reserve(unselected.size());
        for (int idx = 0; idx < (int)unselected.size(); idx++)
            unselectedIndexMap[unselected[idx]] = idx;

        std::unordered_map<StoredMove, int, StoredMoveHash> newMemory;
        newMemory.reserve(moveMemory_.size());
        
        if (move.type == MoveType::InterNode) {
            for (auto &m : moveMemory_) {
                StoredMove updatedMove = m.first;
             
                if (m.first.i1 == i1 || m.first.i2 == i2) {
                    if (updatedMove.i1 == i1)
                        updatedMove.i1 = i2;
                    
                    if (updatedMove.i2 == i2)
                        updatedMove.i2 = i1;
                    
                    if (m.first.type == MoveType::InterNode) {
                        newMemory[updatedMove] = calculateDeltaInter(solution, solutionIndexMap[updatedMove.i1], unselectedIndexMap[updatedMove.i2], unselected);
                        continue;
                        // updatedMove.delta = calculateDeltaInter(solution, std::distance(solution.begin(), std::find(solution.begin(), solution.end(), updatedMove.i1)), std::distance(unselected.begin(), std::find(unselected.begin(), unselected.end(), updatedMove.i2)), unselected);
                    } else {
                        // newMemory[updatedMove] = calculateDeltaIntraTwoEdge(solution, solutionIndexMap[updatedMove.i1], solutionIndexMap[updatedMove.i2]);
                        continue;
                        // updatedMove.delta = calculateDeltaIntraTwoEdge(solution, std::distance(solution.begin(), std::find(solution.begin(), solution.end(), updatedMove.i1)), std::distance(solution.begin(), std::find(solution.begin(), solution.end(), updatedMove.i2)));
                    }
                }
                else {
                    if (m.first.type == MoveType::InterNode) {
                        if (m.first.i1 == i1prev || m.first.i1 == i1next) {
                            newMemory[updatedMove] = m.second;
                            newMemory[updatedMove] = calculateDeltaInter(solution, solutionIndexMap[updatedMove.i1], unselectedIndexMap[updatedMove.i2], unselected);
                            continue;
                            // updatedMove.delta = calculateDeltaInter(solution, std::distance(solution.begin(), std::find(solution.begin(), solution.end(), updatedMove.i1)), std::distance(unselected.begin(), std::find(unselected.begin(), unselected.end(), updatedMove.i2)), unselected);
                        }
                    }
                    else {
                        if (m.first.i1 == i1next || m.first.i2 == i1next) {
                            // newMemory[updatedMove] = calculateDeltaIntraTwoEdge(solution, solutionIndexMap[updatedMove.i1], solutionIndexMap[updatedMove.i2]);
                            continue;
                            // updatedMove.delta = calculateDeltaIntraTwoEdge(solution, std::distance(solution.begin(), std::find(solution.begin(), solution.end(), updatedMove.i1)), std::distance(solution.begin(), std::find(solution.begin(), solution.end(), updatedMove.i2)));
                        }
                    }
                }

                newMemory[updatedMove] = m.second;
            }
        } else {
            for (auto &m : moveMemory_) {
                // Recreate the chosen move
                StoredMove updatedMove = m.first;

                if (m.first.type == MoveType::InterNode) {
                    if (m.first.i1 == i1 || m.first.i1 == i1next || m.first.i1 == i2 || m.first.i1 == i2next) {
                        newMemory[updatedMove] = calculateDeltaInter(solution, solutionIndexMap[updatedMove.i1], unselectedIndexMap[updatedMove.i2], unselected);
                    }
                    else {
                        newMemory[updatedMove] = m.second;
                    }
                }
                else {
                    if (!(solutionIndexMap.count(updatedMove.i1) && solutionIndexMap.count(updatedMove.i2) && solutionIndexMap.count(updatedMove.i1next) && solutionIndexMap.count(updatedMove.i2next)))
                        continue;

                    if (updatedMove.i1 == i1 && updatedMove.i1next == i1next && updatedMove.i2 == i2 && updatedMove.i2next == i2next) {
                        continue;
                    }

                    if (solutionIndexMap[updatedMove.i1] > solutionIndexMap[updatedMove.i1next] && solutionIndexMap[updatedMove.i2] > solutionIndexMap[updatedMove.i2next]) {
                        std::swap(updatedMove.i1, updatedMove.i1next);
                        std::swap(updatedMove.i2, updatedMove.i2next);
                    }

                    if (solutionIndexMap[updatedMove.i1] > solutionIndexMap[updatedMove.i2]) {
                        std::swap(updatedMove.i1, updatedMove.i2);
                        std::swap(updatedMove.i1next, updatedMove.i2next);
                    }

                    if (!adjacent(solutionIndexMap[updatedMove.i1], solutionIndexMap[updatedMove.i1next], solution.size()) || !adjacent(solutionIndexMap[updatedMove.i2], solutionIndexMap[updatedMove.i2next], solution.size())) {
                        continue;
                    }

                    // newMemory[updatedMove] = calculateDeltaIntraTwoEdge(solution, solutionIndexMap[updatedMove.i1], solutionIndexMap[updatedMove.i2]);
                    newMemory[updatedMove] = m.second;
                }
            }

            for (int i = 0; i < solution.size(); i++) {
                for (int j = i + 1; j < solution.size(); j++) {
                    StoredMove moveToCheck = StoredMove{MoveType::IntraEdge, solution[i], solution[j], solution[(i + 1) % solution.size()], solution[(j + 1) % solution.size()]};
                    if (!newMemory.count(moveToCheck)) {
                        int delta = calculateDeltaIntraTwoEdge(solution, i, j);
                        if (delta < INT_MAX)
                            newMemory[moveToCheck] = delta;
                    }
                }
            }
        }

        std::unordered_map<StoredMove, int, StoredMoveHash> checkMemory = newMemory;

        for (auto &m : checkMemory) {
            if (m.first.type == MoveType::InterNode) {
                if (!solutionIndexMap.count(m.first.i1) || !unselectedIndexMap.count(m.first.i2)) {
                    newMemory.erase(m.first);
                }
            }
            else {
                if (!solutionIndexMap.count(m.first.i1) || !solutionIndexMap.count(m.first.i2)) {
                    newMemory.erase(m.first);
                }
            }
        }

        moveMemory_.swap(newMemory);

        bool deltaCorrect = (problem.FullDistanceAndCost(solution) - lastCost == move.delta);

        if (!deltaCorrect) {
            std::cout << "Iteration " << moveNum << ", current cost: " << problem.FullDistanceAndCost(solution) << ", last cost: " << lastCost << ", move delta: " << move.delta << ", delta correct: " << (deltaCorrect ? "yes" : "no") << ", move type: " << (move.type == MoveType::IntraEdge ? "Intra Edge" : "Inter Node") << ", move: " << move.i1 << "-" << move.i2 << ", last move type: " << lastMoveType << std::endl;
        }

        lastCost = problem.FullDistanceAndCost(solution);
        lastMoveType = move.type == MoveType::IntraEdge ? "Intra Edge" : "Inter Node";
    }
    AssertHamiltonian(solution,problem.GetNumberCitiesInCycle());
    return solution;
}
