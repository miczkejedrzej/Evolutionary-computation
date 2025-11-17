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
#include <array>

LocalSearchSolver::LocalSearchSolver(const ProblemInstance& prob,
                                     int randomSeed)
    : Solver(prob),
      rng_(randomSeed)
{}

// return previous index modulo n
inline int prevIdx(int idx, int n) { return (idx == 0) ? n - 1 : idx - 1; }
inline int nextIdx(int idx, int n) { return (idx + 1 == n) ? 0 : idx + 1; }

// collect indices from start to end inclusive, walking forward with wrap
static std::vector<int> collectRangeInclusive(int start, int end, int n) {
    std::vector<int> v;
    if (n == 0) return v;
    int k = start;
    while (true) {
        v.push_back(k);
        if (k == end) break;
        k = (k + 1) % n;
    }
    return v;
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
    // ---- Choose steepest descent move ----
    MoveDelta best{-1, -1, 0, MoveType::InterNode};
    best.delta = INT_MAX;

    int n = solution.size();
    int m = unselected.size();

    // ----- Inter-node -----
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            if (deltaInter[i][j] < best.delta)
                best = {i, j, deltaInter[i][j], MoveType::InterNode};

    // ----- 2-opt -----
    for (int i = 0; i < n; i++)
        for (int j = i + 2; j < n; j++)
            if (delta2opt[i][j] < best.delta)
                best = {i, j, delta2opt[i][j], MoveType::IntraEdge};

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

    int n = solution.size();
    int m = unselected.size();

    deltaInter.assign(n, std::vector<int>(m, INT_MAX));
    delta2opt.assign(n, std::vector<int>(n, INT_MAX));

    int moveNum = 0;
    int lastCost = problem.FullDistanceAndCost(solution);
    std::string lastMoveType = "None";

    // Remove selected from unselected
    for (int idx : solution)
        unselected.erase(std::remove(unselected.begin(), unselected.end(), idx), unselected.end());

    // ------------ BUILD INITIAL DELTAS ----------------
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            deltaInter[i][j] = calculateDeltaInter(solution, i, j, unselected);
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = i + 2; j < n; j++) {   // skip adjacency
            delta2opt[i][j] = calculateDeltaIntraTwoEdge(solution, i, j);
        }
    }
    
    while (true) {
        moveNum++;

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

        // ---- After finding the best move, modify all the changed edges deltas in the memory ----
        std::unordered_map<int, int> solutionIndexMap;
        solutionIndexMap.reserve(solution.size());
        for (int idx = 0; idx < (int)solution.size(); idx++)
            solutionIndexMap[solution[idx]] = idx;
        
        std::unordered_map<int, int> unselectedIndexMap;
        unselectedIndexMap.reserve(unselected.size());
        for (int idx = 0; idx < (int)unselected.size(); idx++)
            unselectedIndexMap[unselected[idx]] = idx;
        
        if (move.type == MoveType::InterNode) {
            int p = move.i1;           // position in solution swapped
            int ucol = move.i2;        // column in unselected that was swapped

            // neighbors in tour (indices)
            int A = prevIdx(p, n);
            int B = nextIdx(p, n);

            // 1) Column ucol changed (unselected[u_col] changed identity), so update this column for all selIdx
            for (int sel = 0; sel < n; ++sel) {
                deltaInter[sel][ucol] = calculateDeltaInter(solution, sel, ucol, unselected);
            }

            // 2) The nodes whose prev/next changed: A, p, B.
            // update their rows (they depend on their prev/next)
            std::array<int,3> rowsToUpdate = {A, p, B};
            for (int sel : rowsToUpdate) {
                for (int uj = 0; uj < m; ++uj) {
                    deltaInter[sel][uj] = calculateDeltaInter(solution, sel, uj, unselected);
                }
            }

            // 3) update 2-opt deltas for pairs where either endpoint is one of A,p,B.
            // For each s in {A,p,B} and for every other index t (non-adjacent), recompute delta2opt[min,max].
            for (int s : rowsToUpdate) {
                for (int t = 0; t < n; ++t) {
                    if (t == s) continue;
                    // skip adjacency (including wrap)
                    if (t == nextIdx(s, n) || t == prevIdx(s, n)) continue;
                    int x = std::min(s, t), y = std::max(s, t);
                    delta2opt[x][y] = calculateDeltaIntraTwoEdge(solution, x, y);
                }
            }
        } else { // IntraEdge
            // assume move.i1 <= move.i2 (you already ensured that earlier)
            int i = move.i1;
            int j = move.i2;

            // affected tour indices: the whole reversed segment [i..j] and their immediate neighbors i-1 and j+1
            std::vector<int> affected = collectRangeInclusive(i, j, n);
            int leftNeighbor = prevIdx(i, n);
            int rightNeighbor = nextIdx(j, n);
            // include neighbors (they are affected because their prev/next may change)
            affected.push_back(leftNeighbor);
            affected.push_back(rightNeighbor);

            // make unique (in case of small segments where neighbors overlap)
            std::sort(affected.begin(), affected.end());
            affected.erase(std::unique(affected.begin(), affected.end()), affected.end());

            // 1) Update deltaInter rows for all affected sel indices
            for (int sel : affected) {
                for (int uj = 0; uj < m; ++uj) {
                    deltaInter[sel][uj] = calculateDeltaInter(solution, sel, uj, unselected);
                }
            }

            // 2) Update delta2opt entries where either endpoint is in 'affected'
            // For each a in affected, and for each b not equal to a and not adjacent, recalc.
            for (int a : affected) {
                for (int b = 0; b < n; ++b) {
                    if (b == a) continue;
                    // skip adjacent pairs because those are invalid for 2-opt
                    if (b == nextIdx(a, n) || b == prevIdx(a, n)) continue;
                    int x = std::min(a, b);
                    int y = std::max(a, b);
                    // ensure we only compute upper triangle [x][y] once; setting repeatedly is fine
                    delta2opt[x][y] = calculateDeltaIntraTwoEdge(solution, x, y);
                }
            }
        }

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
