#include "LocalSearchSolver.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <valarray>
#include <iostream>
#include <sstream>
#include <vector>
#include <random>
#include <cassert>
#include <unordered_set>
#include <chrono>
#include <unordered_map>

// ---------------- Utility index helpers ----------------
inline int prevIdxInline(int i, int n) { return (i - 1 + n) % n; }
inline int nextIdxInline(int i, int n) { return (i + 1) % n; }

// ---------------- Constructor ----------------
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
    int n = (int)solution.size();
    if (n < 4) return INT_MAX;

    // skip adjacency
    if (j <= i+1 || ((i==0) && (j == n-1))) return INT_MAX;

    int u1 = solution[i];
    int v1 = solution[(i + 1) % n];
    int u2 = solution[j];
    int v2 = solution[(j + 1) % n];

    int oldCost = problem.GetDistance(u1, v1) + problem.GetDistance(u2, v2);
    int newCost = problem.GetDistance(u1, u2) + problem.GetDistance(v1, v2);

    return newCost - oldCost;
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

// ---------------- Internal helpers for LM generation ----------------

// Return position (index) of node in route, -1 if not present
int LocalSearchSolver::getPositionOfNode(const std::vector<int>& route, int node) const {
    for (int i = 0; i < (int)route.size(); ++i)
        if (route[i] == node) return i;
    return -1;
}

// Generate all 2-opt / edge-exchange moves across the whole route (explicit 4 variants)
void LocalSearchSolver::generateAllEdgeMoves(const std::vector<int>& route, const std::vector<std::vector<int>>& /*distMatrix*/, const std::vector<bool>& /*isSelected*/) {
    int n = (int)route.size();
    if (n < 4) return;

    // local lambda to add a variant
    auto tryAdd = [&](int a1, int b1, int a2, int b2) {
        StoredMove m = StoredMove::MakeEdgeExchange(a1, b1, a2, b2);
        int oldCost = problem.GetDistance(a1, b1) + problem.GetDistance(a2, b2);
        int newCost = problem.GetDistance(a1, a2) + problem.GetDistance(b1, b2);
        int delta = newCost - oldCost;
        if (delta < 0) {
            inLM_[m] = delta;
            LM_.push(PQItem(delta, m));
        }
    };

    for (int i = 0; i < n; ++i) {
        int u1 = route[i];
        int v1 = route[(i + 1) % n];
        for (int j = i + 2; j < n; ++j) {
            // skip adjacency
            if (j == (i + 1) % n) continue;

            int u2 = route[j];
            int v2 = route[(j + 1) % n];

            tryAdd(u1, v1, u2, v2);
            tryAdd(v1, u1, u2, v2);
            tryAdd(u1, v1, v2, u2);
            tryAdd(v1, u1, v2, u2);
        }
    }
}

// Generate all inter-route moves across the whole route
void LocalSearchSolver::generateAllInterMoves(const std::vector<int>& route, const std::vector<int>& unselected) {
    int n = (int)route.size();
    if (n == 0) return;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < (int)unselected.size(); ++j) {
            int delta = calculateDeltaInter(route, i, j, unselected);
            if (delta < 0) {
                StoredMove m = StoredMove::MakeInterNode(route[i], unselected[j]);
                inLM_[m] = delta;
                LM_.push(PQItem(delta, m));
            }
        }
    }
}

// Generate edge moves that involve any of the affectedEdges (node ID pairs).
void LocalSearchSolver::generateEdgeMovesForEdges(const std::vector<int>& route, const std::vector<std::pair<int,int>>& affectedEdges, const std::vector<bool>& /*isSelected*/) {
    int n = (int)route.size();
    if (n < 4) return;

    // Build mapping node->pos for quick lookup
    int maxNode = -1;
    for (int node : route) if (node > maxNode) maxNode = node;
    std::vector<int> posOf(maxNode + 1, -1);
    for (int i = 0; i < n; ++i) posOf[route[i]] = i;

    // local lambda to add all 4 variants
    auto tryAdd = [&](int a1, int b1, int a2, int b2) {
        StoredMove m = StoredMove::MakeEdgeExchange(a1, b1, a2, b2);
        int oldCost = problem.GetDistance(a1, b1) + problem.GetDistance(a2, b2);
        int newCost = problem.GetDistance(a1, a2) + problem.GetDistance(b1, b2);
        int delta = newCost - oldCost;
        if (delta < 0) {
            inLM_[m] = delta;
            LM_.push(PQItem(delta, m));
        }
    };

    for (const auto& e : affectedEdges) {
        int u1 = e.first;
        int v1 = e.second;
        // Ensure both nodes are selected and present in this route
        if (u1 < 0 || v1 < 0) continue;
        if (u1 >= (int)posOf.size() || v1 >= (int)posOf.size()) continue;
        int posU1 = posOf[u1];
        int posV1 = posOf[v1];
        if (posU1 == -1 || posV1 == -1) continue;

        // Pair this edge with all other edges in the route
        for (int i = 0; i < n; ++i) {
            int u2 = route[i];
            int v2 = route[(i + 1) % n];

            // Skip overlapping/adjacent edges
            if (u2 == u1 || u2 == v1 || v2 == u1 || v2 == v1) continue;

            // compute indices for delta function
            int idx1 = posU1;
            int idx2 = i;
            if (idx1 == idx2) continue;
            if (idx1 > idx2) std::swap(idx1, idx2);
            // skip adjacency
            if (idx2 == idx1 + 1 || (idx1 == 0 && idx2 == n - 1)) continue;

            tryAdd(u1, v1, u2, v2);
            tryAdd(v1, u1, u2, v2);
            tryAdd(u1, v1, v2, u2);
            tryAdd(v1, u1, v2, u2);
        }
    }
}

// Generate inter-route moves
void LocalSearchSolver::generateInterMovesForPositions(const std::vector<int>& route, const std::vector<int>& positions, const std::vector<int>& unselected, const std::vector<bool>& /*isSelected*/) {
    int n = (int)route.size();
    if (n == 0) return;

    for (int pos : positions) {
        if (pos < 0 || pos >= n) continue;
        for (int j = 0; j < (int)unselected.size(); ++j) {
            int delta = calculateDeltaInter(route, pos, j, unselected);
            if (delta < 0) {
                StoredMove m = StoredMove::MakeInterNode(route[pos], unselected[j]);
                inLM_[m] = delta;
                LM_.push(PQItem(delta, m));
            }
        }
    }
}

// Generate inter-route moves for a single newly-unselected node (nodeOut)
void LocalSearchSolver::generateInterMovesForUnselectedNode(const std::vector<int>& route, int nodeOut, const std::vector<int>& unselected, const std::vector<bool>& /*isSelected*/) {
    int n = (int)route.size();
    if (n == 0) return;

    // find index of nodeOut in unselected
    int unIdx = -1;
    for (int j = 0; j < (int)unselected.size(); ++j) {
        if (unselected[j] == nodeOut) { unIdx = j; break; }
    }
    if (unIdx == -1) return;

    for (int pos = 0; pos < n; ++pos) {
        int delta = calculateDeltaInter(route, pos, unIdx, unselected);
        if (delta < 0) {
            StoredMove m = StoredMove::MakeInterNode(route[pos], nodeOut);
            inLM_[m] = delta;
            LM_.push(PQItem(delta, m));
        }
    }
}

// ---------------- findBestMove ----------------
// This function pops from LM_ until it finds a move that is still applicable and improving.
MoveDelta LocalSearchSolver::findBestMove(const std::vector<int>& solution, const std::vector<int>& unselected) {
    MoveDelta best{-1, -1, INT_MAX, MoveType::InterNode};

    int n = (int)solution.size();
    if (n == 0) return best;

    // Build posOf table for nodes
    int maxNode = -1;
    for (int v : solution) if (v > maxNode) maxNode = v;
    for (int v : unselected) if (v > maxNode) maxNode = v;
    std::vector<int> posOf(maxNode + 1, -1);
    for (int i = 0; i < n; ++i) posOf[solution[i]] = i;

    // Build map for unselected node -> index
    std::unordered_map<int,int> unselectedIndex;
    for (int j = 0; j < (int)unselected.size(); ++j) unselectedIndex[unselected[j]] = j;

    std::vector<PQItem> tempMoves; // store moves that are reversed in orientation for later requeue

    while (!LM_.empty()) {
        PQItem top = LM_.top();
        LM_.pop();
        StoredMove m = top.move;

        // Validate popped PQ item against authoritative map inLM_
        auto it = inLM_.find(m);
        if (it == inLM_.end()) {
            // Stale: move no longer in map -> skip
            continue;
        }
        if (it->second != top.delta) {
            // Stale pq item (delta changed) -> skip
            continue;
        }

        inLM_.erase(it);

        if (m.type == MoveType::IntraEdge) {
            // check existence & direction for both edges
            if (m.u1 < 0 || m.v1 < 0 || m.u2 < 0 || m.v2 < 0) continue;
            if (m.u1 >= (int)posOf.size() || m.v1 >= (int)posOf.size() || m.u2 >= (int)posOf.size() || m.v2 >= (int)posOf.size()) {
                continue;
            }

            int pos_u1 = posOf[m.u1];
            int pos_v1 = posOf[m.v1];
            int pos_u2 = posOf[m.u2];
            int pos_v2 = posOf[m.v2];

            if (pos_u1 == -1 || pos_v1 == -1 || pos_u2 == -1 || pos_v2 == -1) continue;

            // check direction: 0=absent, 1=forward, 2=reversed
            auto checkDir = [&](int u, int v) -> int {
                int pu = posOf[u], pv = posOf[v];
                if (pu == -1 || pv == -1) return 0;
                if (nextIdxInline(pu, n) == pv) return 1;
                if (nextIdxInline(pv, n) == pu) return 2;
                return 0;
            };

            int dir1 = checkDir(m.u1, m.v1);
            int dir2 = checkDir(m.u2, m.v2);

            if (dir1 == 0 || dir2 == 0) {
                // edge missing -> discard
                continue;
            }
            if (dir1 == 2 || dir2 == 2) {
                // reversed orientation present -> postpone (requeue later)
                tempMoves.push_back(top);
                continue;
            }

            int idx1 = posOf[m.u1];
            int idx2 = posOf[m.u2];
            if (idx1 == -1 || idx2 == -1) continue;
            if (idx1 > idx2) std::swap(idx1, idx2);
            // adjacency check
            if (idx2 == idx1 + 1 || (idx1 == 0 && idx2 == n - 1)) {
                continue;
            }

            int delta = calculateDeltaIntraTwoEdge(solution, idx1, idx2);
            if (delta == INT_MAX) continue;
            if (delta >= 0) continue;

            best.i1 = idx1;
            best.i2 = idx2;
            best.delta = delta;
            best.type = MoveType::IntraEdge;

            // requeue tempMoves
            for (auto &tm : tempMoves) {
                inLM_[tm.move] = tm.delta;
                LM_.push(tm);
            }
            return best;
        }
        else { // InterNode
            // verify nodeIn still selected and nodeOut still unselected
            if (m.nodeIn < 0 || m.nodeOut < 0) continue;
            if (m.nodeIn >= (int)posOf.size()) continue;

            int pos = posOf[m.nodeIn];
            if (pos == -1) continue;

            if (!unselectedIndex.count(m.nodeOut)) {
                continue;
            }
            int unIdx = unselectedIndex[m.nodeOut];

            int delta = calculateDeltaInter(solution, pos, unIdx, unselected);
            if (delta >= 0) continue;

            best.i1 = pos;
            best.i2 = unIdx;
            best.delta = delta;
            best.type = MoveType::InterNode;

            // requeue tempMoves
            for (auto &tm : tempMoves) {
                inLM_[tm.move] = tm.delta;
                LM_.push(tm);
            }
            return best;
        }
    }

    for (auto &tm : tempMoves) {
        inLM_[tm.move] = tm.delta;
        LM_.push(tm);
    }
    return best;
}

// ---------------- STARTING SOLUTION ----------------

std::vector<int> LocalSearchSolver::initializeSolution() {
    std::vector<int> sol = problem.GiveIndices();
    std::shuffle(sol.begin(), sol.end(), rng_);
    int numCitiesInCycle = static_cast<int>(problem.GetNumberCitiesInCycle());
    std::vector<int> result = slicing(sol,numCitiesInCycle - 1);
    return result;
}

void LocalSearchSolver::SetStartingSol(std::vector<int>& ss) {
    startingSol = ss;
}

// ---------------- MAIN SOLVE LOOP ----------------

std::vector<int> LocalSearchSolver::solve() {
    std::vector<int> solution = initializeSolution();
    if (startingSol.size() == problem.GetNumberCitiesInCycle())
        solution = startingSol;
    std::vector<int> unselected = problem.GiveIndices();

    int moveNum = 0;
    int lastCost = problem.FullDistanceAndCost(solution);
    std::string lastMoveType = "None";

    // Remove selected from unselected
    for (int idx : solution)
        unselected.erase(std::remove(unselected.begin(), unselected.end(), idx), unselected.end());

    int nTotal = problem.getNumCities();
    int n = (int)solution.size();

    std::vector<bool> isSelected(nTotal, false);
    for (int node : solution) {
        if (node >= 0 && node < nTotal) isSelected[node] = true;
    }

    // Build posOf map (node id -> index in solution)
    std::vector<int> posOf(nTotal, -1);
    for (int i = 0; i < (int)solution.size(); ++i) posOf[solution[i]] = i;

    // Clear LM_ and inLM_ in case
    while (!LM_.empty()) LM_.pop();
    inLM_.clear();

    // Pre-generate all moves
    std::vector<std::vector<int>> dummyDist;
    generateAllEdgeMoves(solution, dummyDist, isSelected);
    generateAllInterMoves(solution, unselected);

    while (true) {
        moveNum++;

        MoveDelta move = findBestMove(solution, unselected);
        if (move.delta >= 0) break;

        // Prepare pre-move endpoints for correct new-edge identification
        int apply_i1 = move.i1;
        int apply_i2 = move.i2;
        int i1prev = solution[apply_i1 - 1 < 0 ? solution.size() - 1 : apply_i1 - 1];
        int i1next = solution[(apply_i1 + 1) % solution.size()];
        int i1 = solution[apply_i1];
        int i2;
        int i2prev = -1, i2next = -1;

        // For IntraEdge
        int u1 = -1, v1 = -1, u2 = -1, v2 = -1;
        if (move.type == MoveType::IntraEdge) {
            // positions in solution
            int posA = move.i1;
            int posB = move.i2;
            int size = (int)solution.size();
            u1 = solution[posA];
            v1 = solution[(posA + 1) % size];
            u2 = solution[posB];
            v2 = solution[(posB + 1) % size];
        }

        if (move.type == MoveType::InterNode) {
            i2 = unselected[move.i2];
        } else {
            i2 = solution[move.i2];
            i2prev = solution[move.i2 - 1 < 0 ? solution.size() - 1 : move.i2 - 1];
            i2next = solution[(move.i2 + 1) % solution.size()];
        }

        // Apply the move
        switch (move.type) {
            case MoveType::InterNode:
                std::swap(solution[move.i1], unselected[move.i2]);

                // Update isSelected and posOf
                isSelected[i1] = false;
                isSelected[i2] = true;

                // update posOf mapping
                posOf[i2] = move.i1;
                posOf[i1] = -1;
                break;

            case MoveType::IntraEdge:
                if (move.i1 > move.i2) std::swap(move.i1, move.i2);
                std::reverse(solution.begin() + move.i1 + 1, solution.begin() + move.i2 + 1);

                // rebuild posOf
                for (int idx = 0; idx < (int)solution.size(); ++idx) posOf[solution[idx]] = idx;
                break;
        }

        // After applying the move, generate new moves only for affected edges/positions
        if (move.type == MoveType::InterNode) {
            int pos = move.i1;
            int prevNode = (pos == 0) ? solution.back() : solution[pos - 1];
            int nextNode = (pos == (int)solution.size() - 1) ? solution.front() : solution[pos + 1];
            int newNode = solution[pos];

            std::vector<std::pair<int,int>> newEdges;
            newEdges.push_back({prevNode, newNode});
            newEdges.push_back({newNode, nextNode});

            std::vector<int> affectedPositions = { posOf[prevNode], posOf[newNode], posOf[nextNode] };

            generateEdgeMovesForEdges(solution, newEdges, isSelected);
            generateInterMovesForPositions(solution, affectedPositions, unselected, isSelected);

            // Also generate inter-route moves that use the node that was removed from the route (now unselected)
            // i1 is the node removed
            generateInterMovesForUnselectedNode(solution, i1, unselected, isSelected);
        } else {
            std::vector<std::pair<int,int>> newEdges;
            if (u1 != -1 && v1 != -1 && u2 != -1 && v2 != -1) {
                newEdges.push_back({u1, u2});
                newEdges.push_back({v1, v2});
            }

            std::vector<int> affectedPositions;
            if (u1 != -1) affectedPositions.push_back(posOf[u1]);
            if (u2 != -1) affectedPositions.push_back(posOf[u2]);
            if (v1 != -1) affectedPositions.push_back(posOf[v1]);
            if (v2 != -1) affectedPositions.push_back(posOf[v2]);

            generateEdgeMovesForEdges(solution, newEdges, isSelected);
            generateInterMovesForPositions(solution, affectedPositions, unselected, isSelected);
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