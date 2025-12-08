#pragma once
#include "../../src/ProblemInstance.h"
#include "../../src/Solver.h"
#include <vector>
#include <random>
#include <algorithm>
#include <climits>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <functional>
#include <string>

enum class MoveType { InterNode, IntraEdge };

struct Node_List;

struct StoredMove {
    MoveType type;

    // For IntraEdge (EdgeExchange / 2-opt) moves: edges u1->v1 and u2->v2 (node IDs)
    int u1{-1}, v1{-1}, u2{-1}, v2{-1};

    // For InterNode (InterRoute) moves:
    int nodeIn{-1}, nodeOut{-1};

    // Default constructor
    StoredMove() = default;

    // Convenience constructors
    static StoredMove MakeEdgeExchange(int U1, int V1, int U2, int V2) {
        StoredMove m;
        m.type = MoveType::IntraEdge;
        m.u1 = U1; m.v1 = V1; m.u2 = U2; m.v2 = V2;
        return m;
    }

    static StoredMove MakeInterNode(int nodeIn, int nodeOut) {
        StoredMove m;
        m.type = MoveType::InterNode;
        m.nodeIn = nodeIn;
        m.nodeOut = nodeOut;
        return m;
    }

    bool operator==(const StoredMove& other) const {
        if (type != other.type) return false;
        if (type == MoveType::IntraEdge) {
            auto n1a = std::min(u1, v1), n1b = std::max(u1, v1);
            auto n2a = std::min(u2, v2), n2b = std::max(u2, v2);
            // Ensure a canonical order between the two edges
            if (n1a > n2a || (n1a == n2a && n1b > n2b)) {
                std::swap(n1a, n2a);
                std::swap(n1b, n2b);
            }
            auto o1a = std::min(other.u1, other.v1), o1b = std::max(other.u1, other.v1);
            auto o2a = std::min(other.u2, other.v2), o2b = std::max(other.u2, other.v2);
            if (o1a > o2a || (o1a == o2a && o1b > o2b)) {
                std::swap(o1a, o2a);
                std::swap(o1b, o2b);
            }
            return n1a == o1a && n1b == o1b && n2a == o2a && n2b == o2b;
        } else {
            return nodeIn == other.nodeIn && nodeOut == other.nodeOut;
        }
    }
};

// Hash for StoredMove
struct StoredMoveHash {
    std::size_t operator()(const StoredMove& m) const noexcept {
        if (m.type == MoveType::IntraEdge) {
            int a1 = std::min(m.u1, m.v1), b1 = std::max(m.u1, m.v1);
            int a2 = std::min(m.u2, m.v2), b2 = std::max(m.u2, m.v2);
            // canonical edge order
            if (a1 > a2 || (a1 == a2 && b1 > b2)) {
                std::swap(a1, a2);
                std::swap(b1, b2);
            }
            std::size_t h = 1469598103934665603u;
            h ^= std::hash<int>{}(a1) + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
            h ^= std::hash<int>{}(b1) + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
            h ^= std::hash<int>{}(a2) + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
            h ^= std::hash<int>{}(b2) + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
            return h;
        } else {
            std::size_t seed = std::hash<int>{}(m.nodeIn);
            seed ^= std::hash<int>{}(m.nodeOut) + 0x9e3779b9 + (seed<<6) + (seed>>2);
            return seed;
        }
    }
};

// Equality for unordered_set / unordered_map
struct StoredMoveEqual {
    bool operator()(const StoredMove& a, const StoredMove& b) const noexcept {
        return a == b;
    }
};

// Item for the priority queue: stores delta (double) and the move
struct PQItem {
    int delta;
    StoredMove move;

    PQItem() = default;
    PQItem(int d, const StoredMove& m) : delta(d), move(m) {}

    // Invert comparison here for correct ordering.
    bool operator<(const PQItem& other) const {
        return delta > other.delta;
    }
};

struct MoveDelta {
    int i1;
    int i2;
    int delta;
    MoveType type;
};

class LocalSearchSolver : public Solver {
public:
    LocalSearchSolver(const ProblemInstance& prob,
                      int randomSeed = 42);

    std::vector<int> solve();

private:
    // The LM (list of moves) as a priority queue with their deltas
    std::priority_queue<PQItem> LM_;
    // Map move -> current delta stored for that move (used to validate PQ entries and allow updates)
    std::unordered_map<StoredMove, int, StoredMoveHash, StoredMoveEqual> inLM_;

    std::mt19937 rng_;

    // Helper functions
    std::vector<int> slicing(const std::vector<int>& arr, int X);

    std::vector<int> initializeSolution();

    int calculateDeltaInter(const std::vector<int>& solution, int selIdx, int unselIdx, const std::vector<int>& unselected);
    int calculateDeltaIntraTwoEdge(const std::vector<int>& solution, int i, int j);
    MoveDelta findBestMove(const std::vector<int>& solution, const std::vector<int>& unselected);
    void AssertHamiltonian(const std::vector<int>& visited, int citiesNumber);

    // Helpers used by the solver
    int getPositionOfNode(const std::vector<int>& route, int node) const;
    void generateAllEdgeMoves(const std::vector<int>& route, const std::vector<std::vector<int>>& distMatrix, const std::vector<bool>& isSelected);
    void generateAllInterMoves(const std::vector<int>& route, const std::vector<int>& unselected);
    void generateEdgeMovesForEdges(const std::vector<int>& route, const std::vector<std::pair<int,int>>& affectedEdges, const std::vector<bool>& isSelected);
    void generateInterMovesForPositions(const std::vector<int>& route, const std::vector<int>& positions, const std::vector<int>& unselected, const std::vector<bool>& isSelected);
    void generateInterMovesForUnselectedNode(const std::vector<int>& route, int nodeOut, const std::vector<int>& unselected, const std::vector<bool>& isSelected);

    void fillCandidateEdges();
};
