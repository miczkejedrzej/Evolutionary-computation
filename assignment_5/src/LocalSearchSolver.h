#pragma once
#include "../../src/ProblemInstance.h"
#include "../../src/Solver.h"
#include <vector>
#include <random>
#include <algorithm>
#include <climits>
#include <unordered_set>
#include <unordered_map>

enum class MoveType { InterNode, IntraEdge };

struct Node_List;

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
    struct StoredMove {
        MoveType type;
        int i1, i2;
        int i1next, i2next;

        bool operator==(const StoredMove& other) const {
            return type == other.type &&
                i1 == other.i1 &&
                i2 == other.i2 &&
                i1next == other.i1next &&
                i2next == other.i2next;
        }
    };

    struct StoredMoveHash {
        std::size_t operator()(const StoredMove& m) const noexcept {
            std::size_t h1 = std::hash<int>()(static_cast<int>(m.type));
            std::size_t h2 = std::hash<int>()(m.i1);
            std::size_t h3 = std::hash<int>()(m.i2);
            std::size_t h4 = std::hash<int>()(m.i1next);
            std::size_t h5 = std::hash<int>()(m.i2next);

            // Combine hashes — standard pattern
            std::size_t seed = h1;
            seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= h4 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= h5 + 0x9e3779b9 + (seed << 6) + (seed >> 2);

            return seed;
        }
    };

    std::unordered_map<StoredMove, int, StoredMoveHash> moveMemory_;

    std::mt19937 rng_;

    std::vector<int> slicing(const std::vector<int>& arr, int X);

    std::vector<int> initializeSolution();

    int calculateDeltaInter(const std::vector<int>& solution, int selIdx, int unselIdx, const std::vector<int>& unselected);
    int calculateDeltaIntraTwoEdge(const std::vector<int>& solution, int i, int j);
    MoveDelta findBestMove(const std::vector<int>& solution, const std::vector<int>& unselected);
    void AssertHamiltonian(const std::vector<int>& visited, int citiesNumber);

    void fillCandidateEdges();
};