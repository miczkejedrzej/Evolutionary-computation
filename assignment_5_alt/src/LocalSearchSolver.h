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

    struct LMMove {
        MoveType type;        // InterNode or IntraEdge (2-opt)
        // stored as node IDs (not indices)
        int a;  // for Inter: solution node replaced; for 2-opt: first endpoint
        int b;  // for Inter: unselected node inserted; for 2-opt: second endpoint
        // For 2-opt, also store their next nodes in the stored orientation (used for edge checks)
        int a_next;
        int b_next;

        int delta;            // cached delta for this move (when applicable)
        bool alive;           // not removed permanently
        bool applicableNow;   // if edges exist in same relative direction (ready to apply)
        // note: if edges exist but reversed relative direction, alive==true but applicableNow==false

        LMMove() = default;
        LMMove(MoveType t, int aa, int bb, int aan, int bbn, int d)
            : type(t), a(aa), b(bb), a_next(aan), b_next(bbn), delta(d), alive(true), applicableNow(false) {}
    };

    std::vector<LMMove> moves; // local memory
    std::unordered_map<long long, std::vector<int>> edgeToMoves; // map directed edge -> indices in moves
    std::vector<int> posInSol; // posInSol[nodeID] -> index in solution

    // Inter-node: deltaInter[i][j] = gain of replacing solution[i] with unselected[j]
    std::vector<std::vector<int>> deltaInter;

    // Intra-edge (2-opt): delta2opt[i][j] = gain of reversing solution[i..j]
    std::vector<std::vector<int>> delta2opt;

    std::mt19937 rng_;

    std::vector<int> slicing(const std::vector<int>& arr, int X);

    std::vector<int> initializeSolution();

    int calculateDeltaInter(const std::vector<int>& solution, int selIdx, int unselIdx, const std::vector<int>& unselected);
    int calculateDeltaIntraTwoEdge(const std::vector<int>& solution, int i, int j);
    MoveDelta findBestMove(const std::vector<int>& solution, const std::vector<int>& unselected);
    void AssertHamiltonian(const std::vector<int>& visited, int citiesNumber);

    void fillCandidateEdges();
};