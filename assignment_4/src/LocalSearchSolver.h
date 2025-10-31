#pragma once
#include "../../src/ProblemInstance.h"
#include "../../src/Solver.h"
#include <vector>
#include <random>
#include <algorithm>
#include <climits>

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
    std::vector<std::vector<int>> candidateEdges;

    std::mt19937 rng_;

    std::vector<int> slicing(const std::vector<int>& arr, int X);

    std::vector<int> initializeSolution();

    int calculateDeltaInter(const std::vector<int>& solution, int selIdx, int unselIdx, const std::vector<int>& unselected);
    int calculateDeltaIntraTwoEdge(const std::vector<int>& solution, int i, int j);
    MoveDelta findBestMove(const std::vector<int>& solution, const std::vector<int>& unselected);
    void AssertHamiltonian(const std::vector<int>& visited, int citiesNumber);

    void fillCandidateEdges();
};