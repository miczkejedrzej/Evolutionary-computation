#pragma once
#include "../../src/ProblemInstance.h"
#include "../../src/Solver.h"
#include "../../assignment_2/src/GreedySolver.h"
#include <vector>
#include <random>
#include <algorithm>
#include <climits>

enum class LocalSearchType { Steepest, GreedyRandom };
enum class MoveType { InterRoute, IntraNodeSwap, IntraEdgeSwap };
enum class StartSolutionType { Random, Greedy, Given };

struct MoveDelta {
    int i1;
    int i2;
    int delta;
    MoveType type;
};

class LocalSearchSolver : public Solver {
public:
    LocalSearchSolver(const ProblemInstance& prob,
                      LocalSearchType searchType,
                      MoveType intraType,
                      StartSolutionType startType,
                      int startingIndex,
                      int randomSeed = 42);

    std::vector<int> solve();
    void setGivenSolution(const std::vector<int>& solution);
    std::vector<int> solve_from_solution();
    

private:
    LocalSearchType searchType_;
    MoveType intraType_;
    StartSolutionType startType_;
    std::mt19937 rng_;
    int startingIndex;
    std::vector<int> givenSolution;
  

    std::vector<int> slicing(const std::vector<int>& arr, int X);

    std::vector<int> initializeSolution();

    int calculateDeltaInter(const std::vector<int>& solution, int selIdx, int unselIdx, const std::vector<int>& unselected);
    int calculateDeltaIntraTwoNode(const std::vector<int>& solution, int i, int j);
    int calculateDeltaIntraTwoEdge(const std::vector<int>& solution, int i, int j);
    int giveStartingIndex();
    MoveDelta findBestMove(const std::vector<int>& solution, const std::vector<int>& unselected);
    MoveDelta findRandomGreedyMove(const std::vector<int>& solution, const std::vector<int>& unselected);
    void AssertHamiltonian(const std::vector<int>& visited, int citiesNumber);
};
