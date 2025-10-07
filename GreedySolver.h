#pragma once
#include <vector>
#include "Solver.h"

// Enum 
enum class GreedyMode {
    NearestNeighbour,
    NearestNeighbourEnd,
    GreedyCycle
};

// GreedySolver inherits from Solver
class GreedySolver : public Solver {
public:
    explicit GreedySolver(const ProblemInstance& prob, GreedyMode mode = GreedyMode::NearestNeighbour);


    std::vector<int> solve() override;

private:
    GreedyMode mode_;
    // Todo implement the greedy methods 
    std::vector<int> solveNearestNeighbour();
    std::vector<int> solveNearestNeighbourEnd();
    std::vector<int> solveGreedyCycle();
};
