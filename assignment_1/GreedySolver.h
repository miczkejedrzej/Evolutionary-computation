#pragma once
#include <vector>
#include <tuple>
#include <queue>
#include "..\template_classes\Solver.h"

// Enum 
enum class GreedyMode {
    NearestNeighbour,
    NearestNeighbourEnd,
    GreedyCycle
};



// GreedySolver inherits from Solver
class GreedySolver : public Solver {
public:
    explicit GreedySolver(const ProblemInstance& prob, int index, GreedyMode mode = GreedyMode::NearestNeighbour);

    std::vector<int> solve() override;

    int getStartingIndex();

private:
    GreedyMode mode_;
    // Todo implement the greedy methods 
    std::vector<int> solveNearestNeighbour();
    std::vector<int> solveNearestNeighbourEnd();
    std::vector<int> solveGreedyCycle();
    void AssertHamiltonian(std::vector<int> visited,int citiesNumber);
    int starting_index;
};
