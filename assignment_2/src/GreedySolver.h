#pragma once
#include <vector>
#include <tuple>
#include <queue>
#include "../../src/Solver.h"

// Enum 
enum class GreedyMode {
    NearestNeighbour,
    GreedyCycle
};
enum class Heuristic {
    Regret,
    HybridRegretObjective
};



struct Node_List;

class GreedySolver : public Solver {
public:
    explicit GreedySolver(const ProblemInstance& prob, int index, GreedyMode mode = GreedyMode::NearestNeighbour, 
        Heuristic heuristic = Heuristic::Regret, float weightObjective=0.5);

    std::vector<int> solve() override;

    int getStartingIndex();

    void setWeight(float weight);

    float getWeight();

    std::vector<int> complete_solution(std::vector<int> partial_solution);

private:
    GreedyMode mode_;
    Heuristic heuristic;
    int starting_index;
    void AssertHamiltonian(const std::vector<int>& visited, int citiesNumber);
    int ObjectiveEvaluation(const Node_List& node_list, Heuristic heuristic);
    std::vector<int> solveElementNearestNeighbour(int targetSize, const std::vector<int>& visited, const std::vector<int>& unvisited, Heuristic heuristic);
    std::vector<int> solveElementGreedyCycle(int targetSize, const std::vector<int>& visited, const std::vector<int>& unvisited, Heuristic heuristic);
    float weight;
};
