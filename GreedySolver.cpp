#include <vector>
#include <iostream>
#include "Solver.h"

enum class GreedyMode {
    NearestNeighbour,
    NearestNeighbourEnd,
    GreedyCycle
};

class GreedySolver : public Solver {
public:
    GreedySolver(const ProblemInstance& prob, GreedyMode mode = GreedyMode::NearestNeighbour):
    Solver(prob), mode_(mode) {}

    std::vector<int> solve() override {
        switch (mode_) {
            case GreedyMode::NearestNeighbour:
                return solveNearestNeighbour();
            case GreedyMode::NearestNeighbourEnd:
                return solveNearestNeighbourEnd();
            case GreedyMode::GreedyCycle:
                return solveGreedyCycle();
            default:
                return {};
        }
    }

private:
    GreedyMode mode_;

    std::vector<int> solveNearestNeighbour() {
        std::vector<int> result;
        std::cout << "Solving with Nearest Neighbour...\n";
        // TODO: implement
        return result;
    }

    std::vector<int> solveNearestNeighbourEnd() {
        std::vector<int> result;
        std::cout << "Solving with Nearest Neighbour from End...\n";
        // TODO: implement
        return result;
    }

    std::vector<int> solveGreedyCycle() {
        std::vector<int> result;
        std::cout << "Solving with Greedy Cycle...\n";
        // TODO: implement
        return result;
    }
};
