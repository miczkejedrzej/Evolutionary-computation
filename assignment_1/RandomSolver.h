#pragma once
#include "../template_classes/ProblemInstance.h"
#include "../template_classes/Solver.h"
#include <vector>
#include <random>
class RandomSolver: public Solver {
public:
    explicit RandomSolver(const ProblemInstance& prob);
    std::vector<int> solve() override;
private:
    static std::vector<int> slicing(const std::vector<int> arr, int X);
};
