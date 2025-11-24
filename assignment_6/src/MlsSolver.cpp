#pragma once
#include <vector>
#include <random>
#include <chrono>
#include <string>
#include "../../src/ProblemInstance.h"
#include "LocalSearchSolver.h"

struct MSLSResult {
    int64_t best_cost;
    std::vector<int> best_solution;
    double elapsed_ms; // time for the entire MSLS run
    int ls_runs; // number of basic LS runs performed (should equal restarts)
};

class MSLS {
public:
    MSLS(const ProblemInstance& prob, int restarts = 200, int seed = 12345, MoveType intra = MoveType::IntraNodeSwap);

    // Single MSLS run: perform 'restarts' LS runs and return best
    MSLSResult run_once();

    // Repeat 'runs' times and return vector of results
    std::vector<MSLSResult> run_multiple(int runs = 20);

private:
    const ProblemInstance& prob_;
    int restarts_;
    int seed_;
    MoveType intra_;
};
