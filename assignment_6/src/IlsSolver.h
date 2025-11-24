#pragma once
#include <vector>
#include <random>
#include <chrono>
#include <string>
#include <unordered_set>
#include "../../src/ProblemInstance.h"
#include "../../assignment_3/src/LocalSearchSolver.h"
#include "MlsSolver.h"

// struct Result {
//     int64_t best_cost;
//     std::vector<int> best_solution;
//     double elapsed_ms;
//     int ls_runs;
// };

class IlsSolver {
public:
    IlsSolver(const ProblemInstance& prob, double time_budget_ms, int seed = 7777, MoveType intra = MoveType::IntraEdgeSwap, int perturbation_strength = 10, int gene_pool_extension_size = 10, bool use_double_bridge = false); 

    // Run a single ILS (returns best within time budget)
    Result run_once();

    // Repeat multiple times
    std::vector<Result> run_multiple(int runs = 20);

private:
    const ProblemInstance& prob_;
    double time_budget_ms_;
    int seed_;
    MoveType intra_;
    int perturbation_strength_;
    int gene_pool_extension_size_;
    bool use_double_bridge_;

    // perturbation helpers
    std::vector<int> double_bridge(const std::vector<int>& tour, std::mt19937& rng);
    // std::vector<int> scramble(const std::vector<int>& tour, std::mt19937& rng, int swaps = 10);

    // canonicalization for visited detection
    std::vector<int> extend_gene_pool(const std::vector<int>& tour, const ProblemInstance& prob, int desired_size);
    std::vector<int> perturb(const std::vector<int>& tour,const ProblemInstance& problem, std::mt19937& rng,int number_perturbations);
};
