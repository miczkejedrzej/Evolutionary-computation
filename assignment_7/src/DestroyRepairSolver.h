#pragma once
#include <vector>
#include <memory>
#include <random>
#include "../../src/ProblemInstance.h"
#include "../../assignment_3/src/LocalSearchSolver.h"
#include "../../assignment_2/src/GreedySolver.h"
#include "../../assignment_6/src/MlsSolver.h"
enum class DestroyHeuristic
{
    None,
    Node,
    Mixed
};

// struct Result
// {
//     int64_t best_cost;
//     std::vector<int> best_solution;
//     double elapsed_ms; // total time for run
//     int ls_runs;       // number of iterations
// };

class DestroyRepairSolver
{
public:
    DestroyRepairSolver(const ProblemInstance &instance, bool useLocalSearch, double destroyFraction,
                        DestroyHeuristic heuristic, double timeLimit);

    Result run_once(double timeLimit, int seed = 0);
    std::vector<Result> run_multiple(int runs = 20);

    void initializeGreedySolver(int start_index, GreedyMode mode, Heuristic heuristic, float weight);

private:
    const ProblemInstance &problem;
    bool useLocalSearch;
    double destroyFraction;
    DestroyHeuristic heuristic;
    MoveType intra_;
    std::unique_ptr<GreedySolver> greedy_solver_; // <-- unique_ptr now
    double timeLimit;

    std::vector<int> destroy(const std::vector<int> &solution,
                             std::mt19937 &rng,
                             DestroyHeuristic heuristic);
    std::vector<int> repair(const std::vector<int> &partial_solution);
};
