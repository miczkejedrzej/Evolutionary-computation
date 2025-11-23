#include "MlsSolver.h"
#include <limits>
#include <iostream>

MlsSolver::MlsSolver(const ProblemInstance& prob, int restarts, int seed, MoveType intra)
    : prob_(prob), restarts_(restarts), seed_(seed), intra_(intra) {}

Result MlsSolver::run_once(int run_seed) {
    auto t0 = std::chrono::high_resolution_clock::now();

    int64_t best_cost = std::numeric_limits<int64_t>::max();
    std::vector<int> best_solution;
    int ls_runs = 0;

    // Per assignment: use random starting solutions and steepest local search
    for (int i = 0; i < restarts_; ++i) {
        // startingIndex: we can pass 0 since startType Random will shuffle inside solver
        LocalSearchSolver solver(prob_, LocalSearchType::Steepest, intra_, StartSolutionType::Random, run_seed,i );
        auto sol = solver.solve();
        int64_t cost = prob_.FullDistanceAndCost(sol);
        ++ls_runs;
        if (cost < best_cost) {
            best_cost = cost;
            best_solution = sol;
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    double elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    return Result{best_cost, best_solution, elapsed_ms, ls_runs};
}

std::vector<Result> MlsSolver::run_multiple(int runs) {
    std::vector<Result> results;
    results.reserve(runs);
    for (int i = 0; i < runs; ++i) {
        Result r = run_once(i);
        results.push_back(r);
        std::cout << "[MSLS] run " << (i+1) << "/" << runs << " cost=" << r.best_cost
                  << " time(ms)=" << r.elapsed_ms << " ls_runs=" << r.ls_runs << std::endl;
    }
    return results;
}
