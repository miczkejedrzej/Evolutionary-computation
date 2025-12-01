#include "DestroyRepairSolver.h"
#include <algorithm>
#include <random>
#include <chrono>
#include <iostream>
#include <cassert>

DestroyRepairSolver::DestroyRepairSolver(const ProblemInstance &problem, bool useLocalSearch, double destroyFraction,
                                         DestroyHeuristic heuristic, double timeLimit)
    : problem(problem), useLocalSearch(useLocalSearch), destroyFraction(destroyFraction),
      heuristic(heuristic), timeLimit(timeLimit)
{
}

// -------------------- Run multiple --------------------
std::vector<Result> DestroyRepairSolver::run_multiple(int runs)
{
    std::vector<Result> results;
    results.reserve(runs);
    for (int i = 0; i < runs; ++i)
    {
        Result r = run_once(timeLimit, i);
        results.push_back(r);
        std::cout << "destroy repair run " << (i + 1) << "/" << runs
                  << " cost=" << r.best_cost
                  << " time(ms)=" << r.elapsed_ms
                  << " ls_runs=" << r.ls_runs << std::endl;
    }
    return results;
}

// -------------------- Run once --------------------

Result DestroyRepairSolver::run_once(double timeLimit, int seed)
{
    initializeGreedySolver(0, GreedyMode::NearestNeighbour, Heuristic::HybridRegretObjective, 0.5);

    auto t_start = std::chrono::high_resolution_clock::now();
    std::mt19937 rng(seed);

    LocalSearchSolver solver(
        problem, LocalSearchType::Steepest, MoveType::IntraEdgeSwap,
        StartSolutionType::Random, seed, seed);

    std::vector<int> current = solver.solve();
    int64_t current_cost = problem.FullDistanceAndCost(current);
    std::vector<int> best = current;
    int64_t best_cost = current_cost;
    int main_loop_iterations = 0;

    std::vector<int> last_destroyed; // <-- store previous destroyed solution
    int repeat_count = 0;            // <-- count repeats

    while (true)
    {
        double elapsed_ms =
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - t_start)
                .count();

        if (elapsed_ms >= timeLimit)
        {
            std::cout << "Destroyed solution repeated " << repeat_count
                      << " times during run.\n";
            return Result{best_cost, best, elapsed_ms, main_loop_iterations};
        }

        std::vector<int> destroyed = destroy(current, rng, heuristic);

        // std::cout << "Iteration " << main_loop_iterations + 1 << " destroyed nodes: ";
        // for (int node : destroyed)
        // {
        //     std::cout << node << " ";
        // }
        // std::cout << "\n";

        // Check if destroyed solution is identical to previous
        if (!last_destroyed.empty() && destroyed == last_destroyed)
        {
            repeat_count++;
            std::cout << "Warning: destroyed solution identical to previous iteration!\n";
        }

        last_destroyed = destroyed; // update last_destroyed

        std::vector<int> repaired = repair(destroyed);
        assert(!repaired.empty());

        std::vector<int> candidate;
        if (useLocalSearch)
        {
            solver.setGivenSolution(repaired);
            candidate = solver.solve_from_solution();
        }
        else
        {
            candidate = repaired;
        }

        int64_t cand_cost = problem.FullDistanceAndCost(candidate);

        if (cand_cost < best_cost)
        {
            best = candidate;
            best_cost = cand_cost;
            current = candidate;
            current_cost = cand_cost;
        }

        main_loop_iterations++;
        // std::cout << "best cost so far: " << best_cost << "\n";
    }
}

// -------------------- Repair --------------------
std::vector<int> DestroyRepairSolver::repair(const std::vector<int> &partial_solution)
{
    return greedy_solver_->complete_solution(partial_solution); // <-- access via unique_ptr
}

// -------------------- Destroy --------------------
std::vector<int> DestroyRepairSolver::destroy(const std::vector<int> &solution,
                                              std::mt19937 &rng,
                                              DestroyHeuristic heuristic)
{
    std::vector<int> destroyed = solution;
    int n = destroyed.size();
    int numToRemove = std::max(1, int(n * destroyFraction));

    std::vector<double> weights(n, 1.0);

    switch (heuristic)
    {
    case DestroyHeuristic::None:
        break; // all weights = 1
    case DestroyHeuristic::Node:
        for (int i = 0; i < n; i++)
        {
            int node = destroyed[i];
            weights[i] = problem.GetCost(node);
        }
        break;
    case DestroyHeuristic::Mixed:
        for (int i = 0; i < n; i++)
        {
            int node = destroyed[i];
            int prev_node = destroyed[(i - 1 + n) % n];
            int next_node = destroyed[(i + 1) % n];
            weights[i] = problem.GetCost(node) +
                         problem.GetDistance(prev_node, node) +
                         problem.GetDistance(node, next_node);
        }
        break;
    }

    for (int r = 0; r < numToRemove && !destroyed.empty(); ++r)
    {
        std::discrete_distribution<> dist(weights.begin(), weights.end());
        int idx = dist(rng);

        destroyed.erase(destroyed.begin() + idx);
        weights.erase(weights.begin() + idx);
    }

    return destroyed;
}

// -------------------- Initialize GreedySolver --------------------
void DestroyRepairSolver::initializeGreedySolver(int start_index, GreedyMode mode, Heuristic heuristic, float weight)
{
    greedy_solver_ = std::make_unique<GreedySolver>(problem, start_index, mode, heuristic, weight); // <-- dynamically
}

// std::vector<double>  DestroyRepairSolver::subroutineEdgeCalculation(std::vector<int> partial_solution){
//     int n = partial_solution.size();
//     std::vector<double> weights(n, 1.0);
//     for (int i = 0; i < n; i++)
//         {
//             int node = partial_solution[i];
//             int prev_node = partial_solution[(i - 1 + n) % n];
//             int next_node = partial_solution[(i + 1) % n];
//             weights[i] = problem.GetCost(node) +
//                          problem.GetDistance(prev_node, node) +
//                          problem.GetDistance(node, next_node);
//         }
//         std::reverse(std::sort(weights.begin(), weights.end()));
//         return weights[0];
//     }
//     }
