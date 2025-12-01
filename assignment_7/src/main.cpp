#include <iostream>
#include <filesystem>
#include "../../assignment_6/src/MlsSolver.h"
#include "../../src/ProblemInstance.h"
#include "../../assignment_6/src/Summary.h"
#include "./DestroyRepairSolver.h"
#include <fstream> // For file handling
#include <vector>

// NOTE: adjust the include paths if needed

int main()
{

    ProblemInstance prob1("./data/TSPA.csv", 100, "A");
    ProblemInstance prob2("./data/TSPB.csv", 100, "B");
    std::vector<ProblemInstance> probs = {prob1, prob2};

    std::string resultPath = "./assignment_6/results/";
    std::filesystem::create_directories(resultPath);

    const int MSLS_RESTARTS = 200;
    const int REPEATS = 5;

    for (auto &prob : probs)
    {

        std::cout << "\n=== Instance " << prob.name << " ===\n";

        // ---------------------
        // Run MSLS
        // ---------------------
        MlsSolver msls_runner(
            prob,
            MSLS_RESTARTS,
            /*seed=*/1000,
            MoveType::IntraEdgeSwap);

        auto msls_results = msls_runner.run_multiple(REPEATS);

        // Compute summary
        Summary ms = compute_summary(msls_results);
        std::cout << "MSLS SUMMARY:\n"
                  << "   Best cost:  " << ms.best_cost << "\n"
                  << "   Worst cost: " << ms.worst_cost << "\n"
                  << "   Avg cost:   " << ms.avg_cost << "\n"
                  << "   Best time:  " << ms.best_time << " ms\n"
                  << "   Worst time: " << ms.worst_time << " ms\n"
                  << "   Avg time:   " << ms.avg_time << " ms\n";

        // Save best solution
        {
            std::string f = resultPath + prob.name + "_MSLS_best.csv";
            Solver::writePathCsv(ms.best_solution, f);
        }

        // ---------------------
        // Run ILS (time-limited)
        // ---------------------
        for (bool useLocalSearch : {true, false})
        {
            for (DestroyHeuristic heuristic : {DestroyHeuristic::None, DestroyHeuristic::Node, DestroyHeuristic::Mixed})
            {
                std::cout << "\n--- Running Destroy-Repair Solver (useLocalSearch="
                          << (useLocalSearch ? "true" : "false") << ", heuristic="
                          << (heuristic == DestroyHeuristic::Node ? "Node" : heuristic == DestroyHeuristic::Mixed ? "Mixed"
                                                                                                                  : "None")
                          << ") ---\n";
                DestroyRepairSolver drs_runner(
                    prob,
                    useLocalSearch,
                    0.3,
                    heuristic,
                    ms.avg_time);

                auto drs_results = drs_runner.run_multiple(REPEATS);
                Summary is = compute_summary(drs_results);
                std::cout << "\nDRS SUMMARY:\n"
                          << "Local_Search=" << (useLocalSearch ? "true" : "false")
                          << ", heuristic=" << (heuristic == DestroyHeuristic::Node ? "Node" : heuristic == DestroyHeuristic::Mixed ? "Mixed"
                                                                                                                                    : "None")
                          << "\n"
                          << "   Best cost:  " << is.best_cost << "\n"
                          << "   Worst cost: " << is.worst_cost << "\n"
                          << "   Avg cost:   " << is.avg_cost << "\n"
                          << "   Best time:  " << is.best_time << " ms\n"
                          << "   Worst time: " << is.worst_time << " ms\n"
                          << "   Avg time:   " << is.avg_time << " ms\n";

                {
                    std::string f = resultPath + prob.name + "Local_Search_" + (useLocalSearch ? "true" : "false") + (heuristic == DestroyHeuristic::Node ? "Node" : heuristic == DestroyHeuristic::Mixed ? "Mixed"
                                                                                                                                                                                                          : "None") +
                                    "_drs_best.csv";
                    Solver::writePathCsv(is.best_solution, f);
                }
            }
        }
    }
    return 0;
}