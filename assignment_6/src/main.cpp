#include <iostream>
#include <filesystem>
#include "MlsSolver.h"
#include "IlsSolver.h"
#include "../../src/ProblemInstance.h"
#include "Summary.h"
#include <fstream> // For file handling

// NOTE: adjust the include paths if needed

int main(){
   

   
    ProblemInstance prob1("./data/TSPA.csv", 100, "A");
    ProblemInstance prob2("./data/TSPB.csv", 100, "B");
    std::vector<ProblemInstance> probs = {prob1, prob2};

    std::string resultPath = "./assignment_6/results/";
    std::filesystem::create_directories(resultPath);

    const int MSLS_RESTARTS = 200;
    const int REPEATS = 20;

    for (auto& prob : probs) {

        std::cout << "\n=== Instance " << prob.name << " ===\n";

        // ---------------------
        // Run MSLS
        // ---------------------
        MlsSolver msls_runner(
            prob,
            MSLS_RESTARTS,
            /*seed=*/1000,
            MoveType::IntraEdgeSwap
        );

        auto msls_results = msls_runner.run_multiple(REPEATS);

        // Compute summary
        Summary ms = compute_summary(msls_results);
        std::cout << "MSLS SUMMARY:\n"
                  << "   Best cost:  " << ms.best_cost  << "\n"
                  << "   Worst cost: " << ms.worst_cost << "\n"
                  << "   Avg cost:   " << ms.avg_cost   << "\n"
                  << "   Best time:  " << ms.best_time  << " ms\n"
                  << "   Worst time: " << ms.worst_time << " ms\n"
                  << "   Avg time:   " << ms.avg_time   << " ms\n";

        // Save best solution
        {
            std::string f = resultPath + prob.name + "_MSLS_best.csv";
            Solver::writePathCsv(ms.best_solution, f);
        }

        // ---------------------
        // Run ILS (time-limited)
        // ---------------------
        IlsSolver ils_runner(
            prob,
            /*time_limit_ms=*/ms.avg_time,
            /*seed=*/2000,
            MoveType::IntraEdgeSwap,
            10,
            10, 1
        );

        auto ils_results = ils_runner.run_multiple(REPEATS);

        Summary is = compute_summary(ils_results);
        std::cout << "\nILS SUMMARY:\n"
                  << "   Best cost:  " << is.best_cost  << "\n"
                  << "   Worst cost: " << is.worst_cost << "\n"
                  << "   Avg cost:   " << is.avg_cost   << "\n"
                  << "   Best time:  " << is.best_time  << " ms\n"
                  << "   Worst time: " << is.worst_time << " ms\n"
                  << "   Avg time:   " << is.avg_time   << " ms\n";
        {
            std::string f = resultPath + prob.name + "_ILS_best.csv";
            Solver::writePathCsv(is.best_solution, f);
        }
    }

    return 0;
}