#include "../../src/ProblemInstance.h"
#include "../../src/Solver.h"
#include "LocalSearchSolver.h"
#include <iostream>
#include <climits>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>

int main() {
    // --- Load problem instances ---
    ProblemInstance prob1("./data/TSPA.csv", 100, "A");
    ProblemInstance prob2("./data/TSPB.csv", 100, "B");
    std::vector<ProblemInstance> probs = {prob1, prob2};

    std::string resultPath = "./assignment_5_alt/results/";

    // --- Loop through configurations ---
    for (auto& prob : probs) {
        std::cout << "\n=== Testing " << prob.name << " ===\n";

        // Tracking metrics
        int64_t bestCost = INT64_MAX;
        int64_t worstCost = INT64_MIN;
        int64_t sumCost = 0;

        // Timing variables
        int bestTime = INT_MAX;
        int worstTime = INT_MIN;
        int64_t sumTime = 0;

        int numCities = prob.getNumCities();
        std::vector<int> bestSolution;

        // --- Iterate over all starting cities ---
        for (int startIdx = 0; startIdx < numCities; ++startIdx) {
            auto start = std::chrono::high_resolution_clock::now();

            LocalSearchSolver solver(prob, startIdx);
            std::vector<int> solution = solver.solve();
            int64_t cost = prob.FullDistanceAndCost(solution);
            sumCost += cost;

            auto end = std::chrono::high_resolution_clock::now();
            int duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            sumTime += duration;

            if (duration < bestTime)
                bestTime = duration;
            
            if (duration > worstTime)
                worstTime = duration;

            if (cost < bestCost) {
                bestCost = cost;
                bestSolution = solution;

                std::string filename = resultPath + "LS_best_" +
                    prob.name + ".csv";

                if (!solver.writePathCsv(bestSolution, filename))
                    std::cerr << "Failed to write " << filename << std::endl;
            }

            if (cost > worstCost)
                worstCost = cost;
        }

        double avgCost = static_cast<double>(sumCost) / numCities;

        // --- Print results ---
        std::cout << "  Best:    " << bestCost << "\n";
        std::cout << "  Worst:   " << worstCost << "\n";
        std::cout << "  Average: " << avgCost << "\n";

        std::cout << "\n" << "  TIME (ms)" << "\n"
                    << "  Best: " << bestTime << "\n"
                    << "  Worst: " << worstTime << "\n"
                    << "  Average: " << static_cast<double>(sumTime) / numCities << "\n";
    }

    return 0;
}
