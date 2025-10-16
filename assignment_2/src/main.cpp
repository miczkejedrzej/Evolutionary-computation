#include "../../src/ProblemInstance.h"
#include "GreedySolver.h"
#include "../../src/Solver.h"
#include <iostream>
#include <climits>
#include <string>
#include <vector>
#include <algorithm>

int main() {
    // Load problem instances
    ProblemInstance prob1("./data/TSPA.csv", 100, "A");
    ProblemInstance prob2("./data/TSPB.csv", 100, "B");
    std::vector<ProblemInstance> probs = {prob1, prob2};

    std::string resultPath = "./assignment_2/results/";

    for (auto& prob : probs) {
        // Loop over the two heuristics
        for (int heuristicInt = 0; heuristicInt <= 1; ++heuristicInt) {
            Heuristic heuristic = static_cast<Heuristic>(heuristicInt);
            std::string heuristicName = (heuristic == Heuristic::Regret) ? "REGRET" : "HYBRID";

            std::cout << "=== Testing " << prob.name << " with heuristic: " << heuristicName << " ===\n";

            // Initialize best, worst, and median trackers
            int64_t bestNN = INT64_MAX, bestGC = INT64_MAX;
            int64_t worstNN = INT64_MIN, worstGC = INT64_MIN;
            int64_t sumNN = 0, sumGC = 0;

            int numCities = prob.getNumCities()/20;

            // --- iterate over all starting cities ---
            for (int startIdx = 0; startIdx < numCities; ++startIdx) {
                // --- Nearest Neighbour ---
                GreedySolver solverNN(prob, startIdx, GreedyMode::NearestNeighbour, heuristic);
                std::vector<int> greedyNN = solverNN.solve();
                int64_t nnCost = prob.FullDistanceAndCost(greedyNN);
                sumNN += nnCost;

                if (nnCost < bestNN) {
                    bestNN = nnCost;
                    std::string filename = resultPath + "NN_best_" + heuristicName + "_" + prob.name + ".csv";
                    if (!solverNN.writePathCsv(greedyNN, filename))
                        std::cerr << "Failed to write " << filename << std::endl;
                }
                if (nnCost > worstNN) worstNN = nnCost;

                // --- Greedy Cycle ---
                GreedySolver solverGC(prob, startIdx, GreedyMode::GreedyCycle, heuristic);
                std::vector<int> greedyGC = solverGC.solve();
                int64_t gcCost = prob.FullDistanceAndCost(greedyGC);
                sumGC += gcCost;

                if (gcCost < bestGC) {
                    bestGC = gcCost;
                    std::string filename = resultPath + "GC_best_" + heuristicName + "_" + prob.name + ".csv";
                    if (!solverGC.writePathCsv(greedyGC, filename))
                        std::cerr << "Failed to write " << filename << std::endl;
                }
                if (gcCost > worstGC) worstGC = gcCost;
            }

            // Compute average (median approximation)
            double avgNN = static_cast<double>(sumNN) / numCities;
            double avgGC = static_cast<double>(sumGC) / numCities;

            // --- Print results ---
            std::cout << "\n--- Results for " << prob.name << " (" << heuristicName << ") ---\n";
            std::cout << "Nearest Neighbour:\n";
            std::cout << "  Best:   " << bestNN << "\n";
            std::cout << "  Worst:  " << worstNN << "\n";
            std::cout << "  Average: " << avgNN << "\n\n";

            std::cout << "Greedy Cycle:\n";
            std::cout << "  Best:   " << bestGC << "\n";
            std::cout << "  Worst:  " << worstGC << "\n";
            std::cout << "  Average: " << avgGC << "\n\n";
        }
    }

    return 0;
}
