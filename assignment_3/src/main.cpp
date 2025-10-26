#include "../../src/ProblemInstance.h"
#include "../../assignment_2/src/GreedySolver.h"
#include "../../src/Solver.h"
#include "LocalSearchSolver.h"
#include <iostream>
#include <climits>
#include <string>
#include <vector>
#include <algorithm>


int main() {
    // --- Load problem instances ---
    ProblemInstance prob1("./data/TSPA.csv", 100, "A");
    ProblemInstance prob2("./data/TSPB.csv", 100, "B");
    std::vector<ProblemInstance> probs = {prob1, prob2};

    std::string resultPath = "./assignment_3/results/";

    // --- Loop through configurations ---
    for (auto& prob : probs) {
        for (int startSolution = 0; startSolution <= 1; ++startSolution) {
            StartSolutionType startSolutionType = static_cast<StartSolutionType>(startSolution);
            std::string startSolutionName = (startSolutionType == StartSolutionType::Greedy) ? "GREEDY" : "RANDOM";

            for (int searchType = 0; searchType <= 1; ++searchType) {
                LocalSearchType localSearchType = static_cast<LocalSearchType>(searchType);
                std::string searchTypeName = (localSearchType == LocalSearchType::GreedyRandom) ? "GREEDY" : "STEEPEST";

                // Intra-move types: skip 0 (InterRoute) if not needed
                for (int intraMove = 1; intraMove <= 2; ++intraMove) {
                    MoveType intraMoveType = static_cast<MoveType>(intraMove);
                    std::string intraMoveName = (intraMoveType == MoveType::IntraNodeSwap) ? "NODE_SWAP" : "EDGE_SWAP";

                    std::cout << "\n=== Testing " << prob.name
                              << " | Start: " << startSolutionName
                              << " | Search: " << searchTypeName
                              << " | Move: " << intraMoveName << " ===\n";

                    // Tracking metrics
                    int64_t bestCost = INT64_MAX;
                    int64_t worstCost = INT64_MIN;
                    int64_t sumCost = 0;

                    int numCities = prob.getNumCities();
                    std::vector<int> bestSolution;

                    // --- Iterate over all starting cities ---
                    for (int startIdx = 0; startIdx < numCities; ++startIdx) {
                        LocalSearchSolver solver(prob, localSearchType, intraMoveType, startSolutionType, startIdx, startIdx);
                        std::vector<int> solution = solver.solve();
                        int64_t cost = prob.FullDistanceAndCost(solution);
                        sumCost += cost;

                        if (cost < bestCost) {
                            bestCost = cost;
                            bestSolution = solution;

                            std::string filename = resultPath + "LS_best_" +
                                prob.name + "_" +
                                startSolutionName + "_" +
                                searchTypeName + "_" +
                                intraMoveName + ".csv";

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
                }
            }
        }
    }

    return 0;
}
