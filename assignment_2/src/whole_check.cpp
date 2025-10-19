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
    ProblemInstance prob1("./data/TSPA.csv", 200, "A");
    ProblemInstance prob2("./data/TSPB.csv", 200, "B");
    std::vector<ProblemInstance> probs = {prob1, prob2};

    std::vector<float> weights = {0,0.5,1.0};

    for (auto& prob : probs) {
        // Loop over the two heuristics
        for (int heuristicInt = 0; heuristicInt <= 1; ++heuristicInt) {
            Heuristic heuristic = static_cast<Heuristic>(heuristicInt);
            std::string heuristicName = (heuristic == Heuristic::Regret) ? "REGRET" : "HYBRID";

            std::cout << "=== Testing " << prob.name << " with heuristic: " << heuristicName << " ===\n";

            if(heuristic == Heuristic::HybridRegretObjective){

                for(int i =0; i< weights.size();i++){

                    float weight_element = weights[i];
                    // Initialize best, worst, and median trackers
                    int64_t bestNN = INT64_MAX, bestGC = INT64_MAX;
                    int64_t worstNN = INT64_MIN, worstGC = INT64_MIN;
                    int64_t sumNN = 0, sumGC = 0;

                    int numCities = prob.getNumCities();

                    // --- iterate over all starting cities ---
                    for (int startIdx = 0; startIdx < numCities; ++startIdx) {
                        // --- Nearest Neighbour ---
                        std::string string_weight = std::to_string(weight_element);
                        GreedySolver solverNN(prob, startIdx, GreedyMode::NearestNeighbour, heuristic,weight_element);
                        solverNN.setWeight(weight_element);
                        std::vector<int> greedyNN = solverNN.solve();
                        int64_t nnCost = prob.FullDistanceAndCost(greedyNN);
                        sumNN += nnCost;
                        

                        if (nnCost < bestNN) {
                            bestNN = nnCost;
                            
                        }
                        if (nnCost > worstNN) worstNN = nnCost;

                        // --- Greedy Cycle ---
                        GreedySolver solverGC(prob, startIdx, GreedyMode::GreedyCycle, heuristic);
                        solverGC.setWeight(weight_element);
                        std::vector<int> greedyGC = solverGC.solve();
                        int64_t gcCost = prob.FullDistanceAndCost(greedyGC);
                        sumGC += gcCost;

                        if (gcCost < bestGC) {
                            bestGC = gcCost;
                
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
                    // std::cout << "  Worst:  " << worstNN << "\n";
                    // std::cout << "  Average: " << avgNN << "\n\n";

                    std::cout << "Greedy Cycle:\n";
                    std::cout << "  Best:   " << bestGC << "\n";
                    // std::cout << "  Worst:  " << worstGC << "\n";
                    // std::cout << "  Average: " << avgGC << "\n\n";
                }
        }

        }
    }
     return 0;
}
