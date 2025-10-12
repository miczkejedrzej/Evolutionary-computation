#include "ProblemInstance.h"
#include "GreedySolver.h"
#include "RandomSolver.h"
#include <iostream>

int main() {
    ProblemInstance prob1("./TSPA.csv",100, "A");
    ProblemInstance prob2("./TSPB.csv",100, "B");
    //inside the solver result printed
    std::vector<ProblemInstance> probs = {prob1,prob2};

    for(ProblemInstance prob: probs){
        int64_t bestRandom = INT64_MAX;
        int64_t bestNNE = INT64_MAX;
        int64_t bestNN = INT64_MAX;
        int64_t bestGC = INT64_MAX;

        int64_t medRandom = 0;
        int64_t medNNE = 0;
        int64_t medNN = 0;
        int64_t medGC = 0;

        int64_t worseRandom = 0;
        int64_t worseNNE = 0;
        int64_t worseNN = 0;
        int64_t worseGC = 0;

        for (int i = 0; i < prob.getNumCities(); ++i) {
            RandomSolver solver(prob);
            std::vector<int> result = solver.solve();
            medRandom += prob.FullDistanceAndCost(result);
            if(prob.FullDistanceAndCost(result) < bestRandom){
                bestRandom = prob.FullDistanceAndCost(result);

                std::string filename = "Random_best_" + prob.name + ".csv";
                if (!solver.writePathCsv(result, filename)) {
                    std::cerr << "Failed to write to a file " + filename << std::endl;
                }
            }
            if(prob.FullDistanceAndCost(result) > worseRandom){
                worseRandom = prob.FullDistanceAndCost(result);
            }

            GreedySolver greedySolver(prob,i,GreedyMode::NearestNeighbourEnd);
            std::vector<int> greedyResult = greedySolver.solve();
            // std::cout << "Greedy Result (NearestNeighbourEnd): " << prob.FullDistanceAndCost(greedyResult) << std::endl;
            medNNE += prob.FullDistanceAndCost(greedyResult);
            if(prob.FullDistanceAndCost(greedyResult) < bestNNE){
                bestNNE = prob.FullDistanceAndCost(greedyResult);

                std::string filename = "NNE_best_" + prob.name + ".csv";
                if (!greedySolver.writePathCsv(greedyResult, filename)) {
                    std::cerr << "Failed to write to a file " + filename << std::endl;
                }
            }
            if(prob.FullDistanceAndCost(greedyResult) > worseNNE){
                worseNNE = prob.FullDistanceAndCost(greedyResult);
            }

            GreedySolver greedySolver2(prob,i,GreedyMode::NearestNeighbour);
            greedyResult = greedySolver2.solve();
            // std::cout << "Greedy Result: (NearestNeighbour): " << prob.FullDistanceAndCost(greedyResult) << std::endl;
            medNN += prob.FullDistanceAndCost(greedyResult);
            if(prob.FullDistanceAndCost(greedyResult) < bestNN){
                bestNN = prob.FullDistanceAndCost(greedyResult);

                std::string filename = "NN_best_" + prob.name + ".csv";
                if (!greedySolver.writePathCsv(greedyResult, filename)) {
                    std::cerr << "Failed to write to a file " + filename << std::endl;
                }
            }
            if(prob.FullDistanceAndCost(greedyResult) > worseNN){
                worseNN = prob.FullDistanceAndCost(greedyResult);
            }

            GreedySolver greedySolver3(prob,i,GreedyMode::GreedyCycle);
            greedyResult = greedySolver3.solve();
            // std::cout << "Greedy Result: (Greedy Cycle): " << prob.FullDistanceAndCost(greedyResult) << std::endl;
            medGC += prob.FullDistanceAndCost(greedyResult);
            if(prob.FullDistanceAndCost(greedyResult) < bestGC){
                bestGC = prob.FullDistanceAndCost(greedyResult);

                std::string filename = "GC_best_" + prob.name + ".csv";
                if (!greedySolver.writePathCsv(greedyResult, filename)) {
                    std::cerr << "Failed to write to a file " + filename << std::endl;
                }
            }
            if(prob.FullDistanceAndCost(greedyResult) > worseGC){
                worseGC = prob.FullDistanceAndCost(greedyResult);
            }
        }

        std::cout << "Best Random Result for TSP" + prob.name + ": " << bestRandom << std::endl;
        std::cout << "Worst Random Result for TSP" + prob.name + ": " << worseRandom << std::endl;
        std::cout << "Average Random Result for TSP" + prob.name + ": " << medRandom/prob.getNumCities() << std::endl << std::endl;
        
        std::cout << "Best Nearest Neighbour End Result for TSP" + prob.name + ": " << bestNNE << std::endl;
        std::cout << "Worst Nearest Neighbour End Result for TSP" + prob.name + ": " << worseNNE << std::endl;
        std::cout << "Average Nearest Neighbour End Result for TSP" + prob.name + ": " << medNNE/prob.getNumCities() << std::endl << std::endl;
        
        std::cout << "Best Nearest Neighbour Result for TSP" + prob.name + ": " << bestNN << std::endl;
        std::cout << "Worst Nearest Neighbour Result for TSP" + prob.name + ": " << worseNN << std::endl;
        std::cout << "Average Nearest Neighbour Result for TSP" + prob.name + ": " << medNN/prob.getNumCities() << std::endl << std::endl;

        std::cout << "Best Greedy Cycle Result for TSP" + prob.name + ": " << bestGC << std::endl;
        std::cout << "Worst Greedy Cycle Result for TSP" + prob.name + ": " << worseGC << std::endl;
        std::cout << "Average Greedy Cycle Result for TSP" + prob.name + ": " << medGC/prob.getNumCities() << std::endl << std::endl;
    };
    return 0;
}