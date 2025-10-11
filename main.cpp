#include "ProblemInstance.h"
#include "GreedySolver.h"
#include "RandomSolver.h"
#include <iostream>

int main() {
    ///the pipeline with random
    std::mt19937 g(156064);
    //random starting point for now 
    int startingIdx = g() % 200;
    ProblemInstance prob1("./TSPA.csv",100);
    ProblemInstance prob2("./TSPB.csv",100);
    //inside the solver result printed
    std::vector<ProblemInstance> probs = {prob1,prob2};
    for(ProblemInstance prob: probs){
        RandomSolver solver(prob);
        std::vector<int> result = solver.solve();

        

        GreedySolver greedySolver(prob,startingIdx,GreedyMode::NearestNeighbourEnd);
        std::vector<int> greedyResult = greedySolver.solve();
        std::cout << "Greedy Result (NearestNeighbourEnd): " << prob.FullDistanceAndCost(greedyResult) << std::endl;

        GreedySolver greedySolver2(prob,startingIdx,GreedyMode::NearestNeighbour);
        greedyResult = greedySolver2.solve();
        std::cout << "Greedy Result: (NearestNeighbour): " << prob.FullDistanceAndCost(greedyResult) << std::endl;


        GreedySolver greedySolver3(prob,startingIdx,GreedyMode::GreedyCycle);
        greedyResult = greedySolver3.solve();
        std::cout << "Greedy Result: (Greedy Cycle): " << prob.FullDistanceAndCost(greedyResult) << std::endl;

      
    };
    return 0;
}