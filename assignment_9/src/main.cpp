#include "../../src/ProblemInstance.h"
#include "../../src/Solver.h"
#include "EvolutionarySolver.h"
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

    std::string resultPath = "./assignment_9/results/";

    // --- Loop through configurations ---
    for (auto& prob : probs) {
        std::cout << "\n=== Testing " << prob.name << " ===\n";

        int numCities = prob.getNumCities();
        std::vector<int> bestSolution;

        // --- Iterate over all starting cities ---
        auto start = std::chrono::high_resolution_clock::now();

        EvolutionarySolver solver(prob, 30, true);
        std::vector<int> solution = solver.solve();
        int64_t cost = prob.FullDistanceAndCost(solution);

        auto end = std::chrono::high_resolution_clock::now();
        int duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        std::string filename = resultPath + "LS_best_" + prob.name + ".csv";

        if (!solver.writePathCsv(bestSolution, filename)) std::cerr << "Failed to write " << filename << std::endl;

        // --- Print results ---
        std::cout << "  Best:    " << cost << "\n";

        std::cout << "\n" << "  TIME (ms): " << duration << "\n";
    }

    return 0;
}
