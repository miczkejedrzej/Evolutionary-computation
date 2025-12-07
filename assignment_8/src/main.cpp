#include "../../src/ProblemInstance.h"
#include "../../src/Solver.h"
#include "../../assignment_3/src/LocalSearchSolver.h"
#include <unordered_set>
#include <iostream>
#include <climits>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <sstream>
#include <fstream>
#include <filesystem>

struct PairHash {
    size_t operator()(const std::pair<int,int>& p) const noexcept {
        return (size_t(uint32_t(p.first)) << 32) ^ uint32_t(p.second);
    }
};

void append_or_create(const std::string& path, int a, double b)
{
    // std::ios::app → open for append (creates the file if it doesn't exist)
    std::ofstream file(path, std::ios::app);
    file << std::fixed;

    // Write the line "a;b\n"
    file << (a * 100) << ';' << (b * 100) << '\n';
}

// 2 similarity measures
double NodeSimilarity(std::unordered_set<int>& sol1, std::unordered_set<int>& sol2) {
    int count = 0;
    for (int x : sol1) {
        if (sol2.count(x)) count++;
    }

    return double(count);
}

double EdgeSimilarity(std::unordered_set<std::pair<int, int>, PairHash>& sol1, std::unordered_set<std::pair<int, int>, PairHash>& sol2) {
    int count = 0;
    for (std::pair<int, int> x : sol1) {
        if (sol2.count(x)) count++;
    }

    return double(count);
}

int main() {
    int num_tests = 1000;
    int n = 100;

    using Edge = std::pair<int, int>;

    // --- Load problem instances ---
    ProblemInstance prob1("./data/TSPA.csv", n, "A");
    ProblemInstance prob2("./data/TSPB.csv", n, "B");
    std::vector<ProblemInstance> probs = {prob1, prob2};

    std::string resultPath = "./assignment_8/results/";

    for (const auto& entry : std::filesystem::directory_iterator(resultPath)) {
        if (entry.path().extension() == ".csv")
            std::filesystem::remove(entry.path());
    }

    // --- Loop through configurations ---
    for (ProblemInstance& prob : probs) {
        std::cout << "\n=== Testing " << prob.name << " ===\n";

        // int numCities = prob.getNumCities();
        int64_t best_cost = INT64_MAX;
        int best_i = -1;

        std::vector<std::unordered_set<int>> allSolutionsNodes;
        allSolutionsNodes.reserve(num_tests);
        std::vector<std::unordered_set<std::pair<int, int>, PairHash>> allSolutionsEdges;
        allSolutionsEdges.reserve(num_tests);
        std::vector<int> allSolutionsCosts;
        allSolutionsCosts.reserve(num_tests);

        for (int i = 0; i < num_tests; i++) {
            // 1000 local optima generated from random solutions using greedy local search
            LocalSearchSolver solver(prob, LocalSearchType::GreedyRandom, MoveType::IntraEdgeSwap, StartSolutionType::Random, i);
            std::vector <int> solution = solver.solve();
            int64_t cost = prob.FullDistanceAndCost(solution);
            
            // Add solution to allSolutions
            std::unordered_set<int> solutionSet(solution.begin(), solution.end());
            solutionSet.reserve(n);
            allSolutionsNodes.push_back(solutionSet);

            std::unordered_set<std::pair<int, int>, PairHash> edges;
            edges.reserve(solution.size());
            for (int j = 0; j < solution.size(); j++) {
                int a = solution[j];
                int b = solution[(j + 1) % solution.size()];
                edges.insert(std::make_pair(std::min(a, b), std::max(a, b)));
            }

            allSolutionsEdges.push_back(edges);

            allSolutionsCosts.push_back(cost);

            if (cost < best_cost) {
                best_cost = cost;
                best_i = i;
            }
        }

        // Load best solutions from previous assignments
        std::vector<int> bestPrev;
        std::unordered_set<int> bestPrevNodes;
        bestPrevNodes.reserve(n);
        std::unordered_set<std::pair<int, int>, PairHash> bestPrevEdges;

        std::ifstream file("./assignment_7/results/" + prob.name + "Local_Search_trueMixed_drs_best.csv");
        std::string line;
        if (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string token;

            while (std::getline(ss, token, ';')) {
                if (!token.empty())
                    bestPrev.push_back(std::stoi(token));
            }

            bestPrevNodes.insert(bestPrev.begin(), bestPrev.end());

            bestPrevEdges.reserve(n);
            for (size_t i = 0; i < n; ++i) {
                int a = bestPrev[i];
                int b = bestPrev[(i + 1) % n];
                bestPrevEdges.insert(std::make_pair(std::min(a, b), std::max(a, b)));
            }
        }

        int progress = 1;
        for (int i = 0; i < num_tests; i++) {
            // Similarity for all other (avg)
            double sumSimilarityEdge = 0;
            double sumSimilarityNode = 0;
            for (int j = 0; j < num_tests; j++) {
                if (j == i)
                    continue;
                
                sumSimilarityEdge += EdgeSimilarity(allSolutionsEdges[i], allSolutionsEdges[j]) / n;
                sumSimilarityNode += NodeSimilarity(allSolutionsNodes[i], allSolutionsNodes[j]) / n;
            }

            append_or_create(resultPath + prob.name + "_edge_avgSim.csv", allSolutionsCosts[i], sumSimilarityEdge / (num_tests - 1));
            append_or_create(resultPath + prob.name + "_node_avgSim.csv", allSolutionsCosts[i], sumSimilarityNode / (num_tests - 1));

            // Similarity to best found
            if (i != best_i) {
                append_or_create(resultPath + prob.name + "_edge_bestSim.csv", allSolutionsCosts[i], EdgeSimilarity(allSolutionsEdges[i], allSolutionsEdges[best_i]) / n);
                append_or_create(resultPath + prob.name + "_node_bestSim.csv", allSolutionsCosts[i], NodeSimilarity(allSolutionsNodes[i], allSolutionsNodes[best_i]) / n);
            }
            
            // Similarity to the best method
            append_or_create(resultPath + prob.name + "_edge_prevSim.csv", allSolutionsCosts[i], EdgeSimilarity(allSolutionsEdges[i], bestPrevEdges) / n);
            append_or_create(resultPath + prob.name + "_node_prevSim.csv", allSolutionsCosts[i], NodeSimilarity(allSolutionsNodes[i], bestPrevNodes) / n);
        }
    }

    return 0;
}
