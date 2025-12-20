#pragma once
#include "../../src/ProblemInstance.h"
#include "../../src/Solver.h"
#include "../../assignment_2/src/GreedySolver.h"
#include "../../assignment_5/src/LocalSearchSolver.h"
#include <vector>
#include <random>
#include <algorithm>
#include <climits>
#include <memory>

enum class RecombinationType { CommonEdges, Repair };

class EvolutionarySolver : public Solver {
public:
    EvolutionarySolver(const ProblemInstance& prob, RecombinationType recombinationType, int randomSeed = 42, bool performLS = true,
        int popSize = 21, bool randReplace = false);

    std::vector<int> solve();

private:
    
    struct Solution {
        std::vector<int> genome;
        int fitness;

        bool operator==(const Solution& a) const {
            return fitness == a.fitness;
        }
    };

    std::mt19937 rng_;
    RecombinationType recombinationType;
    bool performLS;
    int popSize;
    bool randReplace;

    std::unique_ptr<GreedySolver> greedy_solver;

    std::vector<Solution> population;

    std::vector<int> slicing(const std::vector<int>& arr, int X);

    std::vector<int> initializeSolution();
    
    void AssertHamiltonian(const std::vector<int>& visited, int citiesNumber);

    bool IsInPopulation(const Solution& s) const;
    void SortPopulation();
    Solution CrossCommonEdges(const Solution& p1, const Solution& p2);
    Solution CrossRepair(const Solution& p1, const Solution& p2);

    LocalSearchSolver lss;
};