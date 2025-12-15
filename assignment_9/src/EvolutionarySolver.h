#pragma once
#include "../../src/ProblemInstance.h"
#include "../../src/Solver.h"
#include "../../assignment_5/src/LocalSearchSolver.h"
#include <vector>
#include <random>
#include <algorithm>
#include <climits>

class EvolutionarySolver : public Solver {
public:
    EvolutionarySolver(const ProblemInstance& prob, int popSize = 21, bool randReplace = false,
                      int randomSeed = 42);

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
    int popSize;
    bool randReplace;

    std::vector<Solution> population;

    std::vector<int> slicing(const std::vector<int>& arr, int X);

    std::vector<int> initializeSolution();
    
    void AssertHamiltonian(const std::vector<int>& visited, int citiesNumber);

    bool IsInPopulation(const Solution& s) const;
    void SortPopulation();
    Solution Cross(const Solution& p1, const Solution& p2);

    LocalSearchSolver lss;
};