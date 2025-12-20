#include "EvolutionarySolver.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
#include <random>
#include <cassert>
#include <unordered_set>
#include <chrono>
#include <unordered_map>
#include "../../src/ProblemInstance.h"
EvolutionarySolver::EvolutionarySolver(const ProblemInstance &prob, RecombinationType recombinationType, int randomSeed, bool performLS, int popSize_, bool randReplace_)
    : Solver(prob),
    recombinationType(recombinationType),
    performLS(performLS),
    rng_(randomSeed),
    popSize(std::max(popSize_, 21)),
    randReplace(randReplace_ && popSize_ > 21),
    lss(prob, rng_())
{}

// ---------------- Helper Functions ----------------

std::vector<int> EvolutionarySolver::slicing(const std::vector<int> &arr, int X)
{
    // Ensure X is within bounds
    assert(X >= 0 && X < static_cast<int>(arr.size()));

    // Copy first X+1 elements
    std::vector<int> result(arr.begin(), arr.begin() + X + 1);
    return result;
}

void EvolutionarySolver::AssertHamiltonian(const std::vector<int> &visited, int citiesNumber)
{
    std::unordered_set<int> uniqueCities(visited.begin(), visited.end());
    assert((int)visited.size() == citiesNumber);
    assert((int)uniqueCities.size() == citiesNumber);
}

// ---------------- STARTING SOLUTION ----------------

std::vector<int> EvolutionarySolver::initializeSolution()
{
    std::vector<int> sol = problem.GiveIndices();
    std::shuffle(sol.begin(), sol.end(), rng_);
    int numCitiesInCycle = static_cast<int>(problem.GetNumberCitiesInCycle());
    std::vector<int> result = slicing(sol, numCitiesInCycle - 1);
    return result;
}

// ---------------- SOLUTION HELPERS ----------------

bool EvolutionarySolver::IsInPopulation(const Solution& s) const {
    for (const Solution& o : population) {
        if (s == o) return true;
    }
    return false;
}

void EvolutionarySolver::SortPopulation() {
    std::sort(population.begin(), population.end(),
    [](const Solution& a, const Solution& b) {
        return a.fitness < b.fitness;
    });
}

EvolutionarySolver::Solution EvolutionarySolver::CrossCommonEdges(const Solution& p1, const Solution& p2)
{
    using Path = std::vector<int>;

    // --- Build adjacency maps (undirected cycles, immutable) ---
    auto buildAdj = [](const std::vector<int>& c) {
        std::unordered_map<int, std::pair<int,int>> adj;
        int n = c.size();
        for (int i = 0; i < n; ++i) {
            int v = c[i];
            int prev = c[(i - 1 + n) % n];
            int next = c[(i + 1) % n];
            adj[v] = {prev, next};
        }
        return adj;
    };

    const auto adj1 = buildAdj(p1.genome);
    const auto adj2 = buildAdj(p2.genome);

    auto hasEdge = [](const auto& adj, int u, int v) {
        auto it = adj.find(u);
        if (it == adj.end()) return false;
        return it->second.first == v || it->second.second == v;
    };

    std::unordered_set<int> used;
    std::vector<Path> subpaths;

    // --- Extract maximal common subpaths (length >= 2) ---
    for (const auto& [v, _] : adj1) {
        if (used.count(v)) continue;
        if (!adj2.count(v)) continue;

        for (int n : {adj1.at(v).first, adj1.at(v).second}) {
            if (used.count(n)) continue;
            if (!hasEdge(adj2, v, n)) continue;

            // ----- Extend backward -----
            Path backward;
            int prev = v;
            int curr = n;

            while (!used.count(curr) &&
                   hasEdge(adj1, prev, curr) &&
                   hasEdge(adj2, prev, curr)) {

                backward.push_back(curr);

                int next1 = (adj1.at(curr).first == prev)
                          ? adj1.at(curr).second
                          : adj1.at(curr).first;
                int next2 = (adj2.at(curr).first == prev)
                          ? adj2.at(curr).second
                          : adj2.at(curr).first;

                if (next1 != next2)
                    break;

                prev = curr;
                curr = next1;
            }

            // ----- Extend forward (opposite direction) -----
            Path forward;
            prev = v;
            curr = (adj1.at(v).first == n)
                 ? adj1.at(v).second
                 : adj1.at(v).first;

            while (!used.count(curr) &&
                   hasEdge(adj1, prev, curr) &&
                   hasEdge(adj2, prev, curr)) {

                forward.push_back(curr);

                int next1 = (adj1.at(curr).first == prev)
                          ? adj1.at(curr).second
                          : adj1.at(curr).first;
                int next2 = (adj2.at(curr).first == prev)
                          ? adj2.at(curr).second
                          : adj2.at(curr).first;

                if (next1 != next2)
                    break;

                prev = curr;
                curr = next1;
            }

            if (!backward.empty() || !forward.empty()) {
                Path path;

                // backward is v → ...
                path.insert(path.end(), backward.rbegin(), backward.rend());
                path.push_back(v);
                path.insert(path.end(), forward.begin(), forward.end());

                if (path.size() >= 2) {
                    for (int x : path)
                        used.insert(x);

                    subpaths.push_back(path);
                }
            }
        }
    }

    // --- Add single common vertices ---
    for (const auto& [v, _] : adj1) {
        if (used.count(v)) continue;
        if (adj2.count(v)) {
            subpaths.push_back({v});
            used.insert(v);
        }
    }

    // --- Add arbitrary vertices until target length ---
    int targetLength = problem.GetNumberCitiesInCycle();
    int vertexCount  = problem.getNumCities();

    std::vector<int> allVertices(vertexCount);
    std::iota(allVertices.begin(), allVertices.end(), 0);
    std::shuffle(allVertices.begin(), allVertices.end(), rng_);

    for (int v : allVertices) {
        if ((int)used.size() == targetLength)
            break;

        if (!used.count(v)) {
            subpaths.push_back({v});
            used.insert(v);
        }
    }

    if ((int)used.size() != targetLength)
        throw std::logic_error("Unable to construct solution of required length");

    // --- Randomize order and orientation ---
    std::shuffle(subpaths.begin(), subpaths.end(), rng_);
    for (auto& p : subpaths)
        if (p.size() > 1 && (rng_() & 1))
            std::reverse(p.begin(), p.end());

    // --- Merge into final solution ---
    Solution offspring;
    offspring.genome.reserve(targetLength);
    for (const auto& p : subpaths)
        offspring.genome.insert(offspring.genome.end(), p.begin(), p.end());

    return offspring;
}

EvolutionarySolver::Solution EvolutionarySolver::CrossRepair(const Solution& p1, const Solution& p2) {
    std::vector<int> v1 = p1.genome;
    std::vector<int> v2 = p2.genome;
    std::unordered_set<int> allowed(v2.begin(), v2.end());

    v1.erase(
        std::remove_if(v1.begin(), v1.end(),
            [&](int x) { return allowed.find(x) == allowed.end(); }),
        v1.end()
    );

    // Repair
    greedy_solver->complete_solution(v1);

    Solution offspring;
    offspring.genome = v1;
    return offspring;
}

// ---------------- MAIN SOLVE LOOP ----------------

std::vector<int> EvolutionarySolver::solve()
{
    greedy_solver = std::make_unique<GreedySolver>(problem, rng_(), GreedyMode::NearestNeighbour, Heuristic::HybridRegretObjective, 0.5);

    // Initialize the set of solutions
    while (population.size() < popSize) {
        std::vector<int> newSol = initializeSolution();
        Solution s{newSol, static_cast<int>(problem.FullDistanceAndCost(newSol))};

        if (!IsInPopulation(s))
            population.push_back(s);
    }
    SortPopulation();

    // Repeat until 57000 ms pass
    auto start = std::chrono::high_resolution_clock::now();
    auto currTime = start;
    while (std::chrono::duration_cast<std::chrono::seconds>(currTime - start).count() < 1) {
        // Cross two random solutions in the population to create an offspring
        std::uniform_int_distribution<size_t> d1(0, popSize - 1);
        std::uniform_int_distribution<size_t> d2(0, popSize - 2);
        size_t i = d1(rng_);
        size_t j = d2(rng_);
        if (j >= i) ++j;
        Solution newSol = recombinationType == RecombinationType::CommonEdges ?
            CrossCommonEdges(population[i], population[j]) : CrossRepair(population[i], population[j]);
        
        if (performLS) {
            // Perform LocalSearch (from assignment 5) on it
            lss.SetStartingSol(newSol.genome);
            newSol.genome = lss.solve();
        }

        newSol.fitness = problem.FullDistanceAndCost(newSol.genome);
        // Add only if not yet in the solution
        if (IsInPopulation(newSol)) {
            currTime = std::chrono::high_resolution_clock::now();
            continue;
        }

        // Replace a solution from population
        if (!randReplace) {
            population.push_back(newSol);
            SortPopulation();
            population.pop_back();
        }
        else {
            std::uniform_int_distribution<size_t> d(20, popSize - 1);
            population[d(rng_)] = newSol;
            SortPopulation();
        }
        
        currTime = std::chrono::high_resolution_clock::now();
    }

    AssertHamiltonian(population[0].genome, problem.GetNumberCitiesInCycle());

    return population[0].genome;
}
