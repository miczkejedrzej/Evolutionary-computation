#include "IlsSolver.h"
#include <algorithm>
#include <numeric>
#include <random>
#include <sstream>
#include <iostream>

IlsSolver::IlsSolver(const ProblemInstance& prob, double time_budget_ms, int seed, MoveType intra, int perturbation_strength, int gene_pool_extension_size, bool use_double_bridge)
    : prob_(prob), time_budget_ms_(time_budget_ms), seed_(seed), intra_(intra), perturbation_strength_(perturbation_strength), gene_pool_extension_size_(gene_pool_extension_size), use_double_bridge_(use_double_bridge) {}





std::vector<int> IlsSolver::perturb(const std::vector<int>& tour,const ProblemInstance& problem, std::mt19937& rng,int number_perturbations) {
   
    std::vector<int> tour_copy = tour;
    std::vector<int> unselected = problem.GiveIndices();

    // Remove selected from unselected
    for (int idx : tour)
        unselected.erase(std::remove(unselected.begin(), unselected.end(), idx), unselected.end());
    for (int pert = 0; pert < number_perturbations; ++pert) {
        int operation = rng() % 3;
        int i = rng() % tour_copy.size();
        int j = rng() % unselected.size();
        if (i == j){operation = 0;}
        else if ( j == i+1 || (i == tour_copy.size() -1 && j ==0)) {operation =operation %2;}

        if (operation == 0)
            std::swap(tour_copy[i], unselected[j]);
        else if (operation == 1)
            std::swap(tour_copy[i], tour_copy[j]);
        else{
            if(i>j){std::swap(i,j);} 
            std::reverse(tour_copy.begin() + i + 1, tour_copy.begin() + j + 1);
        }// else
        // tour_copy = double_bridge(tour_copy, rng);
    }
    return tour_copy;
}

std::vector<int> IlsSolver::extend_gene_pool(
        const std::vector<int>& tour,
        const ProblemInstance& problem,
        int exchange_size)
{
    std::vector<int> selected = tour;
    std::vector<int> unselected = problem.GiveIndices();

    // Remove cities already in the tour
    for (int idx : tour) {
        unselected.erase(std::remove(unselected.begin(),
                                     unselected.end(),
                                     idx),
                         unselected.end());
    }

    std::mt19937 rng(std::random_device{}());

    // Cannot exchange more than min(|tour|, |unselected|)
    int k = std::min<int>(exchange_size,
                          std::min(selected.size(), unselected.size()));

    if (k <= 0)
        return selected;   // nothing to exchange

    // Choose k random indices in the tour
    std::vector<int> tour_indices(selected.size());
    std::iota(tour_indices.begin(), tour_indices.end(), 0);
    std::shuffle(tour_indices.begin(), tour_indices.end(), rng);

    // Choose k random cities from unselected set
    std::shuffle(unselected.begin(), unselected.end(), rng);

    // Perform the exchange
    for (int i = 0; i < k; i++) {
        int tour_pos = tour_indices[i];
        int city_in = unselected[i];        // add to tour
        int city_out = selected[tour_pos];  // remove from tour

        selected[tour_pos] = city_in;       // replace
        // (city_out is now implicitly removed from the tour)
    }

    return selected;
}

std::vector<int> IlsSolver::double_bridge(const std::vector<int>& tour, std::mt19937& rng) {
    int n = (int)tour.size();
    std::uniform_int_distribution<int> cut_dist(1, n-1);
    int a = cut_dist(rng);
    int b = cut_dist(rng);
    int c = cut_dist(rng);
    int d = cut_dist(rng);
    // ensure sorted and distinct
    std::vector<int> cuts = {a,b,c,d};
    std::sort(cuts.begin(), cuts.end());
    a = cuts[0]; b = cuts[1]; c = cuts[2]; d = cuts[3];
    // segments s1:[0,a), s2:[a,b), s3:[b,c), s4:[c,d), s5:[d,n)
    std::vector<int> s1(tour.begin(), tour.begin() + a);
    std::vector<int> s2(tour.begin() + a, tour.begin() + b);
    std::vector<int> s3(tour.begin() + b, tour.begin() + c);
    std::vector<int> s4(tour.begin() + c, tour.begin() + d);
    std::vector<int> s5(tour.begin() + d, tour.end());
    // common reconnection: s1 + s4 + s3 + s2 + s5
    std::vector<int> out;
    out.reserve(n);
    out.insert(out.end(), s1.begin(), s1.end());
    out.insert(out.end(), s4.begin(), s4.end());
    out.insert(out.end(), s3.begin(), s3.end());
    out.insert(out.end(), s2.begin(), s2.end());
    out.insert(out.end(), s5.begin(), s5.end());
    return out;
}

Result IlsSolver::run_once() {
    auto t_start = std::chrono::high_resolution_clock::now();
    std::mt19937 rng(seed_);
    int n = prob_.GetNumberCitiesInCycle();
    int ls_runs = 0;
    // initial solution random -> use LocalSearchSolver with Random start
    LocalSearchSolver solver(prob_, LocalSearchType::Steepest, intra_, StartSolutionType::Random, 0, seed_);
    std::vector<int> current = solver.solve();
    ++ls_runs;
    int64_t current_cost = prob_.FullDistanceAndCost(current);
    std::vector<int> best = current;
    int64_t best_cost = current_cost;
    std::unordered_set<std::string> visited;
    std::unordered_set<int> optima_reached;
    optima_reached.insert(current_cost);
    int how_many_times_optima_reached = 0;
    // iterate until time budget exhausted
    while (true) {
        auto now = std::chrono::high_resolution_clock::now();
        double elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - t_start).count();
        if (elapsed_ms >= time_budget_ms_) break;

        // perturb
        // std::vector<int> pert = double_bridge(current, rng);
        std::vector<int> pert = perturb(current, prob_, rng, perturbation_strength_);
        int ls_seed = seed_ + ls_runs + 1;
        //LocalSearchSolver ls_solver(prob_, LocalSearchType::Steepest, intra_, StartSolutionType::Given, 0, ls_seed);
        solver.setGivenSolution(pert);
        std::vector<int> cand = solver.solve_from_solution();
        ++ls_runs;
        int64_t cand_cost = prob_.FullDistanceAndCost(cand);
        if (optima_reached.find(cand_cost) != optima_reached.end()) {
            how_many_times_optima_reached +=1;
            // already seen optimum, scramble instead
            pert = extend_gene_pool(pert, prob_, gene_pool_extension_size_);
           //even stronger perturbation with the use of the double bridge.
            if (use_double_bridge_) {
                pert = double_bridge(pert, rng);
            }
            solver.setGivenSolution(pert);
            cand = solver.solve_from_solution();
            ++ls_runs;
            cand_cost = prob_.FullDistanceAndCost(cand);
            if(cand_cost < best_cost) {
                best = cand;
                best_cost = cand_cost;
            }
            current = cand;
            current_cost = cand_cost;
        } else {
        
            optima_reached.insert(cand_cost);
            current = cand;
            current_cost = cand_cost;
            if (cand_cost < best_cost) {
                best = cand;
                best_cost = cand_cost;
            }
            continue;
        } 
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed_ms_total = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();
    std::cout << "different [ILS] Optima reached " << optima_reached.size() << " times, of which " << how_many_times_optima_reached << " were repeats." << std::endl;
    return Result{best_cost, best, elapsed_ms_total, ls_runs};
}

std::vector<Result> IlsSolver::run_multiple(int runs) {
    std::vector<Result> results;
    results.reserve(runs);
    for (int i = 0; i < runs; ++i) {
        // vary seed per repetition
        seed_ += i;
        Result r = run_once();
        results.push_back(r);
        std::cout << "[ILS] run " << (i+1) << "/" << runs << " cost=" << r.best_cost
                  << " time(ms)=" << r.elapsed_ms << " ls_runs=" << r.ls_runs << std::endl;
    }
    return results;
}
