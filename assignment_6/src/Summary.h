#include "MlsSolver.h"
#include <limits>
#include <iostream>
struct Summary {
    int64_t best_cost;
    int64_t worst_cost;
    double avg_cost;

    double best_time;
    double worst_time;
    double avg_time;

    std::vector<int> best_solution;
};

Summary compute_summary(const std::vector<Result>& results)
{
    Summary S;
    S.best_cost  = INT64_MAX;
    S.worst_cost = INT64_MIN;
    S.avg_cost   = 0;

    S.best_time  = 1e300;
    S.worst_time = -1e300;
    S.avg_time   = 0;

    for (const auto& r : results) {
        // cost
        S.avg_cost += r.best_cost;
        if (r.best_cost < S.best_cost) {
            S.best_cost = r.best_cost;
            S.best_solution = r.best_solution;
        }
        if (r.best_cost > S.worst_cost)
            S.worst_cost = r.best_cost;

        // time
        S.avg_time += r.elapsed_ms;
        if (r.elapsed_ms < S.best_time)  S.best_time = r.elapsed_ms;
        if (r.elapsed_ms > S.worst_time) S.worst_time = r.elapsed_ms;
    }

    int n = results.size();
    S.avg_cost /= n;
    S.avg_time /= n;

    return S;
}
