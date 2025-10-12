#include <vector>
#include <iostream>
#include <random>
#include <list>
#include "GreedySolver.h"
#include "ProblemInstance.h"
#include "Solver.h"
#include <unordered_map>
#include <unordered_set>
#include <cassert>


GreedySolver::GreedySolver(const ProblemInstance& prob,int starting_index, GreedyMode mode)
    : Solver(prob), mode_(mode), starting_index(starting_index) {}
    
void GreedySolver::AssertHamiltonian(std::vector<int> visited,int citiesNumber){
    std::unordered_set<int> uniqueCities;
    for(int element: visited){
        uniqueCities.insert(element);
    }
    assert(uniqueCities.size() == citiesNumber);
}

std::vector<int> GreedySolver::solve() {
    switch (mode_) {
        case GreedyMode::NearestNeighbour:
            return solveNearestNeighbour();
        case GreedyMode::NearestNeighbourEnd:
            return solveNearestNeighbourEnd();
        case GreedyMode::GreedyCycle:
            return solveGreedyCycle();
        default:
            return {};
    }
}

int GreedySolver::getStartingIndex(){
    return starting_index;
}

std::vector<int> GreedySolver::solveNearestNeighbour() {
    std::cout << "Solving with Nearest Neighbour";

    int targetSize = problem.GetNumberCitiesInCycle();
    std::vector<int> visited;
    visited.reserve(targetSize);
    std::vector<int> unvisited = problem.GiveIndices();
    int starting_index = getStartingIndex();
    visited.push_back(unvisited[starting_index]);
    unvisited.erase(unvisited.begin() + starting_index);

    // add to the path until the proper size 
    while (visited.size() < targetSize) {
        int bestInsertPos = -1;
        int bestCityIdx = -1;
        int64_t bestCost = INT64_MAX;

        // Try inserting each unvisited city at each possible position
        for (int cityIdx = 0; cityIdx < unvisited.size(); ++cityIdx) {
            int city = unvisited[cityIdx];
            
            // Case for extending the beggining - no need for rerouting the previous path
            int64_t costBegin = problem.GetCostAndDistance(visited[0],city);
            if (costBegin < bestCost) {
                bestCost = costBegin;
                bestInsertPos = 0;
                bestCityIdx = cityIdx;
            }

            // Try inserting at end position, no need to rerouting 
            int64_t costEnd = problem.GetCostAndDistance(visited.back(), city);
            if (costEnd < bestCost) {
                bestCost = costEnd;
                bestInsertPos = static_cast<int>(visited.size());
                bestCityIdx = cityIdx;
            }

            // Try inserting in middle positions ( - account for the rerouting to keep the solution
            // a proper hamiltonian
            for (int pos = 1; pos < visited.size(); ++pos) {
                int64_t oldEdge = problem.GetCostAndDistance(visited[pos-1], visited[pos]);
                int64_t newEdges = problem.GetCostAndDistance(visited[pos-1], city) + 
                                    problem.GetCostAndDistance(city, visited[pos]);
                int64_t deltaCost = newEdges - oldEdge;

                if (deltaCost < bestCost) {
                    bestCost = deltaCost;
                    bestInsertPos = pos;
                    bestCityIdx = cityIdx;
                }
            }
        }
    

        if (bestCityIdx == -1) break;

        // Insert the best city at the best position
        visited.insert(visited.begin() + bestInsertPos, unvisited[bestCityIdx]);
        unvisited.erase(unvisited.begin() + bestCityIdx);
    }
    AssertHamiltonian(visited,targetSize);
    return visited;
}

std::vector<int> GreedySolver::solveNearestNeighbourEnd() {
    std::cout << "Solving with Nearest Neighbour from End...\n";
    std::vector<int> unvisited = problem.GiveIndices();
    std::vector<int> visited;
    visited.reserve(problem.GetNumberCitiesInCycle());
    int currIdx = getStartingIndex();
    visited.push_back(unvisited[currIdx]);
    unvisited.erase(unvisited.begin() + currIdx);
    int targetSize = problem.GetNumberCitiesInCycle();
    for (int i = 1; i < targetSize; ++i) {
        int nearestIdx = -1;
        int64_t nearestDist = INT64_MAX;
        for (int j = 0; j < unvisited.size(); ++j) {
            int64_t dist = problem.GetCostAndDistance(visited[i - 1], unvisited[j]);
            if (dist < nearestDist) {
                nearestDist = dist;
                nearestIdx = j;
            }
        }
        visited.push_back(unvisited[nearestIdx]);
        unvisited.erase(unvisited.begin() + nearestIdx);
    }
    AssertHamiltonian(visited,targetSize);
    return visited;
}

std::vector<int> GreedySolver::solveGreedyCycle() {
    std::cout << "Solving with Greedy Cycle";
    int targetSize = problem.GetNumberCitiesInCycle();
    std::vector<int> visited;
    visited.reserve(targetSize);
    std::vector<int> unvisited = problem.GiveIndices();
    int starting_index = getStartingIndex();
    visited.push_back(unvisited[starting_index]);
    unvisited.erase(unvisited.begin() + starting_index);

    // Continue until targeted size achieved 
    while (visited.size() < targetSize) {
        int bestInsertPos = -1;
        int bestCityIdx = -1;
        int64_t bestCost = INT64_MAX;

        //case for the first edge - no cycle closing calculation since
        // only two nodes present 
        if(visited.size() ==1){
             for (int cityIdx = 0; cityIdx < unvisited.size(); ++cityIdx) {
                int city = unvisited[cityIdx];
                //can not form a cycle with 2 nodes so dont account for this here
                int64_t costBegin = problem.GetCostAndDistance(visited[0],city);
                if (costBegin < bestCost) {
                    bestCost = costBegin;
                    bestInsertPos = 1;
                    bestCityIdx = cityIdx;
                }

            }
        }
        else{
        // Try inserting each unvisited city at each possible position
            for (int cityIdx = 0; cityIdx < unvisited.size(); ++cityIdx) {
                int city = unvisited[cityIdx];
                // Try inserting at beginning (position 0)
                // no rerouting but change of the cycle closing with respect to teh previous iteration
                // state 
                int64_t costBegin = problem.GetCostAndDistance(visited[0],city)
                + problem.GetRawDistance(city,visited.back());
                int64_t costBeginDelta = costBegin -  problem.GetRawDistance(visited[0],visited.back());

                if (costBeginDelta < bestCost) {
                    bestCost = costBeginDelta;
                    bestInsertPos = 0;
                    bestCityIdx = cityIdx;
                }

                 //try inserting at the end, no need for rerouting but change of the cycle closing edge
                int64_t costEnd = problem.GetCostAndDistance(visited.back(), city)
                + problem.GetRawDistance(city,visited[0]);
                int64_t costEndDelta = costEnd - problem.GetRawDistance(visited[0],visited.back());
                if (costEndDelta < bestCost) {
                    bestCost = costEndDelta;
                    bestInsertPos = static_cast<int>(visited.size());
                    bestCityIdx = cityIdx;
                }

                // Try inserting in the middle of path
                for (int pos = 1; pos < visited.size(); ++pos) {
                    //needs to account for rerouting of the old path 
                    // no change of the clsoing of the cycle with respect to previous state

                    int64_t oldEdge = problem.GetCostAndDistance(visited[pos-1], visited[pos]);
                    int64_t newEdges = problem.GetCostAndDistance(visited[pos-1], city) + 
                                        problem.GetCostAndDistance(city, visited[pos]);
                    int64_t deltaCost = newEdges - oldEdge;
                    // no change for the cycle closing so do not, include it into the delta function
                    if (deltaCost < bestCost) {
                        bestCost = deltaCost;
                        bestInsertPos = pos;
                        bestCityIdx = cityIdx;
                    }
                }
               
            }
        }
    

        if (bestCityIdx == -1) break;
        // Insert the best city at the best position
        visited.insert(visited.begin() + bestInsertPos, unvisited[bestCityIdx]);
        unvisited.erase(unvisited.begin() + bestCityIdx);
    }
    AssertHamiltonian(visited,targetSize);
    return visited;
}
