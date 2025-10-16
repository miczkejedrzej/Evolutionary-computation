#include "ProblemInstance.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>


ProblemInstance::ProblemInstance(const std::string& file, int citiesInCycle, std::string name)
    : name(name),
citiesInCycle(citiesInCycle){
    setProblem(file);
}

void ProblemInstance::setProblem(const std::string& file) {
    std::ifstream f(file);
    if (!f) {
        std::cerr << "Error: cannot open file " << file << std::endl;
        return;
    }


    while (std::getline(f, line)) {
        std::stringstream ss(line);
        Row row;
        std::string value;

        std::getline(ss, value, ';'); row.x = std::stoll(value);
        std::getline(ss, value, ';'); row.y = std::stoll(value);
        std::getline(ss, value, ';'); row.cost = std::stoll(value);

        rows.push_back(row);
    }
    distances_no_cost.resize(rows.size(), std::vector<int64_t>(rows.size()));
    distances.resize(rows.size(), std::vector<int64_t>(rows.size()));
    for (size_t i = 0; i < rows.size(); ++i) {
        indices.push_back(i);
        for (size_t j = 0; j < rows.size(); ++j) {
            distances[i][j] = DistCalculation(rows[i], rows[j]);
            distances_no_cost[i][j] = RawDistCalculation(rows[i], rows[j]);
        }
    }
}

int64_t ProblemInstance::DistCalculation(const Row& row_start, const Row& row_end, DistanceType d) {
    int64_t dist = -1;
    
    switch (d)
    {
    case Euclidean:
        dist = int64_t(sqrt(pow(row_start.x - row_end.x, 2) + pow(row_start.y - row_end.y, 2)) + 0.5);
        break;
    
    default:
        break;
    }

    return dist + row_end.cost;
}

int64_t ProblemInstance::RawDistCalculation(const Row& row_start, const Row& row_end, DistanceType d) {
    int64_t dist = -1;
    
    switch (d)
    {
    case Euclidean:
        dist = int64_t(sqrt(pow(row_start.x - row_end.x, 2) + pow(row_start.y - row_end.y, 2)) + 0.5);
        break;
    
    default:
        break;
    }

    return dist ;
}




int64_t ProblemInstance::FullDistanceAndCost(std::vector<int> order) const {
    int64_t sum = 0;
    for (size_t i = 0; i < order.size() - 1; ++i)
        sum += distances[order[i]][order[i + 1]];

    sum += distances[order.back()][order.front()];
    return sum;
}

int64_t ProblemInstance::GetCostAndDistance(size_t from, size_t to) const {
    return distances[from][to];
}

int64_t ProblemInstance::GetRawDistance(size_t from, size_t to) const {
    return  distances_no_cost[from][to];
}
int64_t ProblemInstance::GetX(size_t index) const {
    return rows[index].x;
}
int64_t ProblemInstance::GetY(size_t index) const {
    return rows[index].y;
}
int64_t ProblemInstance::GetCost(size_t index) const {
    return rows[index].cost;
}

std::vector<int> ProblemInstance::GiveIndices() const{
    return indices;
}

int64_t ProblemInstance::GetNumberCitiesInCycle() const {
    return citiesInCycle;
}


