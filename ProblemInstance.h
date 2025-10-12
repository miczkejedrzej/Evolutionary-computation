


#pragma once
#include <vector>
#include <string>
#include <cstdint>

struct Row {
    int64_t x;
    int64_t y;
    int64_t cost;
};

enum DistanceType {
    Euclidean
};

class ProblemInstance {
public:

    std::string name;

    //ProblemInstance() = default;                     
    explicit ProblemInstance(const std::string& file, int citiesInCycle, std::string name); 

    void setProblem(const std::string& file);        

 
    int64_t FullDistanceAndCost(std::vector<int> order) const;

     size_t getNumCities() const { return distances.size(); }

    int64_t GetCostAndDistance(size_t from, size_t to) const;
    
    int64_t GetRawDistance(size_t from, size_t to) const;

    int64_t GetX(size_t index) const;

    int64_t GetY(size_t index) const;

    int64_t GetCost(size_t index) const;

    int64_t GetNumberCitiesInCycle() const;

    //Row yieldRow(size_t index) const;
    //std::vector<Row> yieldRows(size_t index) const;

    std::vector<int> GiveIndices() const;

    static int64_t DistCalculation(const Row& a, const Row& b, DistanceType d = Euclidean);

    static int64_t RawDistCalculation(const Row& row_start, const Row& row_end, DistanceType d = Euclidean);

private:
    std::vector<std::vector<int64_t>> distances;
    std::vector<std::vector<int64_t>> distances_no_cost;
    int64_t citiesInCycle;
    std::vector<int> indices;
    std::string line;
    std::vector<Row> rows;
    int64_t HamiltonianNumber;
};
