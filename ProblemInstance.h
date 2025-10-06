#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdint>

using namespace std;

vector<vector<int64_t>> distances;
// distances[nr from][nr to]

struct Row {
    int64_t x;
    int64_t y;
    int64_t cost;
};

enum DistanceType {
    Euclidean
};

int VectorSize (vector<Row> v) {
    return v.size() / sizeof(v[0]);
}

int64_t DistCalculation (Row row_start, Row row_end, DistanceType d = Euclidean) {
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

int64_t FullDistance (vector<int> order) {
    int64_t sum = 0;
    int length = order.size() / sizeof(order[0]);
    for (int i = 0; i < length - 1; i++) {
        sum += distances[order[i]][order[i + 1]];
    }

    return sum + distances[order[length - 1]][order[0]];
}

void ProblemInstance(string filepath) {
    ifstream file("data.csv");
    string line;
    vector<Row> rows;

    while (getline(file, line)) {
        stringstream ss(line);
        Row row;
        string value = "";
        
        getline(ss, value, ';');
        row.x = stoi(value);
        
        getline(ss, value, ';');
        row.y = stoi(value);
        
        getline(ss, value, ';');
        row.cost = stoi(value);
        
        rows.push_back(row);
    }

    for (int i = 0; i < VectorSize(rows); i++) {
        distances.push_back({});
        for (int j = 0; j < VectorSize(rows); j++) {
            distances[i].push_back(DistCalculation(rows[i], rows[j]));
        }
    }
}