#pragma once
#include "ProblemInstance.h"
#include <vector>
#include <random>
#include <iostream>

class Solver {
public:
    // Constructor: takes a reference to the problem instance
    explicit Solver(const ProblemInstance& prob);
    // Virtual destructor 
    virtual ~Solver() = default;
    // function to overwrite 
    virtual std::vector<int> solve() = 0;
protected:
    // for the Solvers to use 
    const ProblemInstance& problem;
};
