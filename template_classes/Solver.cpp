#include <fstream>
#include <filesystem>
#include "Solver.h"
Solver::Solver(const ProblemInstance& prob)
    : problem(prob) {}

bool Solver::writePathCsv(const std::vector<int>& path,
                          const std::string& filename,
                          char sep,
                          bool append) {
    // Ensure target directory exists if filename contains a path
    try {
        auto p = std::filesystem::path(filename).parent_path();
        if (!p.empty()) std::filesystem::create_directories(p);
    } catch (...) {
        // ignore filesystem errors here; we'll still try to open the file
    }

    std::ofstream ofs;
    ofs.open(filename, append ? std::ios::app : std::ios::trunc);
    if (!ofs) return false;

    for (size_t i = 0; i < path.size(); ++i) {
        ofs << path[i];
        if (i + 1 < path.size()) ofs << sep;
    }
    ofs << '\n';
    return static_cast<bool>(ofs);
}

