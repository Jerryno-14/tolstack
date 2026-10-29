// stackup.hpp
// Core data types and analysis functions for a 1D tolerance stack-up.
//
// A "stack-up" is a chain of part dimensions that together decide one
// critical result, usually a gap or clearance in an assembly. Each part has a
// nominal size and a +/- tolerance, so the real gap varies from assembly to
// assembly. This module predicts that variation three ways:
//   1. Worst case  - every part at its extreme at the same time (very safe, pessimistic)
//   2. RSS         - root-sum-square, a statistical estimate (assumes tolerance = 3 sigma)
//   3. Monte Carlo - simulate many random assemblies and count the bad ones

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace tolstack {

// One dimension in the chain.
struct Dimension {
    std::string name;
    double nominal = 0.0;    // mm
    double tolerance = 0.0;  // symmetric +/- tolerance, mm
    int direction = +1;      // +1 if a bigger part makes the gap bigger, -1 if smaller
};

// Allowed range for the gap.
struct Spec {
    double lower = 0.0;  // mm
    double upper = 0.0;  // mm
};

struct StackInput {
    std::vector<Dimension> dims;
    Spec spec;
};

enum class Distribution { Normal, Uniform };

// Summary returned by every method so they can be printed side by side.
struct Result {
    std::string method;
    double nominal = 0.0;        // nominal gap, mm
    double low = 0.0;            // predicted minimum gap, mm
    double high = 0.0;           // predicted maximum gap, mm
    double fraction_out = -1.0;  // estimated share of assemblies outside spec (-1 = not estimated)
    double millis = 0.0;         // run time, ms
};

struct MonteCarloStats {
    Result result;
    double mean = 0.0;
    double stddev = 0.0;
    std::vector<double> samples;  // only filled when keep_samples is true
};

// Reads a CSV file. Throws std::runtime_error with the line number on bad input.
StackInput load_csv(const std::string& path);

// Sum of direction * nominal.
double nominal_gap(const std::vector<Dimension>& dims);

Result worst_case(const StackInput& in);
Result rss(const StackInput& in);
MonteCarloStats monte_carlo(const StackInput& in, std::size_t samples, std::uint64_t seed,
                            Distribution dist, bool keep_samples);

}  // namespace tolstack
