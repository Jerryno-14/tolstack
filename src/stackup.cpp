// stackup.cpp
// Implementation of the three stack-up methods and the CSV loader.

#include "stackup.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>

namespace tolstack {

namespace {

// Removes spaces and tabs from both ends of a string.
std::string trim(const std::string& s) {
    const auto first = s.find_first_not_of(" \t\r");
    if (first == std::string::npos) return "";
    const auto last = s.find_last_not_of(" \t\r");
    return s.substr(first, last - first + 1);
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string field;
    while (std::getline(ss, field, ',')) fields.push_back(trim(field));
    return fields;
}

double to_double(const std::string& text, int line_no) {
    try {
        std::size_t used = 0;
        const double value = std::stod(text, &used);
        if (used != text.size()) throw std::invalid_argument("extra characters");
        return value;
    } catch (const std::exception&) {
        throw std::runtime_error("line " + std::to_string(line_no) + ": '" + text +
                                 "' is not a number");
    }
}

// Share of a normal distribution N(mean, sigma) that falls outside [lower, upper].
double normal_fraction_outside(double mean, double sigma, double lower, double upper) {
    if (sigma <= 0.0) return (mean < lower || mean > upper) ? 1.0 : 0.0;
    const double below = 0.5 * std::erfc((mean - lower) / (sigma * std::sqrt(2.0)));
    const double above = 0.5 * std::erfc((upper - mean) / (sigma * std::sqrt(2.0)));
    return below + above;
}

double elapsed_ms(std::chrono::steady_clock::time_point start) {
    const auto end = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(end - start).count();
}

}  // namespace

StackInput load_csv(const std::string& path) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("cannot open file: " + path);

    StackInput input;
    bool has_spec = false;
    std::string line;
    int line_no = 0;

    while (std::getline(file, line)) {
        ++line_no;
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;  // blank line or comment

        const auto f = split_csv(line);

        if (f[0] == "spec") {  // spec,lower,upper
            if (f.size() != 3) throw std::runtime_error("line " + std::to_string(line_no) +
                                                        ": spec needs 'spec,lower,upper'");
            input.spec.lower = to_double(f[1], line_no);
            input.spec.upper = to_double(f[2], line_no);
            if (input.spec.lower >= input.spec.upper)
                throw std::runtime_error("line " + std::to_string(line_no) +
                                         ": spec lower must be below upper");
            has_spec = true;
            continue;
        }

        // name,nominal,tolerance,direction
        if (f.size() != 4) throw std::runtime_error("line " + std::to_string(line_no) +
                                                    ": expected 'name,nominal,tolerance,direction'");
        Dimension d;
        d.name = f[0];
        d.nominal = to_double(f[1], line_no);
        d.tolerance = to_double(f[2], line_no);
        const double dir = to_double(f[3], line_no);
        if (d.tolerance < 0.0)
            throw std::runtime_error("line " + std::to_string(line_no) + ": tolerance must be >= 0");
        if (dir != 1.0 && dir != -1.0)
            throw std::runtime_error("line " + std::to_string(line_no) + ": direction must be +1 or -1");
        d.direction = static_cast<int>(dir);
        input.dims.push_back(d);
    }

    if (!has_spec) throw std::runtime_error("missing 'spec,lower,upper' line");
    if (input.dims.empty()) throw std::runtime_error("no dimensions found");
    return input;
}

double nominal_gap(const std::vector<Dimension>& dims) {
    double gap = 0.0;
    for (const auto& d : dims) gap += d.direction * d.nominal;
    return gap;
}

Result worst_case(const StackInput& in) {
    const auto start = std::chrono::steady_clock::now();
    Result r;
    r.method = "Worst case";
    r.nominal = nominal_gap(in.dims);

    double total_tol = 0.0;
    for (const auto& d : in.dims) total_tol += d.tolerance;  // tolerances simply add up
    r.low = r.nominal - total_tol;
    r.high = r.nominal + total_tol;

    // Worst case gives a guarantee, not a probability: 0 if always in spec, else 1.
    r.fraction_out = (r.low < in.spec.lower || r.high > in.spec.upper) ? 1.0 : 0.0;
    r.millis = elapsed_ms(start);
    return r;
}

Result rss(const StackInput& in) {
    const auto start = std::chrono::steady_clock::now();
    Result r;
    r.method = "RSS (3-sigma)";
    r.nominal = nominal_gap(in.dims);

    // Each tolerance is treated as 3 standard deviations of a normal distribution.
    // Independent variances add, so the total sigma is the root of the summed squares.
    double sum_sq = 0.0;
    for (const auto& d : in.dims) sum_sq += d.tolerance * d.tolerance;
    const double total_tol = std::sqrt(sum_sq);
    r.low = r.nominal - total_tol;
    r.high = r.nominal + total_tol;
    r.fraction_out = normal_fraction_outside(r.nominal, total_tol / 3.0, in.spec.lower, in.spec.upper);
    r.millis = elapsed_ms(start);
    return r;
}

MonteCarloStats monte_carlo(const StackInput& in, std::size_t samples, std::uint64_t seed,
                            Distribution dist, bool keep_samples) {
    if (samples == 0) throw std::invalid_argument("samples must be > 0");
    const auto start = std::chrono::steady_clock::now();

    std::mt19937_64 rng(seed);  // fixed seed -> same results every run (reproducible)
    std::normal_distribution<double> unit_normal(0.0, 1.0);
    std::uniform_real_distribution<double> unit_uniform(-1.0, 1.0);

    MonteCarloStats stats;
    if (keep_samples) stats.samples.reserve(samples);

    double sum = 0.0, sum_sq = 0.0;
    double lo = 1e300, hi = -1e300;
    std::size_t out_of_spec = 0;

    for (std::size_t i = 0; i < samples; ++i) {
        // Build one random assembly: each part gets a random actual size.
        double gap = 0.0;
        for (const auto& d : in.dims) {
            const double deviation = (dist == Distribution::Normal)
                                         ? unit_normal(rng) * d.tolerance / 3.0
                                         : unit_uniform(rng) * d.tolerance;
            gap += d.direction * (d.nominal + deviation);
        }
        sum += gap;
        sum_sq += gap * gap;
        lo = std::min(lo, gap);
        hi = std::max(hi, gap);
        if (gap < in.spec.lower || gap > in.spec.upper) ++out_of_spec;
        if (keep_samples) stats.samples.push_back(gap);
    }

    const double n = static_cast<double>(samples);
    stats.mean = sum / n;
    stats.stddev = std::sqrt(std::max(0.0, sum_sq / n - stats.mean * stats.mean));

    Result& r = stats.result;
    r.method = (dist == Distribution::Normal) ? "Monte Carlo (normal)" : "Monte Carlo (uniform)";
    r.nominal = nominal_gap(in.dims);
    r.low = lo;
    r.high = hi;
    r.fraction_out = static_cast<double>(out_of_spec) / n;
    r.millis = elapsed_ms(start);
    return stats;
}

}  // namespace tolstack
