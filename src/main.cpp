// main.cpp
// Command-line front end: reads a stack-up CSV, runs all three methods,
// prints a comparison table, and can export a histogram or run a benchmark.
//
// Usage:
//   tolstack <file.csv> [--samples N] [--seed S] [--dist normal|uniform]
//                       [--hist out.csv] [--bins B] [--benchmark]

#include "stackup.hpp"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace tolstack;

namespace {

struct Options {
    std::string input;
    std::size_t samples = 1'000'000;
    std::uint64_t seed = 42;
    Distribution dist = Distribution::Normal;
    std::string hist_path;
    int bins = 40;
    bool benchmark = false;
};

void print_usage() {
    std::cout << "Usage: tolstack <file.csv> [--samples N] [--seed S] [--dist normal|uniform]\n"
                 "                [--hist out.csv] [--bins B] [--benchmark]\n";
}

Options parse_args(int argc, char** argv) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) throw std::runtime_error(arg + " needs a value");
            return argv[++i];
        };
        if (arg == "--samples") opt.samples = std::stoull(next());
        else if (arg == "--seed") opt.seed = std::stoull(next());
        else if (arg == "--dist") {
            const std::string d = next();
            if (d == "normal") opt.dist = Distribution::Normal;
            else if (d == "uniform") opt.dist = Distribution::Uniform;
            else throw std::runtime_error("--dist must be normal or uniform");
        }
        else if (arg == "--hist") opt.hist_path = next();
        else if (arg == "--bins") opt.bins = std::stoi(next());
        else if (arg == "--benchmark") opt.benchmark = true;
        else if (arg == "--help" || arg == "-h") { print_usage(); std::exit(0); }
        else if (!arg.empty() && arg[0] == '-') throw std::runtime_error("unknown option " + arg);
        else opt.input = arg;
    }
    if (opt.input.empty()) throw std::runtime_error("no input file given");
    if (opt.bins < 1) throw std::runtime_error("--bins must be at least 1");
    return opt;
}

std::string percent(double fraction) {
    if (fraction < 0.0) return "n/a";
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(4) << fraction * 100.0 << " %";
    return ss.str();
}

void print_table(const StackInput& in, const std::vector<Result>& results) {
    std::cout << "\nDimension chain (" << in.dims.size() << " parts):\n";
    for (const auto& d : in.dims) {
        std::cout << "  " << (d.direction > 0 ? "+ " : "- ") << std::left << std::setw(22) << d.name
                  << std::right << std::fixed << std::setprecision(3) << std::setw(9) << d.nominal
                  << " +/- " << d.tolerance << " mm\n";
    }
    std::cout << "Gap spec: " << in.spec.lower << " to " << in.spec.upper << " mm\n\n";

    std::cout << std::left << std::setw(24) << "Method" << std::right << std::setw(10) << "Min (mm)"
              << std::setw(10) << "Max (mm)" << std::setw(16) << "Out of spec" << std::setw(12)
              << "Time (ms)" << "\n";
    std::cout << std::string(72, '-') << "\n";
    for (const auto& r : results) {
        std::cout << std::left << std::setw(24) << r.method << std::right << std::fixed
                  << std::setprecision(3) << std::setw(10) << r.low << std::setw(10) << r.high
                  << std::setw(16) << percent(r.fraction_out) << std::setw(12) << std::setprecision(2)
                  << r.millis << "\n";
    }
    std::cout << "\nNote: for Worst case, 'Out of spec' is 100 % if any extreme combination fails\n"
                 "the spec, and 0 % if every possible assembly passes.\n";
}

void write_histogram(const std::string& path, const std::vector<double>& samples, int bins,
                     const Spec& spec) {
    const auto [min_it, max_it] = std::minmax_element(samples.begin(), samples.end());
    const double lo = *min_it, hi = *max_it;
    const double width = (hi > lo) ? (hi - lo) / bins : 1.0;

    std::vector<std::size_t> counts(bins, 0);
    for (double s : samples) {
        int b = static_cast<int>((s - lo) / width);
        counts[std::clamp(b, 0, bins - 1)]++;
    }

    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path);
    out << "bin_start_mm,bin_end_mm,count,in_spec\n";
    for (int b = 0; b < bins; ++b) {
        const double start = lo + b * width, end = start + width;
        const bool in_spec = start >= spec.lower && end <= spec.upper;
        out << std::fixed << std::setprecision(4) << start << "," << end << "," << counts[b] << ","
            << (in_spec ? "yes" : "no") << "\n";
    }
    std::cout << "Histogram written to " << path << " (" << bins << " bins)\n";
}

// Shows how run time and throughput scale with the number of simulated assemblies.
void run_benchmark(const StackInput& in, std::uint64_t seed, Distribution dist) {
    std::cout << "\nBenchmark: Monte Carlo run time vs. sample count\n";
    std::cout << std::right << std::setw(12) << "Samples" << std::setw(12) << "Time (ms)"
              << std::setw(22) << "Samples / second" << std::setw(16) << "Out of spec" << "\n";
    std::cout << std::string(62, '-') << "\n";
    for (std::size_t n : {10'000ULL, 100'000ULL, 1'000'000ULL, 10'000'000ULL}) {
        const auto stats = monte_carlo(in, n, seed, dist, false);
        const double per_sec = n / (stats.result.millis / 1000.0);
        std::cout << std::setw(12) << n << std::setw(12) << std::fixed << std::setprecision(2)
                  << stats.result.millis << std::setw(22) << std::setprecision(0) << per_sec
                  << std::setw(16) << percent(stats.result.fraction_out) << "\n";
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options opt = parse_args(argc, argv);
        const StackInput input = load_csv(opt.input);

        const bool keep = !opt.hist_path.empty();
        const auto mc = monte_carlo(input, opt.samples, opt.seed, opt.dist, keep);
        const std::vector<Result> results = {worst_case(input), rss(input), mc.result};

        print_table(input, results);
        std::cout << "Monte Carlo: " << opt.samples << " assemblies, seed " << opt.seed
                  << ", mean " << std::setprecision(4) << mc.mean << " mm, std dev " << mc.stddev
                  << " mm\n";

        if (keep) write_histogram(opt.hist_path, mc.samples, opt.bins, input.spec);
        if (opt.benchmark) run_benchmark(input, opt.seed, opt.dist);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        print_usage();
        return 1;
    }
}
