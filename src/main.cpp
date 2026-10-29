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

using namespace std;
using namespace tolstack;

namespace {

struct Options {
    string input;
    size_t samples = 1'000'000;
    uint64_t seed = 42;
    Distribution dist = Distribution::Normal;
    string hist_path;
    int bins = 40;
    bool benchmark = false;
};

void print_usage() {
    cout << "Usage: tolstack <file.csv> [--samples N] [--seed S] [--dist normal|uniform]\n"
                 "                [--hist out.csv] [--bins B] [--benchmark]\n";
}

Options parse_args(int argc, char** argv) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const string arg = argv[i];
        auto next = [&]() -> string {
            if (i + 1 >= argc) throw runtime_error(arg + " needs a value");
            return argv[++i];
        };
        if (arg == "--samples") opt.samples = stoull(next());
        else if (arg == "--seed") opt.seed = stoull(next());
        else if (arg == "--dist") {
            const string d = next();
            if (d == "normal") opt.dist = Distribution::Normal;
            else if (d == "uniform") opt.dist = Distribution::Uniform;
            else throw runtime_error("--dist must be normal or uniform");
        }
        else if (arg == "--hist") opt.hist_path = next();
        else if (arg == "--bins") opt.bins = stoi(next());
        else if (arg == "--benchmark") opt.benchmark = true;
        else if (arg == "--help" || arg == "-h") { print_usage(); exit(0); }
        else if (!arg.empty() && arg[0] == '-') throw runtime_error("unknown option " + arg);
        else opt.input = arg;
    }
    if (opt.input.empty()) throw runtime_error("no input file given");
    if (opt.bins < 1) throw runtime_error("--bins must be at least 1");
    return opt;
}

string percent(double fraction) {
    if (fraction < 0.0) return "n/a";
    ostringstream ss;
    ss << fixed << setprecision(4) << fraction * 100.0 << " %";
    return ss.str();
}

void print_table(const StackInput& in, const vector<Result>& results) {
    cout << "\nDimension chain (" << in.dims.size() << " parts):\n";
    for (const auto& d : in.dims) {
        cout << "  " << (d.direction > 0 ? "+ " : "- ") << left << setw(22) << d.name
                  << right << fixed << setprecision(3) << setw(9) << d.nominal
                  << " +/- " << d.tolerance << " mm\n";
    }
    cout << "Gap spec: " << in.spec.lower << " to " << in.spec.upper << " mm\n\n";

    cout << left << setw(24) << "Method" << right << setw(10) << "Min (mm)"
              << setw(10) << "Max (mm)" << setw(16) << "Out of spec" << setw(12)
              << "Time (ms)" << "\n";
    cout << string(72, '-') << "\n";
    for (const auto& r : results) {
        cout << left << setw(24) << r.method << right << fixed
                  << setprecision(3) << setw(10) << r.low << setw(10) << r.high
                  << setw(16) << percent(r.fraction_out) << setw(12) << setprecision(2)
                  << r.millis << "\n";
    }
    cout << "\nNote: for Worst case, 'Out of spec' is 100 % if any extreme combination fails\n"
                 "the spec, and 0 % if every possible assembly passes.\n";
}

void write_histogram(const string& path, const vector<double>& samples, int bins,
                     const Spec& spec) {
    const auto [min_it, max_it] = minmax_element(samples.begin(), samples.end());
    const double lo = *min_it, hi = *max_it;
    const double width = (hi > lo) ? (hi - lo) / bins : 1.0;

    vector<size_t> counts(bins, 0);
    for (double s : samples) {
        int b = static_cast<int>((s - lo) / width);
        counts[clamp(b, 0, bins - 1)]++;
    }

    ofstream out(path);
    if (!out) throw runtime_error("cannot write " + path);
    out << "bin_start_mm,bin_end_mm,count,in_spec\n";
    for (int b = 0; b < bins; ++b) {
        const double start = lo + b * width, end = start + width;
        const bool in_spec = start >= spec.lower && end <= spec.upper;
        out << fixed << setprecision(4) << start << "," << end << "," << counts[b] << ","
            << (in_spec ? "yes" : "no") << "\n";
    }
    cout << "Histogram written to " << path << " (" << bins << " bins)\n";
}

// Shows how run time and throughput scale with the number of simulated assemblies.
void run_benchmark(const StackInput& in, uint64_t seed, Distribution dist) {
    cout << "\nBenchmark: Monte Carlo run time vs. sample count\n";
    cout << right << setw(12) << "Samples" << setw(12) << "Time (ms)"
              << setw(22) << "Samples / second" << setw(16) << "Out of spec" << "\n";
    cout << string(62, '-') << "\n";
    for (size_t n : {10'000ULL, 100'000ULL, 1'000'000ULL, 10'000'000ULL}) {
        const auto stats = monte_carlo(in, n, seed, dist, false);
        const double per_sec = n / (stats.result.millis / 1000.0);
        cout << setw(12) << n << setw(12) << fixed << setprecision(2)
                  << stats.result.millis << setw(22) << setprecision(0) << per_sec
                  << setw(16) << percent(stats.result.fraction_out) << "\n";
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options opt = parse_args(argc, argv);
        const StackInput input = load_csv(opt.input);

        const bool keep = !opt.hist_path.empty();
        const auto mc = monte_carlo(input, opt.samples, opt.seed, opt.dist, keep);
        const vector<Result> results = {worst_case(input), rss(input), mc.result};

        print_table(input, results);
        cout << "Monte Carlo: " << opt.samples << " assemblies, seed " << opt.seed
                  << ", mean " << setprecision(4) << mc.mean << " mm, std dev " << mc.stddev
                  << " mm\n";

        if (keep) write_histogram(opt.hist_path, mc.samples, opt.bins, input.spec);
        if (opt.benchmark) run_benchmark(input, opt.seed, opt.dist);
        return 0;
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << "\n";
        print_usage();
        return 1;
    }
}
