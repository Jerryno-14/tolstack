// test_stackup.cpp
// Small self-contained tests (no external framework needed).
// Build and run with: ctest  (or ./test_stackup)

#include "../src/stackup.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace tolstack;

static int failures = 0;

#define CHECK(cond)                                                              \
    do {                                                                         \
        if (!(cond)) {                                                           \
            std::cerr << "FAIL line " << __LINE__ << ": " #cond "\n";            \
            ++failures;                                                          \
        }                                                                        \
    } while (0)

static bool near(double a, double b, double eps = 1e-9) { return std::fabs(a - b) <= eps; }

// Simple chain: 10 +/- 0.3 minus 9 +/- 0.4, spec 0.5 to 1.5.
static StackInput simple_input() {
    StackInput in;
    in.dims = {{"outer", 10.0, 0.3, +1}, {"inner", 9.0, 0.4, -1}};
    in.spec = {0.5, 1.5};
    return in;
}

static void test_nominal_and_worst_case() {
    const auto in = simple_input();
    CHECK(near(nominal_gap(in.dims), 1.0));
    const auto r = worst_case(in);
    CHECK(near(r.low, 0.3));   // 1.0 - (0.3 + 0.4)
    CHECK(near(r.high, 1.7));
    CHECK(near(r.fraction_out, 1.0));  // extremes break the spec
}

static void test_rss() {
    const auto r = rss(simple_input());
    CHECK(near(r.low, 0.5));   // sqrt(0.3^2 + 0.4^2) = 0.5
    CHECK(near(r.high, 1.5));
    // Spec edges sit exactly at +/- 3 sigma -> about 0.27 % outside.
    CHECK(std::fabs(r.fraction_out - 0.0027) < 0.0001);
}

static void test_monte_carlo_matches_rss() {
    const auto in = simple_input();
    const auto mc = monte_carlo(in, 2'000'000, 7, Distribution::Normal, false);
    CHECK(std::fabs(mc.mean - 1.0) < 0.001);
    CHECK(std::fabs(mc.stddev - 0.5 / 3.0) < 0.001);
    CHECK(std::fabs(mc.result.fraction_out - 0.0027) < 0.0003);
}

static void test_monte_carlo_is_reproducible() {
    const auto in = simple_input();
    const auto a = monte_carlo(in, 10'000, 123, Distribution::Normal, true);
    const auto b = monte_carlo(in, 10'000, 123, Distribution::Normal, true);
    CHECK(a.samples == b.samples);  // same seed -> identical samples
}

static void test_uniform_stays_inside_worst_case() {
    const auto in = simple_input();
    const auto mc = monte_carlo(in, 100'000, 1, Distribution::Uniform, false);
    CHECK(mc.result.low >= 0.3 - 1e-12);
    CHECK(mc.result.high <= 1.7 + 1e-12);
}

static void test_csv_loading_and_errors() {
    const std::string good = "test_good.csv";
    {
        std::ofstream f(good);
        f << "# comment\nspec,0.5,1.5\nouter,10,0.3,+1\ninner, 9 ,0.4,-1\n";
    }
    const auto in = load_csv(good);
    CHECK(in.dims.size() == 2);
    CHECK(in.dims[1].direction == -1);
    CHECK(near(in.spec.upper, 1.5));
    std::remove(good.c_str());

    const std::string bad = "test_bad.csv";
    {
        std::ofstream f(bad);
        f << "spec,0.5,1.5\nouter,abc,0.3,+1\n";
    }
    bool threw = false;
    try { load_csv(bad); } catch (const std::runtime_error& e) {
        threw = std::string(e.what()).find("line 2") != std::string::npos;
    }
    CHECK(threw);
    std::remove(bad.c_str());
}

int main() {
    test_nominal_and_worst_case();
    test_rss();
    test_monte_carlo_matches_rss();
    test_monte_carlo_is_reproducible();
    test_uniform_stays_inside_worst_case();
    test_csv_loading_and_errors();

    if (failures == 0) std::cout << "All tests passed.\n";
    else std::cout << failures << " check(s) failed.\n";
    return failures == 0 ? 0 : 1;
}
