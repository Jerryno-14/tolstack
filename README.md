# tolstack: Tolerance Stack-up Analyzer (C++17)

A command-line tool that predicts whether an assembly gap stays inside its
specification, given each part's nominal size and +/- tolerance. It compares
three methods on the same input:

| Method | Idea | Strength | Weakness |
|---|---|---|---|
| Worst case | Every part at its extreme at once | Guarantees 100 % of assemblies fit | Very pessimistic, forces tight (expensive) tolerances |
| RSS | Root-sum-square, tolerance = 3 sigma | Instant, realistic for well-controlled processes | Assumes normal, centred, independent parts |
| Monte Carlo | Simulates millions of random assemblies | Works for any distribution; gives an out-of-spec rate | Slower; result depends on sample count |

## Build

```bash
cmake -S . -B build
cmake --build build
cd build && ctest          # runs the unit tests
```

Requires a C++17 compiler (GCC, Clang, or MSVC) and CMake 3.16+.

## Run

```bash
./build/tolstack data/gear_shaft_axial.csv
./build/tolstack data/gear_shaft_axial_tight.csv --dist uniform
./build/tolstack data/gear_shaft_axial.csv --hist hist.csv --bins 40
./build/tolstack data/gear_shaft_axial.csv --benchmark
```

Options: `--samples N` (default 1,000,000), `--seed S` (default 42, fixed for
reproducible results), `--dist normal|uniform`, `--hist out.csv`, `--bins B`,
`--benchmark`.

## Input format

```
# comment lines start with #
spec,lower_mm,upper_mm
name,nominal_mm,tolerance_mm,direction
```

`direction` is `+1` if a larger part makes the gap larger and `-1` if it makes
the gap smaller. Bad input is reported with its line number.

## Example result

Axial end play of a gear shaft in a reducer housing, 6 parts (illustrative
values, not measured from a real design). With a tight 0.30 to 0.55 mm spec
(`data/gear_shaft_axial_tight.csv`, seed 42, 1,000,000 assemblies):

| Method | Min (mm) | Max (mm) | Out of spec |
|---|---|---|---|
| Worst case | 0.030 | 0.770 | fails |
| RSS (3-sigma) | 0.243 | 0.557 | 3.02 % |
| Monte Carlo, normal parts | 0.142 | 0.645 | 3.03 % |
| Monte Carlo, uniform parts | 0.064 | 0.729 | 19.1 % |

What this shows:
- Monte Carlo with normal parts agrees with RSS to within 0.01 percentage
  points, which validates both implementations against each other.
- If the parts are instead spread evenly across their tolerance band (a
  less controlled process), the real reject rate is about 6 times higher than
  RSS predicts. RSS is only safe when its assumptions hold.
- Worst case rejects the design outright, even though about 97 % of normally
  distributed assemblies would pass.

Monte Carlo throughput on one CPU core (6 parts per assembly, `--benchmark`,
Release build) was roughly 5 to 6 million assemblies per second in the test
environment. Run `--benchmark` on your own machine and quote your own number.

## Tests

`tests/test_stackup.cpp` checks the worst-case and RSS maths against hand
calculations, confirms Monte Carlo converges to the RSS result, confirms the
same seed gives identical samples, confirms uniform samples never leave the
worst-case band, and checks CSV error reporting.

## Project layout

```
src/stackup.hpp   data types and function declarations
src/stackup.cpp   CSV loader and the three methods
src/main.cpp      command-line interface, table, histogram, benchmark
tests/            unit tests
data/             example inputs
```
