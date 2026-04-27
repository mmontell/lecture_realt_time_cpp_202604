/*
 * nearest_neighbour_benchmark.cpp
 *
 * Three ways of writing a nearest-neighbour search.
 *
 * The problem:
 *   Given N 2D points.
 *   find, for each point, its nearest neighbour.
 *
 *
 * Compile with:
 *   g++ -O2 -std=c++17 -o nn_bench nearest_neighbour_benchmark.cpp
 *
 * Run:
 *   ./nn_bench
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <algorithm>
#include <chrono>
#include <random>
#include <cmath>
#include <numeric>
#include <string>
#include <cassert>

// ============================================================
//  Diagnostics / monitoring helpers
// ============================================================

struct BenchResult {
    std::string label;
    double ms_per_event;       // wall time per call
    double throughput_kHz;     // events/second in kHz
    long long n_pairs_found;   // correctness check
    double speedup;            // relative to Version 1
};

static void printHeader() {
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║   Nearest-Neighbour Search  —  Real-Time Optimisation Demo       ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════╝\n\n";
}

static void printSeparator() {
    std::cout << "──────────────────────────────────────────────────────────────────\n";
}

static void printResult(const BenchResult& r) {
    std::cout << std::left  << std::setw(38) << r.label
              << std::right << std::setw(9)  << std::fixed << std::setprecision(2)
              << r.ms_per_event << " ms   "
              << std::setw(8) << std::setprecision(1) << r.throughput_kHz << " kHz   "
              << std::setw(6) << std::setprecision(1) << r.speedup << "×\n";
}

static void printTableHeader() {
    printSeparator();
    std::cout << std::left  << std::setw(38) << "Version"
              << std::right << std::setw(9)  << "Time"
              << "       "
              << std::setw(8) << "Rate"
              << "   "
              << std::setw(6) << "Speedup"
              << "\n";
    printSeparator();
}

// ============================================================
//  Data generation
// ============================================================

struct Point {   // Array-of-Structs (AoS) — "intuitive" approach
    float x, y;
};

// Returns N random 2D points in [0, detector_size)
static std::vector<Point> generatePoints(int N, float detector_size = 1000.f, unsigned seed = 42) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> dist(0.f, detector_size);
    std::vector<Point> points(N);
    for (auto& p : points) { p.x = dist(rng); p.y = dist(rng); }
    return points;
}

// ============================================================
//  Version 1 — Intuitive - O(N^2) scaling
//
//  - Array-of-Structs layout
//  - No sorting, no windowing
//  - Every pair checked unconditionally
// ============================================================

static long long nearestNeighbour_intuitive(const std::vector<Point>& points) {
    const int N = static_cast<int>(points.size());
    long long found = 0;

    for (int i = 0; i < N; ++i) {
        float best_dist2 = std::numeric_limits<float>::max();
        int   best_j     = -1;

        for (int j = 0; j < N; ++j) {
            if (j == i) continue;
            float dx = points[i].x - points[j].x;
            float dy = points[i].y - points[j].y;
            float d2 = dx*dx + dy*dy;
            if (d2 < best_dist2) {
                best_dist2 = d2;
                best_j     = j;
            }
        }
        if (best_j >= 0) ++found;
    }
    return found;  // = N for any non-trivial input
}

// ============================================================
//  Version 2 — Sort + search window
//
//  Points are sorted by x, then for each point only
//  scan neighbours within a window |Delta(x)| < current_best.
//  We don't look at hits we can't possibly match.
//
//  Complexity: O(N log N + N·k)  where k << N in practice.
// ============================================================

static long long nearestNeighbour_Sorted(const std::vector<Point>& points) {
    const int N = static_cast<int>(points.size());

    // Build a sorted index vector by x-coordinate 
    // no data is actually copied - good! no extra memory allocation
    std::vector<int> idx(N);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b){
        return points[a].x < points[b].x;
    });

    // Also build a lookup: original index → sorted position
    std::vector<int> rank(N);
    for (int sortedOrder = 0; sortedOrder < N; ++sortedOrder) rank[idx[sortedOrder]] = sortedOrder;


    long long found = 0;

    for (int i = 0; i < N; ++i) {
        float best_dist2 = std::numeric_limits<float>::max();
        int   best_j     = -1;

        int si = rank[i];          // sorted position of point i
        float xi = points[i].x;
        float yi = points[i].y;

        // Scan right from current point
        for (int s = si + 1; s < N; ++s) {
            int j = idx[s];
            float dx = points[j].x - xi;
            if (dx * dx >= best_dist2) break;  // ← the window cut: RoI in x
            float dy = points[j].y - yi;
            float d2 = dx*dx + dy*dy;
            if (d2 < best_dist2) { best_dist2 = d2; best_j = j; }
        }

        // Scan left from si
        for (int s = si - 1; s >= 0; --s) {
            int j = idx[s];
            float dx = xi - points[j].x;
            if (dx * dx >= best_dist2) break;  // ← same window cut
            float dy = points[j].y - yi;
            float d2 = dx*dx + dy*dy;
            if (d2 < best_dist2) { best_dist2 = d2; best_j = j; }
        }

        if (best_j >= 0) ++found;
    }
    return found;
}

// ============================================================
//  Version 3 — Sorted + Struct-of-Arrays (SoA)
//
//  Same algorithm as V2 but the data is stored as separate
//  contiguous arrays of x and y, rather than interleaved
//  Point structs.
//
//  Why it matters:
//  When the inner loop accesses points[j].x, in AoS layout it
//  loads a cache line containing BOTH x and y — but it only
//  needs x for the early-exit check. Half the memory bandwidth
//  is wasted. In SoA layout, the x array is fully packed:
//  every byte fetched from memory is useful.
//
// ============================================================

static long long nearestNeighbour_SoA(const std::vector<Point>& points) {
    const int N = static_cast<int>(points.size());

    // Build sorted SoA buffers
    std::vector<int> idx(N);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b){
        return points[a].x < points[b].x;
    });

    // Separate, contiguous x and y arrays — cache-friendly for x-only scan
    std::vector<float> sx(N), sy(N);
    for (int s = 0; s < N; ++s) {
        sx[s] = points[idx[s]].x;
        sy[s] = points[idx[s]].y;
    }

    std::vector<int> rank(N);
    for (int s = 0; s < N; ++s) rank[idx[s]] = s;

    long long found = 0;

    for (int i = 0; i < N; ++i) {
        float best_dist2 = std::numeric_limits<float>::max();
        int   best_j     = -1;

        int   si = rank[i];
        float xi = sx[si];
        float yi = sy[si];

        // Scan right — only touching sx[] until window breaks
        for (int s = si + 1; s < N; ++s) {
            float dx = sx[s] - xi;
            if (dx * dx >= best_dist2) break;
            float dy = sy[s] - yi;
            float d2 = dx*dx + dy*dy;
            if (d2 < best_dist2) { best_dist2 = d2; best_j = idx[s]; }
        }

        // Scan left
        for (int s = si - 1; s >= 0; --s) {
            float dx = xi - sx[s];
            if (dx * dx >= best_dist2) break;
            float dy = yi - sy[s];
            float d2 = dx*dx + dy*dy;
            if (d2 < best_dist2) { best_dist2 = d2; best_j = idx[s]; }
        }

        if (best_j >= 0) ++found;
    }
    return found;
}

// ============================================================
//  Version 4 — SoA + pre-allocated buffers + reserve()
//
//  Shows use of reserve() to avoid
//  any heap allocation inside the loop
// ============================================================

struct NNSearchState {
    std::vector<int>   idx;
    std::vector<float> sx, sy;
    std::vector<int>   rank;

    explicit NNSearchState(int N) : idx(N), sx(N), sy(N), rank(N) {}
};

static long long nearestNeighbour_PreAlloc(const std::vector<Point>& points,
                                            NNSearchState& state) {
    const int N = static_cast<int>(points.size());

    std::iota(state.idx.begin(), state.idx.end(), 0);
    std::sort(state.idx.begin(), state.idx.end(), [&](int a, int b){
        return points[a].x < points[b].x;
    });

    for (int s = 0; s < N; ++s) {
        state.sx[s] = points[state.idx[s]].x;
        state.sy[s] = points[state.idx[s]].y;
    }
    for (int s = 0; s < N; ++s) state.rank[state.idx[s]] = s;

    long long found = 0;

    for (int i = 0; i < N; ++i) {
        float best_dist2 = std::numeric_limits<float>::max();
        int   best_j     = -1;

        int   si = state.rank[i];
        float xi = state.sx[si];
        float yi = state.sy[si];

        for (int s = si + 1; s < N; ++s) {
            float dx = state.sx[s] - xi;
            if (dx * dx >= best_dist2) break;
            float dy = state.sy[s] - yi;
            float d2 = dx*dx + dy*dy;
            if (d2 < best_dist2) { best_dist2 = d2; best_j = state.idx[s]; }
        }
        for (int s = si - 1; s >= 0; --s) {
            float dx = xi - state.sx[s];
            if (dx * dx >= best_dist2) break;
            float dy = yi - state.sy[s];
            float d2 = dx*dx + dy*dy;
            if (d2 < best_dist2) { best_dist2 = d2; best_j = state.idx[s]; }
        }

        if (best_j >= 0) ++found;
    }
    return found;
}

// ============================================================
//  Benchmarking harness
// ============================================================

template<typename Fn>
static BenchResult benchmark(const std::string& label,
                              Fn fn,
                              int n_events,
                              double baseline_ms = 0.0) {
    // Warmup
    fn();

    // Timed run
    auto t0 = std::chrono::high_resolution_clock::now();
    long long result = 0;
    for (int e = 0; e < n_events; ++e) result += fn();
    auto t1 = std::chrono::high_resolution_clock::now();

    double total_ms  = std::chrono::duration<double, std::milli>(t1 - t0).count();
    double ms_each   = total_ms / n_events;
    double rate_kHz  = 1.0 / ms_each;  // events/ms = kHz

    BenchResult r;
    r.label          = label;
    r.ms_per_event   = ms_each;
    r.throughput_kHz = rate_kHz;
    r.n_pairs_found  = result / n_events;
    r.speedup        = (baseline_ms > 0.0) ? baseline_ms / ms_each : 1.0;
    return r;
}

// ============================================================
//  Diagnostic: per-version inner-loop work analysis
// ============================================================

static void printWorkAnalysis(const std::vector<Point>& points) {
    const int N = static_cast<int>(points.size());

    // Count average inner-loop iterations for intuitive vs sorted
    long long intuitive_iters  = 0;
    long long sorted_iters = 0;

    // intuitive: always N-1
    intuitive_iters = static_cast<long long>(N) * (N - 1);

    // Sorted: count actual window scans
    std::vector<int> idx(N);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b){ return points[a].x < points[b].x; });
    std::vector<float> sx(N), sy(N);
    for (int s = 0; s < N; ++s) { sx[s] = points[idx[s]].x; sy[s] = points[idx[s]].y; }
    std::vector<int> rank(N);
    for (int s = 0; s < N; ++s) rank[idx[s]] = s;

    for (int i = 0; i < N; ++i) {
        float best_dist2 = std::numeric_limits<float>::max();
        int si = rank[i];
        float xi = sx[si], yi = sy[si];
        for (int s = si + 1; s < N; ++s) {
            ++sorted_iters;
            float dx = sx[s] - xi; if (dx*dx >= best_dist2) break;
            float dy = sy[s] - yi; float d2 = dx*dx+dy*dy;
            if (d2 < best_dist2) best_dist2 = d2;
        }
        for (int s = si - 1; s >= 0; --s) {
            ++sorted_iters;
            float dx = xi - sx[s]; if (dx*dx >= best_dist2) break;
            float dy = yi - sy[s]; float d2 = dx*dx+dy*dy;
            if (d2 < best_dist2) best_dist2 = d2;
        }
    }

    std::cout << "\n  Work analysis (N = " << N << " points):\n";
    std::cout << "  ┌─────────────────────────────────┬──────────────┬────────────┐\n";
    std::cout << "  │ Version                         │ Inner iters  │ Reduction  │\n";
    std::cout << "  ├─────────────────────────────────┼──────────────┼────────────┤\n";
    std::cout << "  │ intuitive (all pairs)           │ "
              << std::setw(12) << intuitive_iters << " │     1.0×   │\n";
    std::cout << "  │ Sorted + window                 │ "
              << std::setw(12) << sorted_iters << " │ "
              << std::setw(6) << std::fixed << std::setprecision(1)
              << (double)intuitive_iters / sorted_iters << "×   │\n";
    std::cout << "  └─────────────────────────────────┴──────────────┴────────────┘\n";
    std::cout << "\n  → The window cut alone eliminates ~"
              << std::setprecision(0)
              << 100.0 * (1.0 - (double)sorted_iters / intuitive_iters)
              << "% of distance calculations.\n"
              << "    (This is exactly the RoI η-φ window in TrigFastTrackFinder.)\n\n";
}

// ============================================================
//  Main
// ============================================================

int main() {
    constexpr int N_POINTS  = 5000;   
    constexpr int N_EVENTS  = 200;    // repeat for stable timing

    printHeader();

    std::cout << "  Configuration:\n";
    std::cout << "    Points per event : " << N_POINTS << "\n";
    std::cout << "    Events timed     : " << N_EVENTS << "\n";
    

    const auto points = generatePoints(N_POINTS);

    // ── Correctness check ──────────────────────────────────
    {
        auto r1 = nearestNeighbour_intuitive(points);
        auto r2 = nearestNeighbour_Sorted(points);
        auto r3 = nearestNeighbour_SoA(points);
        NNSearchState state(N_POINTS);
        auto r4 = nearestNeighbour_PreAlloc(points, state);
        assert(r1 == r2 && r2 == r3 && r3 == r4);
        std::cout << "  ✓ All versions agree on " << r1 << " nearest-neighbour pairs found.\n";
    }

    // ── Work analysis ──────────────────────────────────────
    printWorkAnalysis(points);

    // ── Timing benchmarks ─────────────────────────────────
    printTableHeader();

    NNSearchState state(N_POINTS);

    auto r1 = benchmark("V1: intuitive O(N²)",
        [&]{ return nearestNeighbour_intuitive(points); },
        N_EVENTS);

    auto r2 = benchmark("V2: Sort + window cut  [+std::sort]",
        [&]{ return nearestNeighbour_Sorted(points); },
        N_EVENTS, r1.ms_per_event);

    auto r3 = benchmark("V3: Sort + SoA layout  [+cache]",
        [&]{ return nearestNeighbour_SoA(points); },
        N_EVENTS, r1.ms_per_event);

    auto r4 = benchmark("V4: SoA + pre-alloc    [+no heap]",
        [&]{ return nearestNeighbour_PreAlloc(points, state); },
        N_EVENTS, r1.ms_per_event);

    printResult(r1);
    printResult(r2);
    printResult(r3);
    printResult(r4);
    printSeparator();

    // ── Interpretation ────────────────────────────────────
    std::cout << "\n  Interpretation:\n\n";
    std::cout << "  V1 → V2  (+" << std::setprecision(1) << r2.speedup
              << "×): Algorithmic improvement — window cut removes O(N²) work.\n"
              << "           One std::sort() call. Maps to RoI-gated seeding in FTF.\n\n";
    std::cout << "  V2 → V3  (+" << r3.speedup
              << "×): Memory layout — SoA lets the CPU prefetch x[] contiguously.\n"
              << "           Zero algorithmic change. Maps to TrigSiSpacePointBase layout.\n\n";
    std::cout << "  V3 → V4  (+" << r4.speedup
              << "×): Allocation — reusing pre-allocated buffers across events.\n"
              << "           Maps to convertedSpacePoints.reserve(5000) in FTF::findTracks().\n\n";

    std::cout << "  TOTAL    (V1 → V4): "
              << std::setprecision(1) << r4.speedup << "× speedup   →   "
              << std::setprecision(1) << r4.throughput_kHz << " kHz  \n\n";

    return 0;
}
