/*
 * nearest_neighbour_benchmark.cpp
 *
 * A pedagogical example for: "C++ in Real-Time Analysis"
 * Demonstrates the interplay between memory optimisations and
 * algorithmic improvement in a nearest-neighbour search.
 *
 * The problem:
 *   Given N 2D points (e.g. detector hits in the transverse plane),
 *   find, for each point, its nearest neighbour.
 *
 * This is structurally identical to doublet-seed finding in a tracker:
 *   - the RoI window  → the search window cutoff
 *   - hit collections → space point collections
 *   - AoS vs SoA     → layout of TrigSiSpacePointBase collections
 *
 * Narrative arc (the pedagogical point):
 *   V1 → V2  Memory layout:  AoS  →  SoA               (modest gain)
 *   V2 → V3  Allocation:     heap per event → pre-alloc (marginal gain)
 *   V3 → V4  Algorithm:      O(N²) → sort + window cut  (dominant gain)
 *
 *   Good memory hygiene is worth doing — but it cannot rescue a bad algorithm.
 *   The sort turns an O(N²) problem into O(N log N + N·k), k ≪ N.
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
#include <algorithm>
#include <chrono>
#include <random>
#include <numeric>
#include <string>
#include <cassert>
#include <limits>

// ============================================================
//  Diagnostics / monitoring helpers
// ============================================================

struct BenchResult {
    std::string label;
    double ms_per_event;
    double throughput_kHz;
    long long n_pairs_found;
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
    std::cout << std::left  << std::setw(44) << r.label
              << std::right << std::setw(9)  << std::fixed << std::setprecision(2)
              << r.ms_per_event << " ms   "
              << std::setw(8) << std::setprecision(1) << r.throughput_kHz << " kHz   "
              << std::setw(6) << std::setprecision(1) << r.speedup << "×\n";
}

static void printTableHeader() {
    printSeparator();
    std::cout << std::left  << std::setw(44) << "Version"
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

struct Point {   // Array-of-Structs (AoS) — intuitive layout
    float x, y;
};

// Generate N random points; seed varies per event to prevent
// the compiler from memoising results across repeated calls.
static std::vector<Point> generatePoints(int N, float detector_size = 1000.f,
                                          unsigned seed = 42) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> dist(0.f, detector_size);
    std::vector<Point> pts(N);
    for (auto& p : pts) { p.x = dist(rng); p.y = dist(rng); }
    return pts;
}

// ============================================================
//  Version 1 — Intuitive O(N²), Array-of-Structs
//
//  What a student just out of their C++ course would write.
//  Correct. Clear. Slow.
//  - AoS layout: Point{x, y} packed together in memory.
//    Loading pts[j].x fetches a cache line that also contains
//    pts[j].y — only half the data we paid for is used right now.
//  - No sorting, no windowing: every pair checked unconditionally.
// ============================================================

__attribute__((noinline))
static long long nearestNeighbour_Intuitive(const std::vector<Point>& pts) {
    const int N = static_cast<int>(pts.size());
    long long found = 0;

    for (int i = 0; i < N; ++i) {
        float best_dist2 = std::numeric_limits<float>::max();
        int   best_j     = -1;

        for (int j = 0; j < N; ++j) {
            if (j == i) continue;
            float dx = pts[i].x - pts[j].x;
            float dy = pts[i].y - pts[j].y;
            float d2 = dx*dx + dy*dy;
            if (d2 < best_dist2) {
                best_dist2 = d2;
                best_j     = j;
            }
        }
        if (best_j >= 0) ++found;
    }
    return found;
}

// ============================================================
//  Version 2 — Intuitive O(N²), Struct-of-Arrays
//
//  Same brute-force algorithm as V1. One change only: x and y
//  are stored in separate contiguous arrays (SoA).
//
//  Why it helps:
//  The inner loop checks dx² against best_dist² before it needs
//  dy at all. In AoS, loading pts[j].x drags pts[j].y along for
//  the ride — we pay for a full cache line even when the x check
//  immediately rejects the candidate. In SoA, sx[] and sy[] are
//  separate arrays: when dx² ≥ best_dist², we skip the sy[] load
//  entirely. Every cache line fetched from sx[] contains only
//  useful x values, and sy[] bandwidth is never wasted on
//  rejected candidates.
//
//  This mirrors the SoA layout in TrigTrackSeedGenerator, where
//  the φ coordinate is stored separately so the η window check
//  can reject hits before the full 3D calculation.
//
//  Note: sx[] and sy[] are still allocated fresh each call.
//  V3 removes that cost.
// ============================================================

__attribute__((noinline))
static long long nearestNeighbour_SoA(const std::vector<Point>& pts) {
    const int N = static_cast<int>(pts.size());

    // Split AoS into two contiguous arrays — heap allocated here
    std::vector<float> sx(N), sy(N);
    for (int i = 0; i < N; ++i) { sx[i] = pts[i].x; sy[i] = pts[i].y; }

    long long found = 0;

    for (int i = 0; i < N; ++i) {
        float best_dist2 = std::numeric_limits<float>::max();
        int   best_j     = -1;
        float xi = sx[i], yi = sy[i];

        for (int j = 0; j < N; ++j) {
            if (j == i) continue;
            float dx = sx[j] - xi;
            if (dx * dx >= best_dist2) continue;  // ← skip sy[j] fetch entirely
            float dy = sy[j] - yi;
            float d2 = dx*dx + dy*dy;
            if (d2 < best_dist2) { best_dist2 = d2; best_j = j; }
        }
        if (best_j >= 0) ++found;
    }
    return found;
}

// ============================================================
//  Version 3 — Intuitive O(N²), SoA + pre-allocated buffers
//
//  Same O(N²) algorithm and SoA layout as V2. One further
//  change: the working buffers sx[] and sy[] are allocated once
//  outside the hot path and reused across events.
//
//  Why it helps:
//  Every call to V2 invokes operator new twice (for sx, sy).
//  On a busy HLT node handling thousands of RoIs per second,
//  that allocator pressure accumulates. Pre-allocating eliminates
//  it entirely — the only work on entry is a linear copy from
//  the input AoS, which is both predictable and cache-warm.
//
//  This mirrors convertedSpacePoints.reserve(5000) in
//  TrigFastTrackFinder::findTracks(), which avoids repeated
//  reallocation across the event loop.
//
//  After V3 we have applied every reasonable memory-level trick.
//  The timing gap between V1 and V3 represents the maximum benefit
//  available from engineering alone. The bottleneck is now purely
//  algorithmic: we are still doing O(N²) work.
// ============================================================

struct NNWorkspace {
    std::vector<float> sx, sy;
    explicit NNWorkspace(int N) : sx(N), sy(N) {}
};

__attribute__((noinline))
static long long nearestNeighbour_SoA_PreAlloc(const std::vector<Point>& pts,
                                                NNWorkspace& ws) {
    const int N = static_cast<int>(pts.size());

    // Reuse pre-allocated buffers — no heap allocation on the hot path
    for (int i = 0; i < N; ++i) { ws.sx[i] = pts[i].x; ws.sy[i] = pts[i].y; }

    long long found = 0;

    for (int i = 0; i < N; ++i) {
        float best_dist2 = std::numeric_limits<float>::max();
        int   best_j     = -1;
        float xi = ws.sx[i], yi = ws.sy[i];

        for (int j = 0; j < N; ++j) {
            if (j == i) continue;
            float dx = ws.sx[j] - xi;
            if (dx * dx >= best_dist2) continue;
            float dy = ws.sy[j] - yi;
            float d2 = dx*dx + dy*dy;
            if (d2 < best_dist2) { best_dist2 = d2; best_j = j; }
        }
        if (best_j >= 0) ++found;
    }
    return found;
}

// ============================================================
//  Version 4 — Sort + window cut, SoA, pre-allocated buffers
//
//  Now the algorithmic insight: sort all points by x, then for
//  each query point scan outward in sorted order and break the
//  moment |Δx| ≥ best_dist. Every point beyond that cut is
//  provably farther away — no full distance calculation needed.
//
//  This is exactly the RoI η-φ window in TrigFastTrackFinder:
//  hits outside the window cannot form a valid seed, so they
//  are never evaluated.
//
//  Complexity: O(N log N + N·k)  where k ≪ N in practice.
//  For N = 5000 uniformly distributed points, the window cut
//  eliminates ~98% of inner-loop iterations before a single
//  multiply is wasted on them.
//
//  The SoA layout and pre-allocated buffers from V3 are
//  retained — they still contribute, but are now a minor
//  fraction of the total gain.
// ============================================================

struct NNSearchState {
    std::vector<int>   idx;    // sorted index array
    std::vector<float> sx, sy; // coordinates in sorted-x order
    std::vector<int>   rank;   // original index → position in sorted order

    explicit NNSearchState(int N) : idx(N), sx(N), sy(N), rank(N) {}
};

__attribute__((noinline))
static long long nearestNeighbour_Sorted_SoA(const std::vector<Point>& pts,
                                              NNSearchState& state) {
    const int N = static_cast<int>(pts.size());

    // Sort by x — O(N log N), run once per event
    std::iota(state.idx.begin(), state.idx.end(), 0);
    std::sort(state.idx.begin(), state.idx.end(), [&](int a, int b) {
        return pts[a].x < pts[b].x;
    });

    // Fill SoA buffers in sorted order and build rank lookup
    for (int s = 0; s < N; ++s) {
        state.sx[s] = pts[state.idx[s]].x;
        state.sy[s] = pts[state.idx[s]].y;
    }
    for (int s = 0; s < N; ++s) state.rank[state.idx[s]] = s;

    long long found = 0;

    for (int i = 0; i < N; ++i) {
        float best_dist2 = std::numeric_limits<float>::max();
        int   best_j     = -1;

        int   si = state.rank[i];
        float xi = state.sx[si];
        float yi = state.sy[si];

        // Scan right in x — break as soon as |Δx|² ≥ best_dist²
        for (int s = si + 1; s < N; ++s) {
            float dx = state.sx[s] - xi;
            if (dx * dx >= best_dist2) break;   // ← the window cut: the key line
            float dy = state.sy[s] - yi;
            float d2 = dx*dx + dy*dy;
            if (d2 < best_dist2) { best_dist2 = d2; best_j = state.idx[s]; }
        }

        // Scan left in x — same window cut in the opposite direction
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
//
//  Events are pre-generated with varying seeds so the compiler
//  cannot memoise results across iterations. All four versions
//  process the same event data, preserving the correctness check.
// ============================================================

template<typename Fn>
static BenchResult benchmark(const std::string& label,
                              Fn fn,
                              const std::vector<std::vector<Point>>& events,
                              double baseline_ms = 0.0) {
    const int n_events = static_cast<int>(events.size());

    fn(events[0]);  // warmup

    auto t0 = std::chrono::high_resolution_clock::now();
    long long result = 0;
    for (int e = 0; e < n_events; ++e) result += fn(events[e]);
    auto t1 = std::chrono::high_resolution_clock::now();

    double total_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    double ms_each  = total_ms / n_events;

    BenchResult r;
    r.label          = label;
    r.ms_per_event   = ms_each;
    r.throughput_kHz = 1.0 / ms_each;
    r.n_pairs_found  = result / n_events;
    r.speedup        = (baseline_ms > 0.0) ? baseline_ms / ms_each : 1.0;
    return r;
}

// ============================================================
//  Diagnostic: inner-loop work analysis
// ============================================================

static void printWorkAnalysis(const std::vector<Point>& pts) {
    const int N = static_cast<int>(pts.size());

    long long intuitive_iters  = static_cast<long long>(N) * (N - 1);
    long long sorted_iters = 0;

    std::vector<int>   idx(N);
    std::vector<float> sx(N), sy(N);
    std::vector<int>   rank(N);

    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b){ return pts[a].x < pts[b].x; });
    for (int s = 0; s < N; ++s) { sx[s] = pts[idx[s]].x; sy[s] = pts[idx[s]].y; }
    for (int s = 0; s < N; ++s) rank[idx[s]] = s;

    for (int i = 0; i < N; ++i) {
        float best_dist2 = std::numeric_limits<float>::max();
        int   si = rank[i];
        float xi = sx[si], yi = sy[si];
        for (int s = si + 1; s < N; ++s) {
            ++sorted_iters;
            float dx = sx[s] - xi; if (dx*dx >= best_dist2) break;
            float dy = sy[s] - yi;
            float d2 = dx*dx + dy*dy;
            if (d2 < best_dist2) best_dist2 = d2;
        }
        for (int s = si - 1; s >= 0; --s) {
            ++sorted_iters;
            float dx = xi - sx[s]; if (dx*dx >= best_dist2) break;
            float dy = yi - sy[s];
            float d2 = dx*dx + dy*dy;
            if (d2 < best_dist2) best_dist2 = d2;
        }
    }

    std::cout << "\n  Work analysis (N = " << N << " points):\n";
    std::cout << "  ┌─────────────────────────────────┬──────────────┬────────────┐\n";
    std::cout << "  │ Version                         │ Inner iters  │ Reduction  │\n";
    std::cout << "  ├─────────────────────────────────┼──────────────┼────────────┤\n";
    std::cout << "  │ V1–V3: Intuitive O(N²)          │ "
              << std::setw(12) << intuitive_iters << " │     1.0×   │\n";
    std::cout << "  │ V4:    Sort + window cut        │ "
              << std::setw(12) << sorted_iters << " │ "
              << std::setw(6) << std::fixed << std::setprecision(1)
              << (double)intuitive_iters / sorted_iters << "×   │\n";
    std::cout << "  └─────────────────────────────────┴──────────────┴────────────┘\n";
    std::cout << "\n  → The window cut eliminates ~"
              << std::setprecision(0)
              << 100.0 * (1.0 - (double)sorted_iters / intuitive_iters)
              << "% of distance calculations.\n"
              << "    V1–V3 differ only in memory behaviour — identical arithmetic work.\n"
              << "    (Window cut = RoI η-φ gate in TrigFastTrackFinder.)\n\n";
}

// ============================================================
//  Main
// ============================================================

int main() {
    constexpr int N_POINTS = 5000;   // realistic HLT RoI spacepoint count
    constexpr int N_EVENTS = 200;    // events to time over

    printHeader();

    std::cout << "  Configuration:\n";
    std::cout << "    Points per event : " << N_POINTS << "\n";
    std::cout << "    Events timed     : " << N_EVENTS << "\n";
    std::cout << "    HLT rate target  : ~1 kHz sustained per core\n\n";

    // Pre-generate N_EVENTS datasets with distinct seeds.
    // All four versions process the same sequence of events so
    // the correctness check and timing are directly comparable.
    std::vector<std::vector<Point>> events(N_EVENTS);
    for (int e = 0; e < N_EVENTS; ++e)
        events[e] = generatePoints(N_POINTS, 1000.f, 42 + e);

    // ── Correctness check ─────────────────────────────────────
    {
        NNWorkspace   ws(N_POINTS);
        NNSearchState st(N_POINTS);
        auto r1 = nearestNeighbour_Intuitive(events[0]);
        auto r2 = nearestNeighbour_SoA(events[0]);
        auto r3 = nearestNeighbour_SoA_PreAlloc(events[0], ws);
        auto r4 = nearestNeighbour_Sorted_SoA(events[0], st);
        assert(r1 == r2 && r2 == r3 && r3 == r4);
        std::cout << "  ✓ All versions agree on " << r1
                  << " nearest-neighbour pairs found.\n";
    }

    // ── Work analysis (on the first event) ───────────────────
    printWorkAnalysis(events[0]);

    // ── Timing benchmarks ─────────────────────────────────────
    printTableHeader();

    NNWorkspace   ws(N_POINTS);
    NNSearchState st(N_POINTS);

    auto r1 = benchmark("V1: Intuitive O(N²), AoS",
        [&](const std::vector<Point>& pts){ return nearestNeighbour_Intuitive(pts); },
        events);

    auto r2 = benchmark("V2: Intuitive O(N²), SoA            [+cache]",
        [&](const std::vector<Point>& pts){ return nearestNeighbour_SoA(pts); },
        events, r1.ms_per_event);

    auto r3 = benchmark("V3: Intuitive O(N²), SoA+pre-alloc  [+no heap]",
        [&](const std::vector<Point>& pts){ return nearestNeighbour_SoA_PreAlloc(pts, ws); },
        events, r1.ms_per_event);

    auto r4 = benchmark("V4: Sort+window, SoA+pre-alloc  [+O(Nk)]",
        [&](const std::vector<Point>& pts){ return nearestNeighbour_Sorted_SoA(pts, st); },
        events, r1.ms_per_event);

    printResult(r1);
    printResult(r2);
    printResult(r3);
    printResult(r4);
    printSeparator();

    // ── Interpretation ────────────────────────────────────────
    double memory_speedup = r1.ms_per_event / r3.ms_per_event;
    double algo_speedup   = r3.ms_per_event / r4.ms_per_event;

    std::cout << "\n  Interpretation:\n\n";

    std::cout << "  V1 → V2  (+" << std::setprecision(1) << r2.speedup
              << "×): Memory layout — SoA separates x[] and y[] so the inner\n"
              << "           loop can skip the sy[] fetch when dx² ≥ best_dist².\n"
              << "           Zero change to algorithmic complexity.\n"
              << "           Maps to TrigSiSpacePointBase SoA layout.\n\n";

    std::cout << "  V2 → V3  (+" << r3.speedup
              << "×): Allocation — pre-allocated buffers remove operator new\n"
              << "           from the hot path. Marginal but free.\n"
              << "           Maps to convertedSpacePoints.reserve(5000) in\n"
              << "           TrigFastTrackFinder::findTracks().\n\n";

    std::cout << "  V3 → V4  (+" << r4.speedup
              << "×): Algorithm — sort by x + window cut reduces O(N²) inner\n"
              << "           iterations to O(N·k), k ≪ N. This single change\n"
              << "           dwarfs all memory optimisations combined.\n"
              << "           Maps to RoI-gated seeding in TrigFastTrackFinder.\n\n";

    std::cout << "  Summary:\n";
    std::cout << "    Memory gains   (V1 → V3): +" << std::setprecision(1)
              << memory_speedup << "×\n";
    std::cout << "    Algorithm gain (V3 → V4): +" << algo_speedup << "×\n";
    std::cout << "    Total          (V1 → V4): +" << r4.speedup << "×"
              << "   →   " << r4.throughput_kHz << " kHz"
              << "   (target: ~1 kHz/core)\n\n";

    std::cout << "  Takeaway: good memory practice is necessary but not sufficient.\n"
              << "  Algorithmic complexity sets the ceiling that everything else\n"
              << "  works under. Fix the algorithm first; then tune the constants.\n\n";

    return 0;
}
