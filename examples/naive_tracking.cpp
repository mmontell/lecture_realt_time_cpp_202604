/*
 * naive_tracking.cpp
 *
 * More realistic example: simplified 2D track reconstruction
 *
 * Models a simplified ATLAS ITk-like barrel tracker in the
 * transverse (r-φ) plane:
 *   - 9 concentric cylindrical silicon layers (5 pixel + 4 strip)
 *   - Charged particles curve in a 4 Tesla solenoidal B-field
 *   - Tracks are circles in the transverse plane
 *
 * Physics recap:
 *   In a uniform B-field along the z axis
 *   (the direction orthogonal to the (r-φ) plane)
 *   a charged particle follows
 *   a helix. Projected into the transverse (x-y) plane, this
 *   is a circle with radius:
 *
 *       R [m] = pT [GeV] / (0.3 × B [T])
 *
 *   For B = 4T:  R = pT / 1.2  (in metres)
 *
 * The code is written intuitively, without specific care 
 * for optimization.
 * Comments marked with ⚠ PERFORMANCE highlight places where
 * the code could be made more efficient, becoming more suited to 
 * a real time usage.
 *
 * Compile:
 *   g++ -O2 -std=c++17 -o naive_tracking naive_tracking.cpp
 *
 * Run:
 *   ./naive_tracking
 *
 * Expected runtime: O(10s) before optimizations
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <map>
#include <string>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <memory>
#include <cassert>

// ============================================================
//  Constants
// ============================================================

constexpr double B_FIELD  = 4.0;      // Tesla — ATLAS solenoid
constexpr double K_FACTOR = 0.3;      // pT = K × B × R  (natural units: GeV, T, m)
constexpr double PI       = 3.14159;
constexpr double TWO_PI   = 2.0 * PI;

// ITk-inspired barrel layer radii [mm]
// 5 pixel layers (inner) + 4 strip layers (outer)
// Approximate values from the ITk TDR layout
constexpr int    N_LAYERS = 9;
constexpr double LAYER_RADII[N_LAYERS] = {
    34.0,   // Pixel L0 (innermost)
    70.0,   // Pixel L1
    116.0,  // Pixel L2
    172.0,  // Pixel L3
    230.0,  // Pixel L4 (outermost pixel)
    405.0,  // Strip L0
    562.0,  // Strip L1
    762.0,  // Strip L2
    1000.0  // Strip L3 (outermost)
};

// Hit resolution per layer [mm] — pixel vs strip
// Using position resolution in the bending plane
constexpr double HIT_SIGMA[N_LAYERS] = {
    0.010,  // Pixel — 10 µm
    0.010,
    0.010,
    0.010,
    0.010,
    0.080,  // Strip — 80 µm
    0.080,
    0.080,
    0.080
};

// Effective sigma for chi2 calculation [mm]
// Includes alignment, multiple scattering, etc.
constexpr double FIT_SIGMA[N_LAYERS] = {
    0.10,   // Pixel — inflated to ~100 µm effective
    0.10,
    0.10,
    0.10,
    0.10,
    0.50,   // Strip — inflated to ~500 µm effective
    0.50,
    0.50,
    0.50
};


// ============================================================
//  Hit class — deliberately uses virtual functions
//
//  ⚠ PERFORMANCE: virtual dispatch in a base class that will
//  be called millions of times inside the seeding/fitting loop.
//  In real-time code, templates or concrete types are preferred.
// ============================================================

class IHit {
public:
    virtual ~IHit() = default;
    virtual double x() const = 0;
    virtual double y() const = 0;
    virtual double r() const = 0;
    virtual double phi() const = 0;
    virtual int    layer() const = 0;
    virtual double sigma() const = 0;
};

// ⚠ PERFORMANCE: each PixelHit / StripHit is individually
// heap-allocated via std::make_unique later. This means
// every hit is at a random address — pointer-chasing when
// iterating over the hit collection.

class PixelHit : public IHit {
    double m_x, m_y, m_r, m_phi;
    int    m_layer;
    double m_sigma;
public:
    PixelHit(double x, double y, int layer, double sigma)
        : m_x(x), m_y(y),
          m_r(std::sqrt(x*x + y*y)),
          m_phi(std::atan2(y, x)),
          m_layer(layer), m_sigma(sigma) {}

    double x()     const override { return m_x; }
    double y()     const override { return m_y; }
    double r()     const override { return m_r; }
    double phi()   const override { return m_phi; }
    int    layer() const override { return m_layer; }
    double sigma() const override { return m_sigma; }
};

class StripHit : public IHit {
    double m_x, m_y, m_r, m_phi;
    int    m_layer;
    double m_sigma;
public:
    StripHit(double x, double y, int layer, double sigma)
        : m_x(x), m_y(y),
          m_r(std::sqrt(x*x + y*y)),
          m_phi(std::atan2(y, x)),
          m_layer(layer), m_sigma(sigma) {}

    double x()     const override { return m_x; }
    double y()     const override { return m_y; }
    double r()     const override { return m_r; }
    double phi()   const override { return m_phi; }
    int    layer() const override { return m_layer; }
    double sigma() const override { return m_sigma; }
};


// ============================================================
//  Track candidate — result of fitting
// ============================================================

struct TrackCandidate {
    double pT;          // GeV
    double phi0;        // direction at closest approach
    double d0;          // transverse impact parameter [mm]
                        // i.e. minimum track approach distance to origin
                        // usually 0, but some physics processes may result 
                        // in large d0 tracks -- Displace tracks!
                        // you will learn about those if you continue in particle physics!
    double chi2;        // fit chi-square
    int    nHits;       // number of hits on track
    int    charge;      // +1 or -1 
    std::vector<const IHit*> hits;  // pointers to associated hits
};


// ============================================================
//  Event generation — simulate tracks + noise hits
// ============================================================

struct TruthTrack {
    double pT, phi0, d0, charge;
};

struct SimulatedEvent {
    std::vector<std::unique_ptr<IHit>> hits;

    // ⚠ PERFORMANCE: using std::map for layer-indexed access.
    // Each map node is a separate heap allocation with pointers.
    // A flat std::array<std::vector<IHit*>, N_LAYERS> would be
    // far more cache-friendly.
    std::map<int, std::vector<const IHit*>> hitsByLayer;

    std::vector<TruthTrack> truth;
};

// Generate a circle (track) and produce hits at each layer crossing
static void generateTrackHits(
    double pT, double phi0, double d0, double charge,
    std::mt19937& rng, SimulatedEvent& event)
{
    // Circle radius in mm
    double R = (pT / (K_FACTOR * B_FIELD)) * 1000.0;  // GeV → mm

    // Centre of curvature
    // For a track with d0 ≈ 0, the centre is at distance R
    // perpendicular to the initial direction
    double xc = (d0 + charge * R) * std::sin(phi0);
    double yc = -(d0 + charge * R) * std::cos(phi0);

    for (int iL = 0; iL < N_LAYERS; ++iL) {
        double rL = LAYER_RADII[iL];

        // Does the circle cross this layer?
        double dc = std::sqrt(xc*xc + yc*yc);
        if (std::abs(dc - R) > rL || dc + R < rL) continue;

        // Solve circle–circle intersection
        // Layer circle: x² + y² = rL²
        // Track circle: (x - xc)² + (y - yc)² = R²
        double cos_alpha = (dc*dc + rL*rL - R*R) / (2.0 * dc * rL);
        if (std::abs(cos_alpha) > 1.0) continue;

        double phi_centre = std::atan2(yc, xc);
        double alpha = std::acos(cos_alpha);

        // Take the intersection closest to the initial direction
        double phi_hit = phi_centre + alpha;  // could also be - alpha

        // Smear by resolution
        std::normal_distribution<double> smear(0.0, HIT_SIGMA[iL]);
        double dphi = smear(rng) / rL;  // angular smearing
        phi_hit += dphi;

        double hx = rL * std::cos(phi_hit);
        double hy = rL * std::sin(phi_hit);

        // ⚠ PERFORMANCE: each hit is individually heap-allocated.
        // In a hot path this means one malloc() per hit per event.
        // A pre-allocated pool or flat buffer would eliminate this.
        std::unique_ptr<IHit> hit;
        if (iL < 5) {
            hit = std::make_unique<PixelHit>(hx, hy, iL, HIT_SIGMA[iL]);
        } else {
            hit = std::make_unique<StripHit>(hx, hy, iL, HIT_SIGMA[iL]);
        }

        event.hitsByLayer[iL].push_back(hit.get());
        event.hits.push_back(std::move(hit));
    }
}

static SimulatedEvent generateEvent(int nTracks, int nNoise, unsigned seed) {
    std::mt19937 rng(seed);

    SimulatedEvent event;

    // pT spectrum: exponential fall-off, minimum 0.5 GeV
    std::exponential_distribution<double> ptDist(2.0);  // mean ~ 0.5 GeV
    std::uniform_real_distribution<double> phiDist(0, TWO_PI);
    std::normal_distribution<double> d0Dist(0.0, 0.1);  // mm
    std::bernoulli_distribution chargeDist(0.5);

    for (int i = 0; i < nTracks; ++i) {
        double pT  = 0.5 + ptDist(rng);       // min 0.5 GeV
        double phi = phiDist(rng);
        double d0  = d0Dist(rng);
        double q   = chargeDist(rng) ? 1.0 : -1.0;

        event.truth.push_back({pT, phi, d0, q});
        generateTrackHits(pT, phi, d0, q, rng, event);
    }

    // Add noise hits — random positions on random layers
    std::uniform_int_distribution<int> layerDist(0, N_LAYERS - 1);
    for (int i = 0; i < nNoise; ++i) {
        int iL = layerDist(rng);
        double rL = LAYER_RADII[iL];
        double phi = phiDist(rng);
        double hx = rL * std::cos(phi);
        double hy = rL * std::sin(phi);

        std::unique_ptr<IHit> hit;
        if (iL < 5) {
            hit = std::make_unique<PixelHit>(hx, hy, iL, HIT_SIGMA[iL]);
        } else {
            hit = std::make_unique<StripHit>(hx, hy, iL, HIT_SIGMA[iL]);
        }
        event.hitsByLayer[iL].push_back(hit.get());
        event.hits.push_back(std::move(hit));
    }

    return event;
}


// ============================================================
//  Track seeding — triplet finder (NAIVE)
//
//  For each combination of hits on three different layers,
//  check if they are compatible with a circle (helix in 2D).
//
//  ⚠ PERFORMANCE: this is O(N³) in the number of hits.
//  No phi-sorting, no windowing, no early rejection.
//  Every triplet is checked even if it's geometrically
//  impossible. This is the combinatorial bottleneck.
// ============================================================

struct Seed {
    const IHit* inner;
    const IHit* middle;
    const IHit* outer;
    double curvature;     // 1/R estimated from the 3 points
    double phi_estimate;
};

// Compute the signed curvature (1/R) of the circle through 3 points
static double circleCurvature(double x1, double y1,
                                double x2, double y2,
                                double x3, double y3) {
    // Area of triangle × 4 / product of side lengths
    double ax = x2 - x1, ay = y2 - y1;
    double bx = x3 - x1, by = y3 - y1;
    double cross = ax * by - ay * bx;

    double a = std::sqrt(ax*ax + ay*ay);
    double b = std::sqrt(bx*bx + by*by);
    double cx = x3 - x2, cy = y3 - y2;
    double c = std::sqrt(cx*cx + cy*cy);

    double denom = a * b * c;
    if (denom < 1e-12) return 0.0;
    return 2.0 * cross / denom;  // signed: +ve = clockwise
}

// ⚠ PERFORMANCE: returns a std::vector of Seed by value —
// each call may trigger reallocation.
// ⚠ PERFORMANCE: uses std::list internally for "flexibility"
// but this means scattered heap nodes during iteration.

static std::vector<Seed> findSeeds(const SimulatedEvent& event,
                                    double maxCurvature = 0.01,   // 1/mm
                                    double phiWindow    = 0.3)    // rad
{
    // ⚠ PERFORMANCE: using std::list to collect seeds,
    // then copying to a vector at the end.
    // This is how a student might "keep things flexible."
    std::list<Seed> seedList;

    // Use layers 0, 2, 4 as the seed triplet layers
    // (inner pixel, middle pixel, outer pixel)
    const int seedLayers[3] = {0, 2, 4};

    // ⚠ PERFORMANCE: hitsByLayer is a std::map — each [] access
    // involves a tree traversal. A flat array would be O(1).
    auto it0 = event.hitsByLayer.find(seedLayers[0]);
    auto it1 = event.hitsByLayer.find(seedLayers[1]);
    auto it2 = event.hitsByLayer.find(seedLayers[2]);

    if (it0 == event.hitsByLayer.end() ||
        it1 == event.hitsByLayer.end() ||
        it2 == event.hitsByLayer.end()) {
        return {};
    }

    const auto& hitsL0 = it0->second;
    const auto& hitsL1 = it1->second;
    const auto& hitsL2 = it2->second;

    // ⚠ PERFORMANCE: Triple nested loop — O(N₀ × N₁ × N₂).
    // No sorting by phi, no windowing, no early exit.
    // Every single combination is tested.
    for (const IHit* h0 : hitsL0) {
        for (const IHit* h1 : hitsL1) {

            // ⚠ PERFORMANCE: this phi check could be done BEFORE
            // entering the inner loop if hits were sorted by phi.
            double dphi01 = h1->phi() - h0->phi();
            // Wrap to [-π, π]
            while (dphi01 >  PI) dphi01 -= TWO_PI;
            while (dphi01 < -PI) dphi01 += TWO_PI;
            if (std::abs(dphi01) > phiWindow) continue;

            for (const IHit* h2 : hitsL2) {

                double dphi12 = h2->phi() - h1->phi();
                while (dphi12 >  PI) dphi12 -= TWO_PI;
                while (dphi12 < -PI) dphi12 += TWO_PI;
                if (std::abs(dphi12) > phiWindow) continue;

                // ⚠ PERFORMANCE: virtual dispatch ×3 per triplet
                // for x() and y() calls.
                double curv = circleCurvature(
                    h0->x(), h0->y(),
                    h1->x(), h1->y(),
                    h2->x(), h2->y());

                if (std::abs(curv) > maxCurvature) continue;

                double phi_est = std::atan2(
                    h1->y() - h0->y(),
                    h1->x() - h0->x());

                // ⚠ PERFORMANCE: push_back on a std::list —
                // one heap allocation per seed.
                seedList.push_back({h0, h1, h2, curv, phi_est});
            }
        }
    }

    // ⚠ PERFORMANCE: copying from list to vector — O(N) copy
    // that wouldn't be needed if we used a vector from the start.
    std::vector<Seed> seeds(seedList.begin(), seedList.end());
    return seeds;
}


// ============================================================
//  Track fitting — simple circle fit (NAIVE)
//
//  Given a seed, extend it by collecting compatible hits on
//  other layers, then fit a circle to all associated hits.
//  Uses a basic χ² minimisation (algebraic circle fit).
//
//  ⚠ PERFORMANCE: re-scans all hits on every layer for every
//  seed — no caching, no bookkeeping of used hits.
// ============================================================

// Circle through exactly 3 points — numerically stable for small arcs
static std::tuple<double, double, double> circleFrom3Points(
    double x1, double y1, double x2, double y2, double x3, double y3)
{
    double ax = x1, ay = y1;
    double bx = x2, by = y2;
    double cx = x3, cy = y3;

    double D = 2.0 * (ax*(by - cy) + bx*(cy - ay) + cx*(ay - by));
    if (std::abs(D) < 1e-12) return {0, 0, 1e9};

    double ux = ((ax*ax + ay*ay)*(by - cy) + (bx*bx + by*by)*(cy - ay) + (cx*cx + cy*cy)*(ay - by)) / D;
    double uy = ((ax*ax + ay*ay)*(cx - bx) + (bx*bx + by*by)*(ax - cx) + (cx*cx + cy*cy)*(bx - ax)) / D;

    double R = std::sqrt((ax - ux)*(ax - ux) + (ay - uy)*(ay - uy));
    return {ux, uy, R};
}

// Algebraic circle fit: Kåsa method
// Returns (xc, yc, R) of the best-fit circle
static std::tuple<double, double, double> fitCircle(
    const std::vector<const IHit*>& hits)
{
    // ⚠ PERFORMANCE: allocating temporary vectors inside
    // the fit function, which is called once per seed.
    // These should be pre-allocated and reused.
    int n = static_cast<int>(hits.size());
    if (n < 3) return {0, 0, 0};

    // For 3 hits, use exact circle
    if (n == 3) {
        return circleFrom3Points(
            hits[0]->x(), hits[0]->y(),
            hits[1]->x(), hits[1]->y(),
            hits[2]->x(), hits[2]->y());
    }

    // For more hits, use iterative mean of 3-point circles
    // (more stable than Kåsa for small arcs)
    // ⚠ PERFORMANCE: O(n²) approach — could use a proper
    // least-squares fit for production code
    double sum_xc = 0, sum_yc = 0, sum_R = 0;
    int count = 0;

    for (int i = 0; i < n - 2; i += 2) {
        auto [xc, yc, R] = circleFrom3Points(
            hits[i]->x(), hits[i]->y(),
            hits[i+1]->x(), hits[i+1]->y(),
            hits[i+2]->x(), hits[i+2]->y());
        if (R > 0 && R < 1e8) {
            sum_xc += xc; sum_yc += yc; sum_R += R;
            count++;
        }
    }

    if (count == 0) return {0, 0, 1e9};
    return {sum_xc / count, sum_yc / count, sum_R / count};
}

static double computeChi2(const std::vector<const IHit*>& hits,
                           double xc, double yc, double R) {
    double chi2 = 0;
    for (const IHit* h : hits) {
        double dx = h->x() - xc;
        double dy = h->y() - yc;
        double residual = std::sqrt(dx*dx + dy*dy) - R;
        // ⚠ PERFORMANCE: virtual dispatch for h->sigma() and h->layer()
        double sigma = FIT_SIGMA[h->layer()];
        chi2 += (residual * residual) / (sigma * sigma);
    }
    return chi2;
}


// ============================================================
//  Full reconstruction chain — seed → extend → fit
//
//  ⚠ PERFORMANCE: no early rejection, no RoI windowing,
//  no shared-hit removal between tracks.
//  Every seed is fully processed regardless of quality.
// ============================================================

static std::vector<TrackCandidate> reconstructTracks(
    const SimulatedEvent& event,
    double phiMatchWindow = 0.15,    // rad — window for hit extension
    double chi2Cut        = 1000.0)   // χ²/ndf cut (loose — naive code)
{
    // ⚠ PERFORMANCE: no reserve() — vector grows dynamically
    std::vector<TrackCandidate> tracks;

    auto seeds = findSeeds(event);

    for (const auto& seed : seeds) {

        // Start with the 3 seed hits
        // ⚠ PERFORMANCE: allocating a new vector for every seed
        std::vector<const IHit*> trackHits;
        trackHits.push_back(seed.inner);
        trackHits.push_back(seed.middle);
        trackHits.push_back(seed.outer);

        // Extend: look for compatible hits on remaining layers
        // ⚠ PERFORMANCE: iterates over ALL layers and ALL hits
        // on each layer for every seed. No sorting or windowing.
        for (int iL = 0; iL < N_LAYERS; ++iL) {
            // Skip seed layers
            if (iL == 0 || iL == 2 || iL == 4) continue;

            auto layerIt = event.hitsByLayer.find(iL);
            if (layerIt == event.hitsByLayer.end()) continue;

            const IHit* bestHit = nullptr;
            double bestDphi = phiMatchWindow;

            // ⚠ PERFORMANCE: linear scan of all hits on the layer
            for (const IHit* h : layerIt->second) {
                double dphi = h->phi() - seed.phi_estimate;
                while (dphi >  PI) dphi -= TWO_PI;
                while (dphi < -PI) dphi += TWO_PI;

                if (std::abs(dphi) < bestDphi) {
                    bestDphi = std::abs(dphi);
                    bestHit = h;
                }
            }

            if (bestHit) {
                trackHits.push_back(bestHit);
            }
        }

        // Require minimum 5 hits
        if (static_cast<int>(trackHits.size()) < 5) continue;

        // Fit the circle
        auto [xc, yc, R] = fitCircle(trackHits);
        if (R < 1.0) continue;  // degenerate fit

        double chi2 = computeChi2(trackHits, xc, yc, R);
        double ndf  = static_cast<double>(trackHits.size()) - 3;  // 3 params
        if (ndf > 0 && chi2 / ndf > chi2Cut) continue;

        // Convert radius to pT
        double pT = K_FACTOR * B_FIELD * R / 1000.0;  // mm → m → GeV

        // Estimate phi at origin from circle centre
        double phi0 = std::atan2(-xc, yc);

        // d0: signed distance of closest approach to origin
        double dc = std::sqrt(xc*xc + yc*yc);
        double d0 = dc - R;

        // Charge from curvature sign
        int charge = (seed.curvature > 0) ? -1 : +1;

        // ⚠ PERFORMANCE: TrackCandidate contains a vector of
        // hit pointers — each push_back to 'tracks' copies it.
        // std::move or emplace_back would avoid the copy.
        TrackCandidate trk;
        trk.pT     = pT;
        trk.phi0   = phi0;
        trk.d0     = d0;
        trk.chi2   = chi2;
        trk.nHits  = static_cast<int>(trackHits.size());
        trk.charge = charge;
        trk.hits   = trackHits;
        tracks.push_back(trk);
    }

    return tracks;
}


// ============================================================
//  Main — run reconstruction with timing
// ============================================================

int main() {
    constexpr int N_EVENTS  = 20;
    constexpr int N_TRACKS  = 200;    // truth tracks per event
    constexpr int N_NOISE   = 1000;   // noise hits per event

    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║   Naive 2D Track Reconstruction — ITk-like Geometry          ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════╝\n\n";

    std::cout << "  Configuration:\n";
    std::cout << "    Detector layers     : " << N_LAYERS << " (5 pixel + 4 strip)\n";
    std::cout << "    B-field             : " << B_FIELD << " T (along z)\n";
    std::cout << "    Tracks per event    : " << N_TRACKS << "\n";
    std::cout << "    Noise hits per event: " << N_NOISE << "\n";
    std::cout << "    Events              : " << N_EVENTS << "\n\n";

    // ── Timing variables ───────────────────────────────────
    double total_time_ms    = 0;
    double total_seed_ms    = 0;
    double total_fit_ms     = 0;
    int    total_seeds      = 0;
    int    total_tracks     = 0;
    int    total_hits       = 0;

    // ── Per-stage timing (like TrigFTF's Monitored::Timer) ──
    std::cout << "  Running reconstruction...\n\n";

    auto t_total_start = std::chrono::high_resolution_clock::now();

    for (int iEvt = 0; iEvt < N_EVENTS; ++iEvt) {

        // Generate event
        auto t0 = std::chrono::high_resolution_clock::now();
        auto event = generateEvent(N_TRACKS, N_NOISE, 42 + iEvt);
        auto t1 = std::chrono::high_resolution_clock::now();

        total_hits += static_cast<int>(event.hits.size());

        // Seeding
        auto t2 = std::chrono::high_resolution_clock::now();
        auto seeds = findSeeds(event);
        auto t3 = std::chrono::high_resolution_clock::now();

        total_seeds += static_cast<int>(seeds.size());

        // Full reconstruction
        auto t4 = std::chrono::high_resolution_clock::now();
        auto tracks = reconstructTracks(event);
        auto t5 = std::chrono::high_resolution_clock::now();

        total_tracks += static_cast<int>(tracks.size());

        double evt_ms = std::chrono::duration<double, std::milli>(t5 - t0).count();
        double seed_ms = std::chrono::duration<double, std::milli>(t3 - t2).count();
        double fit_ms = std::chrono::duration<double, std::milli>(t5 - t4).count();
        total_time_ms += evt_ms;
        total_seed_ms += seed_ms;
        total_fit_ms  += fit_ms;

        // Print first event's tracks as a sample
        if (iEvt == 0) {
            std::cout << "  ── Event 0 reconstructed tracks ─────────────────────────\n";
            std::cout << "    " << std::left << std::setw(6) << "Track"
                      << std::setw(12) << "pT [GeV]"
                      << std::setw(12) << "φ₀ [rad]"
                      << std::setw(12) << "d₀ [mm]"
                      << std::setw(10) << "χ²/ndf"
                      << std::setw(6)  << "nHits"
                      << std::setw(8)  << "charge"
                      << "\n";
            std::cout << "    ------------------------------------------------------------------\n";

            int nPrint = std::min(static_cast<int>(tracks.size()), 20);
            for (int i = 0; i < nPrint; ++i) {
                const auto& t = tracks[i];
                double ndf = std::max(1.0, static_cast<double>(t.nHits - 3));
                std::cout << "    " << std::left << std::setw(6) << i
                          << std::fixed
                          << std::setw(12) << std::setprecision(3) << t.pT
                          << std::setw(12) << std::setprecision(4) << t.phi0
                          << std::setw(12) << std::setprecision(4) << t.d0
                          << std::setw(10) << std::setprecision(2) << t.chi2 / ndf
                          << std::setw(6)  << t.nHits
                          << std::setw(8)  << (t.charge > 0 ? "+1" : "-1")
                          << "\n";
            }
            if (static_cast<int>(tracks.size()) > 20) {
                std::cout << "    ... (" << tracks.size() - 20 << " more tracks)\n";
            }
            std::cout << "\n";
        }
    }

    auto t_total_end = std::chrono::high_resolution_clock::now();
    double wall_ms = std::chrono::duration<double, std::milli>(t_total_end - t_total_start).count();

    // ── Summary ────────────────────────────────────────────
    std::cout << "  ── Performance summary ──────────────────────────────────\n\n";
    std::cout << "    Total wall time      : " << std::fixed << std::setprecision(1)
              << wall_ms << " ms\n";
    std::cout << "    Time per event       : " << std::setprecision(1)
              << total_time_ms / N_EVENTS << " ms\n";
    std::cout << "    Event rate           : " << std::setprecision(2)
              << N_EVENTS / (wall_ms / 1000.0) << " Hz\n\n";

    std::cout << "    ┌────────────────────┬──────────────┬──────────────┐\n";
    std::cout << "    │ Stage              │ Time/event   │ Fraction     │\n";
    std::cout << "    ├────────────────────┼──────────────┼──────────────┤\n";
    std::cout << "    │ Seeding            │ "
              << std::setw(8) << std::setprecision(1) << total_seed_ms / N_EVENTS << " ms │ "
              << std::setw(8) << std::setprecision(1) << 100.0 * total_seed_ms / total_time_ms << " %   │\n";
    std::cout << "    │ Fit + extension    │ "
              << std::setw(8) << std::setprecision(1) << total_fit_ms / N_EVENTS << " ms │ "
              << std::setw(8) << std::setprecision(1) << 100.0 * total_fit_ms / total_time_ms << " %   │\n";
    std::cout << "    └────────────────────┴──────────────┴──────────────┘\n\n";

    std::cout << "    Avg hits/event   : " << total_hits / N_EVENTS << "\n";
    std::cout << "    Avg seeds/event  : " << total_seeds / N_EVENTS << "\n";
    std::cout << "    Avg tracks/event : " << total_tracks / N_EVENTS << "\n\n";

    std::cout << "    HLT target rate  : ~1 kHz per core\n";
    std::cout << "    Current rate     : " << std::setprecision(2)
              << N_EVENTS / (wall_ms / 1000.0) << " Hz\n";
    std::cout << "    Gap to target    : "
              << std::setprecision(0) << 1000.0 / (N_EVENTS / (wall_ms / 1000.0))
              << "× too slow\n\n";

    // ── Annotated performance bottlenecks ──────────────────
    std::cout << "  ── Where the time goes (bottleneck summary) ───────────\n\n";
    std::cout << "    1. SEEDING: O(N³) triple loop with no phi-sorting\n";
    std::cout << "       → Fix: sort hits by φ, binary-search the window\n";
    std::cout << "       → FTF: TrigTrackSeedGenerator does exactly this\n\n";
    std::cout << "    2. HIT STORAGE: std::map<int, vector<IHit*>>\n";
    std::cout << "       → Fix: std::array<vector<Hit>, N_LAYERS>\n";
    std::cout << "       → FTF: flat arrays indexed by layer ID\n\n";
    std::cout << "    3. HIT CLASS: virtual dispatch on x(), y(), phi()\n";
    std::cout << "       → Fix: concrete struct, no inheritance\n";
    std::cout << "       → FTF: TrigSiSpacePointBase is non-virtual\n\n";
    std::cout << "    4. ALLOCATION: unique_ptr per hit, list for seeds\n";
    std::cout << "       → Fix: pre-allocated flat buffers, reserve()\n";
    std::cout << "       → FTF: convertedSpacePoints.reserve(5000)\n\n";
    std::cout << "    5. NO EARLY REJECTION: every seed fully processed\n";
    std::cout << "       → Fix: cheap curvature/pT cut before fitting\n";
    std::cout << "       → FTF: if(nSPs < minHits) return SUCCESS\n\n";

    return 0;
}
