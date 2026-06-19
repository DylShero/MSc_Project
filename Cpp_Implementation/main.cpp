// main.cpp
// Command-line options for baseline model.

//   selftest                 verify bit-exact RNG vectors, ICDF, and MC->BS
//   price    [opts]          price one option, print estimate vs Black-Scholes
//   converge [opts]          sweep N, show error shrinking ~ 1/sqrt(N)
//   dump     [opts]          write raw uint32 words (binary) to stdout

// Common opts: --rng {lfsr|xorshift|taus|mt}  --paths N  --seed S
//              --S0 --K --r --sigma --T --put  --count N (dump)

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <string>
#include <vector>
#include <cmath>
#include <unistd.h>     // isatty
#include "rng.hpp"
#include "gaussian.hpp"
#include "blackscholes.hpp"
#include "montecarlo.hpp"

using namespace baseline;

// argument parsing
struct Args {
    std::string rng = "mt";
    uint64_t paths = 10000000ULL;
    uint64_t count = 100000000ULL;   // for dump
    uint64_t seed = 12345ULL;
    OptionSpec opt;
};

static double dbl(const char* s) { return std::strtod(s, nullptr); }
static uint64_t u64(const char* s) { return std::strtoull(s, nullptr, 10); }

static Args parse(int argc, char** argv, int start) {
    Args a;
    for (int i = start; i < argc; ++i) {
        std::string k = argv[i];
        auto need = [&](const char* opt) { return k == opt && i + 1 < argc; };
        if      (need("--rng"))    a.rng = argv[++i];
        else if (need("--paths"))  a.paths = u64(argv[++i]);
        else if (need("--count"))  a.count = u64(argv[++i]);
        else if (need("--seed"))   a.seed = u64(argv[++i]);
        else if (need("--S0"))     a.opt.S0 = dbl(argv[++i]);
        else if (need("--K"))      a.opt.K = dbl(argv[++i]);
        else if (need("--r"))      a.opt.r = dbl(argv[++i]);
        else if (need("--sigma"))  a.opt.sigma = dbl(argv[++i]);
        else if (need("--T"))      a.opt.T = dbl(argv[++i]);
        else if (k == "--put")     a.opt.call = false;
        else if (k == "--call")    a.opt.call = true;
    }
    return a;
}

// test
static int test() {
    int fails = 0;
    auto check = [&](bool ok, const char* what) {
        std::printf("  [%s] %s\n", ok ? "PASS" : "FAIL", what);
        if (!ok) ++fails;
    };
    std::printf("RNG bit-exact reference checks\n");

    //
    {
        uint32_t key[4] = {0x123u, 0x234u, 0x345u, 0x456u};
        MT19937 mt(key, 4);
        const uint32_t expect[5] = {1067595299u, 955945823u, 477289528u,
                                    4107218783u, 4228976476u};
        bool ok = true;
        for (int i = 0; i < 5; ++i) { uint32_t v = mt.next_u32(); ok = ok && (v == expect[i]); }
        check(ok, "MT19937 matches canonical mt19937ar.c vector");
    }
    // The other three are defined by their reference recurrences; verify, determinism, non-degeneracy, and record their first words.
    auto first_words = [](Rng32& g, uint32_t out[4]) {
        for (int i = 0; i < 4; ++i) out[i] = g.next_u32();
    };
    {
        LFSR32 a, b; uint32_t wa[4], wb[4];
        first_words(a, wa); first_words(b, wb);
        bool ok = std::memcmp(wa, wb, sizeof wa) == 0 && (wa[0] | wa[1] | wa[2] | wa[3]);
        check(ok, "LFSR32 deterministic and non-zero");
        std::printf("        first words: %08x %08x %08x %08x\n", wa[0],wa[1],wa[2],wa[3]);
    }
    {
        Xorshift128 a, b; uint32_t wa[4], wb[4];
        first_words(a, wa); first_words(b, wb);
        check(std::memcmp(wa, wb, sizeof wa) == 0, "Xorshift128 deterministic");
        std::printf("        first words: %08x %08x %08x %08x\n", wa[0],wa[1],wa[2],wa[3]);
    }
    {
        Taus88 a, b; uint32_t wa[4], wb[4];
        first_words(a, wa); first_words(b, wb);
        check(std::memcmp(wa, wb, sizeof wa) == 0, "Taus88 deterministic");
        std::printf("        first words: %08x %08x %08x %08x\n", wa[0],wa[1],wa[2],wa[3]);
    }

    std::printf("\nGaussian transform (Box-Muller) checks\n");
    {
        // u1=0.5, u2=0  ->  r = sqrt(-2 ln 0.5), theta = 0  ->  z0 = r, z1 = 0
        double z0, z1;
        box_muller(0.5, 0.0, z0, z1);
        check(std::fabs(z0 - std::sqrt(-2.0 * std::log(0.5))) < 1e-12 &&
              std::fabs(z1) < 1e-12, "box_muller(0.5, 0) == (sqrt(2 ln2), 0)");
    }
    {
        // distribution: sample mean ~ 0, variance ~ 1 over a large batch
        auto rng = make_rng("mt", 4242);
        GaussianBoxMuller g;
        const uint64_t M = 2000000;
        double s = 0.0, s2 = 0.0;
        for (uint64_t i = 0; i < M; ++i) {
            double z = g.next([&]{ return rng->next_unit(); });
            s += z; s2 += z * z;
        }
        double mean = s / M, var = s2 / M - mean * mean;
        std::printf("        sample mean=%.5f  var=%.5f  (N=%llu)\n",
                    mean, var, (unsigned long long)M);
        check(std::fabs(mean) < 0.01 && std::fabs(var - 1.0) < 0.01,
              "Box-Muller: mean ~ 0, variance ~ 1");
    }
    

    std::printf("\nMonte Carlo converges to Black-Scholes\n");
    {
        OptionSpec o; // default ATM call
        double bs = black_scholes(o);
        auto rng = make_rng("mt", 777);
        McResult mc = price_european(*rng, o, 4000000ULL);
        double err_in_se = std::fabs(mc.price - bs) / mc.std_error;
        std::printf("        BS=%.6f  MC=%.6f  SE=%.6f  |err|=%.3f SE\n", bs, mc.price, mc.std_error, err_in_se);
        check(err_in_se < 4.0, "MC within 4 standard errors of BS");
    }

    std::printf("\n%s\n", fails ? "TEST FAILED" : "ALL TESTS PASSED");
    return fails ? 1 : 0;
}

// price 
static int do_price(const Args& a) {
    auto rng = make_rng(a.rng, a.seed);
    double bs = black_scholes(a.opt);
    McResult mc = price_european(*rng, a.opt, a.paths);
    double err = std::fabs(mc.price - bs);
    std::printf("rng=%-11s paths=%llu  %s\n", a.rng.c_str(),
                (unsigned long long)a.paths, a.opt.call ? "CALL" : "PUT");
    std::printf("  Black-Scholes : %.6f\n", bs);
    std::printf("  Monte Carlo   : %.6f  (SE %.6f)\n", mc.price, mc.std_error);
    std::printf("  95%% CI        : [%.6f, %.6f]\n", mc.ci_lo, mc.ci_hi);
    std::printf("  |MC - BS|     : %.6f  (%.2f SE)\n", err, err / mc.std_error);
    return 0;
}

// convergence
static int do_converge(const Args& a) {
    double bs = black_scholes(a.opt);
    std::printf("rng=%s  BS=%.6f\n", a.rng.c_str(), bs);
    std::printf("%12s %12s %12s %12s\n", "paths", "MC", "abs_err", "err*sqrt(N)");
    for (uint64_t N = 1000; N <= a.paths; N *= 10) {
        auto rng = make_rng(a.rng, a.seed);
        McResult mc = price_european(*rng, a.opt, N);
        double err = std::fabs(mc.price - bs);
        std::printf("%12llu %12.6f %12.6f %12.4f\n",
                    (unsigned long long)N, mc.price, err, err * std::sqrt((double)N));
    }
    return 0;
}

// dump 
// Raw little-endian uint32 stream to match FPGA output eg.  ./baseline dump --rng taus | RNG_test stdin32
static int do_dump(const Args& a) {
    if (isatty(fileno(stdout))) {
        std::fprintf(stderr, "refusing to write binary to a terminal; pipe or redirect stdout\n");
        return 2;
    }
    auto rng = make_rng(a.rng, a.seed);
    const size_t BUF = 1 << 16;
    std::vector<uint32_t> buf(BUF);
    uint64_t remaining = a.count;
    while (remaining) {
        size_t n = (remaining < BUF) ? (size_t)remaining : BUF;
        for (size_t i = 0; i < n; ++i) buf[i] = rng->next_u32();
        std::fwrite(buf.data(), sizeof(uint32_t), n, stdout);
        remaining -= n;
    }
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr,
            "usage: %s {selftest|price|converge|dump} [opts]\n"
            "  --rng {lfsr|xorshift|taus|mt}  --paths N  --seed S\n"
            "  --S0 --K --r --sigma --T  --call|--put  --count N(dump)\n", argv[0]);
        return 2;
    }
    std::string cmd = argv[1];
    Args a = parse(argc, argv, 2);
    if (cmd == "selftest") return selftest();
    if (cmd == "price")    return do_price(a);
    if (cmd == "converge") return do_converge(a);
    if (cmd == "dump")     return do_dump(a);
    std::fprintf(stderr, "unknown command: %s\n", cmd.c_str());
    return 2;
}
