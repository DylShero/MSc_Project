// rng.hpp
// Baseline-model reference implementations of the four GF(2)-linear generators.

// Every generator's recurrence is written in explicit fixed-width
// integer arithmetic (uint32_t / uint64_t) so that it maps one-to-one onto RTL
// and so the FPGA output can be verified bit exactly against this model.

// All four generators are linear over GF(2): a single LFSR, xorshift128, the
// combined-Tausworthe taus88, and MT19937 (a twisted GFSR). 

#ifndef BASELINE_RNG_HPP
#define BASELINE_RNG_HPP

#include <cstdint>
#include <array>
#include <vector>
#include <memory>
#include <string>
#include <stdexcept>

namespace baseline {


// Common 32-bit generator interface. Mirrors the intended hardware lane port
// set: a generator is seeded once, then emits one 32-bit word per "cycle".
struct Rng32 {
    virtual ~Rng32() = default;
    virtual uint32_t next_u32() = 0;          // one output word
    virtual const char* name() const = 0;
    // Map a 32-bit word to a double in the open interval (0,1).
    // (r + 0.5) / 2^32  keeps the result strictly inside (0,1) and is symmetric,
    // which matters for the inverse-CDF transform (avoids +/-inf at the tails).
    double next_unit() { return (static_cast<double>(next_u32()) + 0.5) * (1.0 / 4294967296.0); }
};


// SplitMix64 — used ONLY to expand a single user seed into well-mixed,
// decorrelated seed words for the linear generators (and, later, to produce
// independent per-lane substream seeds on the ARM before loading them into PL).
// This is the safe way to seed: never seed lane i with i.

class SplitMix64 {
    uint64_t s_;
public:
    explicit SplitMix64(uint64_t seed) : s_(seed) {}
    uint64_t next() {
        uint64_t z = (s_ += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
    uint32_t next32() { return static_cast<uint32_t>(next() >> 32); }
};

// 1. LFSR — 32-bit maximal-length Galois LFSR.
//    Toggle mask 0x80200003 == taps [32,22,2,1] 
class LFSR32 : public Rng32 {
    uint32_t state_;
    static constexpr uint32_t MASK = 0x80200003u; // taps 32,22,2,1
public:
    explicit LFSR32(uint32_t seed = 0xACE1u) { set_state(seed); }
    void set_state(uint32_t s) { state_ = s ? s : 0xDEADBEEFu; }
    uint32_t state() const { return state_; }

    inline void step_bit() {
        uint32_t lsb = state_ & 1u;
        state_ >>= 1;
        state_ ^= (0u - lsb) & MASK;        // branchless conditional XOR
    }
    uint32_t next_u32() override {
        for (int i = 0; i < 32; ++i) step_bit();
        return state_;
    }
    const char* name() const override { return "lfsr"; }
};

// 2. Xorshift128 — Marsaglia (2003), "Xorshift RNGs".
//    Pure shifts and XORs, 32 bits out per call, naturally one word/cycle.
//    Stays strictly GF(2)-linear.
//    Canonical default seeds; state must be non-zero.

class Xorshift128 : public Rng32 {
    uint32_t x_, y_, z_, w_;
public:
    Xorshift128(uint32_t x = 123456789u, uint32_t y = 362436069u, uint32_t z = 521288629u, uint32_t w = 88675123u) : x_(x), y_(y), z_(z), w_(w) {
        if ((x_ | y_ | z_ | w_) == 0u) x_ = 1u;
    }
    uint32_t next_u32() override {
        uint32_t t = x_ ^ (x_ << 11);
        x_ = y_; 
        y_ = z_; 
        z_ = w_;
        w_ = w_ ^ (w_ >> 19) ^ (t ^ (t >> 8));
        return w_;
    }
    const char* name() const override { return "xorshift128"; }
};


// 3. Taus88 — L'Ecuyer. Three LFSR components of periods
//    ~2^31, 2^29, 2^28 XOR-combined to ~2^88. Component seeds must exceed
//    1, 7, 15 respectively, or the components degrade.

class Taus88 : public Rng32 {
    uint32_t z1_, z2_, z3_;
    static uint32_t force(uint32_t v, uint32_t floor_excl) {
        // ensure v > floor_excl
        return (v > floor_excl) ? v : (floor_excl + 1u);
    }
public:
    Taus88(uint32_t s1 = 123456789u, uint32_t s2 = 362436069u, uint32_t s3 = 521288629u) {
        set_state(s1, s2, s3);
    }
    void set_state(uint32_t s1, uint32_t s2, uint32_t s3) {
        z1_ = force(s1, 1u);
        z2_ = force(s2, 7u);
        z3_ = force(s3, 15u);
    }
    uint32_t next_u32() override {
        uint32_t b;
        b   = (((z1_ << 13) ^ z1_) >> 19);
        z1_ = (((z1_ & 4294967294u) << 12) ^ b);
        b   = (((z2_ << 2)  ^ z2_) >> 25);
        z2_ = (((z2_ & 4294967288u) << 4)  ^ b);
        b   = (((z3_ << 3)  ^ z3_) >> 11);
        z3_ = (((z3_ & 4294967280u) << 17) ^ b);
        return z1_ ^ z2_ ^ z3_;
    }
    const char* name() const override { return "taus88"; }
};


// MT19937 — Matsumoto & Nishimura reference .
class MT19937 : public Rng32 {
    static constexpr int      N = 624;
    static constexpr int      M = 397;
    static constexpr uint32_t MATRIX_A   = 0x9908b0dfu;
    static constexpr uint32_t UPPER_MASK = 0x80000000u;
    static constexpr uint32_t LOWER_MASK = 0x7fffffffu;
    std::array<uint32_t, N> mt_{};
    int mti_ = N + 1;
public:
    explicit MT19937(uint32_t seed = 5489u) { init_genrand(seed); }
    MT19937(const uint32_t* key, int key_length) { init_by_array(key, key_length); }

    void init_genrand(uint32_t s) {
        mt_[0] = s;
        for (int i = 1; i < N; ++i)
            mt_[i] = 1812433253u * (mt_[i-1] ^ (mt_[i-1] >> 30)) + static_cast<uint32_t>(i);
        mti_ = N;
    }
    void init_by_array(const uint32_t* key, int key_length) {
        init_genrand(19650218u);
        int i = 1, j = 0;
        int k = (N > key_length) ? N : key_length;
        for (; k; --k) {
            mt_[i] = (mt_[i] ^ ((mt_[i-1] ^ (mt_[i-1] >> 30)) * 1664525u))
                     + key[j] + static_cast<uint32_t>(j);
            ++i; ++j;
            if (i >= N) { mt_[0] = mt_[N-1]; i = 1; }
            if (j >= key_length) j = 0;
        }
        for (k = N - 1; k; --k) {
            mt_[i] = (mt_[i] ^ ((mt_[i-1] ^ (mt_[i-1] >> 30)) * 1566083941u))
                     - static_cast<uint32_t>(i);
            ++i;
            if (i >= N) { mt_[0] = mt_[N-1]; i = 1; }
        }
        mt_[0] = 0x80000000u;
    }
    uint32_t next_u32() override {
        uint32_t y;
        static const uint32_t mag01[2] = { 0x0u, MATRIX_A };
        if (mti_ >= N) {
            int kk;
            for (kk = 0; kk < N - M; ++kk) {
                y = (mt_[kk] & UPPER_MASK) | (mt_[kk+1] & LOWER_MASK);
                mt_[kk] = mt_[kk+M] ^ (y >> 1) ^ mag01[y & 1u];
            }
            for (; kk < N - 1; ++kk) {
                y = (mt_[kk] & UPPER_MASK) | (mt_[kk+1] & LOWER_MASK);
                mt_[kk] = mt_[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 1u];
            }
            y = (mt_[N-1] & UPPER_MASK) | (mt_[0] & LOWER_MASK);
            mt_[N-1] = mt_[M-1] ^ (y >> 1) ^ mag01[y & 1u];
            mti_ = 0;
        }
        y = mt_[mti_++];
        y ^= (y >> 11);
        y ^= (y << 7)  & 0x9d2c5680u;
        y ^= (y << 15) & 0xefc60000u;
        y ^= (y >> 18);
        return y;
    }
    const char* name() const override { return "mt19937"; }
};

// Build RNG function
inline std::unique_ptr<Rng32> make_rng(const std::string& which, uint64_t seed) {
    SplitMix64 sm(seed ? seed : 0xC0FFEEull);
    if (which == "lfsr")        return std::make_unique<LFSR32>(sm.next32() | 1u);
    if (which == "xorshift" || which == "xorshift128")
        return std::make_unique<Xorshift128>(sm.next32(), sm.next32(), sm.next32(), sm.next32() | 1u);
    if (which == "taus" || which == "taus88")
        return std::make_unique<Taus88>(sm.next32(), sm.next32(), sm.next32());
    if (which == "mt" || which == "mt19937")
        return std::make_unique<MT19937>(sm.next32());
    throw std::invalid_argument("unknown rng: " + which);
}

} // namespace baseline
#endif // BASELINE_RNG_HPP
