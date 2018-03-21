#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

#include "constraints.h"

#include <stack>

#define BLOCK_SIZE 64
#define DIGEST_SIZE 32

static inline uint8_t rotr(uint8_t x, uint8_t n) { return (x >> n) | (x << (8 - n)); }
static inline uint8_t maj(uint8_t x, uint8_t y, uint8_t z) { return (x & y) ^ (x & z) ^ (y & z); }
static inline uint8_t sigma0(uint8_t x) { return rotr(x, 2) ^ rotr(x, 3) ^ rotr(x, 5); }
static inline uint8_t sigma1(uint8_t x) { return rotr(x, 1) ^ rotr(x, 4) ^ (x >> 3); }
static inline uint8_t add(uint8_t x, uint8_t y) { return x + y; }

// Differential/probability pair
struct diff_prob
{
    uint8_t diff;   // Output difference
    float   l2prob; // log2 of the probability
};

// Vector of diff/prob pairs -- similar to std::vector
struct diff_vec
{
    struct diff_prob *pairs;    // Distinct samples
    size_t len;                 // Length of pairs
};

struct prop_state
{
    size_t round;               // Which round we are currently propagating
    size_t step;                // Which step in particular
    float l2prob;               // log2(prob) that our trail holds
    uint8_t sched[16];          // Differences in message schedule
    union 
    {
        struct { uint8_t a, b, c, d; }; // Current registers
        uint32_t diff;          // Registers as a 32-bit int
    };
    uint8_t t1, t2, maj;        // Differences in temporary variables
    uint32_t tmp_diffs;         //
    union
    {
        uint8_t trail[4*16];    // Trail of differences
        uint32_t trail32[16];   //
    };
};

// Linear
struct diff_prob propagate_sigma0(uint8_t d_m)
{
    struct diff_prob result = { (uint8_t) (sigma0(0 ^ d_m) ^ sigma0(0)), 0.0f };
    return result;
}

// Linear
struct diff_prob propagate_sigma1(uint8_t d_m)
{
    struct diff_prob result = { (uint8_t) (sigma1(0 ^ d_m) ^ sigma1(0)), 0.0f };
    return result;
}

// Nonlinear
struct diff_vec propagate_keymix(uint8_t d_x, size_t round, float l2pthresh)
{
    const uint8_t K[16] = 
    {
        0xb7, 0xe1, 0x51, 0x62,
        0x8a, 0xed, 0x2a, 0x6a,
        0xbf, 0x71, 0x58, 0x80,
        0x9c, 0xf4, 0xf3, 0xc7
    };

    const size_t sample_size = 256;
    size_t counts[256];
    memset(counts, 0, 256 * sizeof(size_t));

    // Sample
    for (int x = 0; x < 256; x++)
    {
        counts[add(K[round], x ^ d_x) ^ add(K[round], x)]++;
    }

    // Filter insignificant results
    struct diff_vec samples = { NULL, 0 };
    for (int i = 0; i < 256; i++) if (counts[i])
    {
        if (log2f(counts[i]) - log2f(sample_size) > l2pthresh) samples.len++;
        else counts[i] = 0;
    }

    // Allocate enough space for these results
    samples.pairs = (struct diff_prob *) calloc(samples.len, sizeof(struct diff_prob));

    // Copy the results in
    for (int i = 0, k = 0; i < 256; i++) if (counts[i])
    {
        samples.pairs[k++] = { (uint8_t) i, log2f(counts[i]) - log2f(sample_size) };
    }
    return samples;
}

// Nonlinear
struct diff_vec propagate_add(uint8_t d_x, uint8_t d_y, float l2pthresh)
{
    const size_t sample_size = 256*256;
    size_t counts[256];
    memset(counts, 0, 256 * sizeof(size_t));

    // Sample
    for (int x = 0; x < 256; x++) for (int y = 0; y < 256; y++)
    {
        counts[add(x ^ d_x, y ^ d_y) ^ add(x, y)]++;
    }
    
    // Filter insignificant results
    struct diff_vec samples = { NULL, 0 };
    for (int i = 0; i < 256; i++) if (counts[i])
    {
        if (log2f(counts[i]) - log2f(sample_size) > l2pthresh) samples.len++;
        else counts[i] = 0; 
    }

    // Allocate enough space for these results
    samples.pairs = (struct diff_prob *) calloc(samples.len, sizeof(struct diff_prob));

    // Copy the results in
    for (int i = 0, k = 0; i < 256; i++) if (counts[i])
    {
        samples.pairs[k++] = { (uint8_t) i, log2f(counts[i]) - log2f(sample_size) };
    }
    return samples;
}

// Nonlinear
struct diff_vec propagate_maj(uint8_t d_x, uint8_t d_y, uint8_t d_z, float l2pthresh)
{
    const size_t sample_size = 256*256*256;
    size_t counts[256];
    memset(counts, 0, 256 * sizeof(size_t));
    
    // Sample
    for (int x = 0; x < 256; x++) for (int y = 0; y < 256; y++) for (int z = 0; z < 256; z++)
    {
        counts[maj(x ^ d_x, y ^ d_y, z ^ d_z) ^ maj(x, y, z)]++; 
    }

    // Filter insignificant results
    struct diff_vec samples = { NULL, 0 };
    for (int i = 0; i < 256; i++) if (counts[i])
    {
        if (log2f(counts[i]) - log2f(sample_size) > l2pthresh) samples.len++;
        else counts[i] = 0;
    }

    // Allocate enough space for these results
    samples.pairs = (struct diff_prob *) calloc(samples.len, sizeof(struct diff_prob));

    // Copy the results in
    for (int i = 0, k = 0; i < 256; i++) if (counts[i])
    {
        samples.pairs[k++] = { (uint8_t) i, log2f(counts[i]) - log2f(sample_size) };
    }
    return samples;
}

// Take a differential, and propagate through to round n.
void propagate(const char *msg_diff, const size_t n, const float pthresh)
{
    // 16 rounds max
    if (n > 16)
    {
        fprintf(stderr, "Cannot propagate over more than 16 rounds\n");
        return;
    }
    // For backtracking
    std::stack<std::pair<struct prop_state, struct diff_vec>> stack;

    // Build the default state objects
    struct prop_state sstate;
    memset(&sstate, 0, sizeof(sstate));
    struct diff_vec svec = { NULL, 0 };

    // Set up the message schedule differentials; first 8 are concrete, but the
    // last 8 will be variable (non-linear transforms of the message)
    for (int i = 0; i < 8; i++) 
    {
        uint8_t msg = 0x00;
        for (int j = 0; j < 8; j++)
        {
            msg |= msg_diff[8 * i + j] == 'x'? 1 : 0;
            msg <<= 1;
        }
        sstate.sched[i] = msg;
    }

    // Push our start state onto the stack to to kickstart the propagation process
    stack.push(std::make_pair(sstate, svec));
   
    bool b = false;
    // Stack contains points where we can restart a propagation
    while (!stack.empty())
    {
        // Grab the propagation state on top
        struct prop_state state = stack.top().first;
        // !!! Somehow need for force sstate off the stack
        if (b && !memcmp(&state, &sstate, sizeof(struct prop_state))) break;
        b = true;
        // !!!

        // Run this propagation through to completion
        while (state.round < n) switch (state.step)
        {
            #define PROP_START(...)                                             \
            if (memcmp(&state, &stack.top().first, sizeof(struct prop_state)))  \
            {                                                                   \
                struct diff_vec diffs = __VA_ARGS__;                            \
                if (diffs.len == 0) goto BAILOUT;                               \
                stack.push(std::make_pair(state, diffs));                       \
            }
            #define PROP_INTROS                                   \
            struct diff_vec *diffs = &stack.top().second;         \
            struct diff_prob pair  = diffs->pairs[--diffs->len];
            #define PROP_END        \
            if (diffs->len == 0)    \
            {                       \
                free(diffs->pairs); \
                stack.pop();        \
            }
            #define t (state.round)

            case 0:
                {
                    // t1 = sigma1(b)
                    state.t1    = propagate_sigma1(state.b).diff;
                    state.step += 1;
                    break;
                }

            case 1:
                {
                    // t1 = t1 + d
                    PROP_START(propagate_add(state.t1, state.d, pthresh));
                    PROP_INTROS
                    state.l2prob += pair.l2prob;
                    state.t1      = pair.diff;
                    state.step   += 1;
                    PROP_END
                    break;
                }

            case 2:
                {
                    // t1 = t1 + K[t]
                    PROP_START(propagate_keymix(state.t1, state.round, pthresh));
                    PROP_INTROS
                    state.l2prob += pair.l2prob;
                    state.t1      = pair.diff;
                    state.step   += t < 8? 3 : 1;
                    PROP_END
                    break;
                }

            case 3:
                {
                    // t >= 8: W[t] = sigma0(W[t-3]) + W[t-4]
                    PROP_START(propagate_add(propagate_sigma0(state.sched[t-3]).diff, state.sched[t-4], pthresh));
                    PROP_INTROS
                    state.l2prob  += pair.l2prob;
                    state.sched[t] = pair.diff;
                    state.step    += 1;
                    PROP_END
                    break;
                }

            case 4:
                {
                    // t >= 8: W[t] = sigma1(W[t-8]) + W[t]
                    PROP_START(propagate_add(propagate_sigma1(state.sched[t-8]).diff, state.sched[t], pthresh));
                    PROP_INTROS
                    state.l2prob  += pair.l2prob;
                    state.sched[t] = pair.diff;
                    state.step    += 1;
                    PROP_END
                    break;
                }

            case 5:
                {
                    // t1 = t1 + W[t]
                    PROP_START(propagate_add(state.t1, state.sched[t], pthresh));
                    PROP_INTROS
                    state.l2prob += pair.l2prob;
                    state.t1      = pair.diff;
                    state.step   += 1;
                    PROP_END
                    break;
                }
            case 6:
                {
                    // t2 = sigma0(a)
                    state.t2    = propagate_sigma0(state.a).diff;
                    state.step += 1;
                    break;
                }
            case 7:
                {
                    // maj = maj(a, b, c)
                    PROP_START(propagate_maj(state.a, state.b, state.c, pthresh));
                    PROP_INTROS
                    state.l2prob += pair.l2prob;
                    state.maj     = pair.diff;
                    state.step   += 1;
                    PROP_END
                    break;
                }

            case 8:
                {
                    // t2 = t2 + maj
                    PROP_START(propagate_add(state.t2, state.maj, pthresh));
                    PROP_INTROS
                    state.l2prob += pair.l2prob;
                    state.t2      = pair.diff;
                    state.step   += 1;
                    PROP_END
                    break;
                }

            case 9:
                {
                    // d = c; c = b + t1
                    PROP_START(propagate_add(state.b, state.t1, pthresh));
                    PROP_INTROS
                    state.l2prob += pair.l2prob;
                    state.d       = state.c;
                    state.c       = pair.diff;
                    state.step   += 1;
                    PROP_END
                    break;
                }

            case 10:
                {
                    // b = a; a = t1 + t2
                    PROP_START(propagate_add(state.t1, state.t2, pthresh));
                    PROP_INTROS
                    state.l2prob += pair.l2prob;
                    state.b       = state.a;
                    state.a       = pair.diff;
                   
                    // Copy registers into trail
                    state.trail[4*t+0] = state.a;
                    state.trail[4*t+1] = state.b;
                    state.trail[4*t+2] = state.c;
                    state.trail[4*t+3] = state.d;

                    // Reset temps
                    state.t1     = 0;
                    state.t2     = 0;
                    state.maj    = 0;
                    state.step   = 0;
                    state.round += 1;
                    PROP_END
                    break;
                }
            
            default: // Impossible
                break;

        #undef PROP_START
        #undef PROP_INTROS
        #undef PROP_END
        #undef t
        }
        
        // Finished our search
        if (state.round == n && state.diff == 0)
        {
            printf("0 -> ");
            for (int i = 0; i < 16; i++)
            {
                printf("%x -> ", state.trail32[i]);
            }
            printf("* ~ 2^%f\n", state.l2prob);
        }
BAILOUT: continue;
    }
}

// Create an input differential randomly
char *make_input_diff()
{
    // Allocate a block of memory
    char *diff = (char *) malloc(BLOCK_SIZE + 1);
    // Guess that the sparse area should be the first 4 words
    memset(diff, '-', 4*8);
    // Randomly assign the remaining dense area
    for (size_t idx = 4*8; idx < BLOCK_SIZE; idx++)
    {
        diff[idx] = (rand() & 1)? '-' : 'x';
    }
    diff[BLOCK_SIZE] = '\0';
    return diff;
}

// Entry point
int main(int argc, char **argv)
{
    unsigned int seed = time(NULL);
    srand(seed);
    float pthresh = -4.0f;

    do
    {
        seed = rand();
        
        // Set up the PRNG
        fprintf(stdout, "SEED: %u\n", seed);
        srand(seed);

        // Build an input differential
        char *diff = make_input_diff();
        size_t rounds = 8;

        // Print the differential + round count
        for (int i = 1; i <= BLOCK_SIZE; i++)
        {
            putc(diff[i-1], stdout);
            if (!(i % 8)) putc('\n', stdout);
        }
        fprintf(stdout, "Rounds: %zu\n", rounds);
        fprintf(stdout, "Threshold probability: 2^%f\n\n", pthresh);
        fflush(stdout);
        propagate(diff, rounds, pthresh);
        free(diff);
    }
    while (1);
    return 0;
}
