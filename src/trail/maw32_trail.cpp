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

#define LOG2(x) log2f(x)

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
    size_t current_round;           // Which round we are currently propagating
    size_t current_step;            // Which step in particular
    float l2prob;                   // Probability that our differential holds in log2
    union 
    {
        uint8_t a, b, c, d;         // Differences in every register
        uint32_t differential;      // Difference as a 32-bit int
    };
    union
    {
        uint8_t d0, d1, d2, d3, d4, d5; // Differences in temporary variables
        uint64_t tmp_diffs;             // Differences as a 64-bit int
    };
    uint8_t W[16];                  // Differences in message schedule
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
    // For backtracking
    std::stack<std::pair<struct prop_state, struct diff_vec>> stack;

    // Build the default prop_state object
    struct prop_state sstate;
    sstate.current_round = 0;
    sstate.current_step  = 0;
    sstate.l2prob        = 0.0f;
    sstate.differential  = 0;
    sstate.tmp_diffs     = 0;
    // Set up the message schedule differentials; first 8 are concrete, but the
    // last 8 can be variable
    memset(sstate.W, 0, sizeof(sstate.W));
    for (int i = 0; i < 8; i++) 
    {
        uint8_t m = 0x00;
        for (int j = 0; j < 8; j++)
        {
            m |= msg_diff[8 * i + j] == 'x'? 1 : 0;
            m <<= 1;
        }
        sstate.W[i] = m;
    }
    
    struct diff_vec svec = { NULL, 0 };
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
        while (state.current_round < n) switch (state.current_step)
        {
            case 0:
                {
                    // d0 = propagate_sigma1(b)
                    struct diff_prob pair = propagate_sigma1(state.b);
                    state.d0              = pair.diff;
                    state.current_step   += 1;
                    break;
                }
            case 1:
                {
                    // d1 = propagate_add(d0, W[t])
                    // Fresh run
                    if (memcmp(&state, &stack.top().first, sizeof(struct prop_state)))
                    {
                        struct diff_vec diffs = propagate_add(state.d0, state.W[state.current_round], pthresh);
                        if (diffs.len == 0) goto BAILOUT; // Dead end
                        stack.push(std::make_pair(state, diffs));
                    }

                    // Pop the back element from the vec on the stack
                    struct diff_vec *diffs = &stack.top().second;
                    struct diff_prob pair  = diffs->pairs[--diffs->len];
                    state.l2prob          += pair.l2prob;
                    state.d1               = pair.diff;
                    state.current_step    += 1;

                    // If that was the last element, pop from stack
                    if (diffs->len == 0)
                    {
                        free(diffs->pairs);
                        stack.pop();
                    }
                    break;
                }
            case 2:
                {
                    // d2 = propagate_add(d, d1)
                    // Fresh run
                    if (memcmp(&state, &stack.top().first, sizeof(struct prop_state)))
                    {
                        struct diff_vec diffs = propagate_add(state.d, state.d1, pthresh);
                        if (diffs.len == 0) goto BAILOUT; // Dead end
                        stack.push(std::make_pair(state, diffs));
                    }

                    // Pop the back element from the vec on the stack
                    struct diff_vec *diffs = &stack.top().second;
                    struct diff_prob pair  = diffs->pairs[--diffs->len];
                    state.l2prob          += pair.l2prob;
                    state.d2               = pair.diff;
                    state.current_step    += 1;

                    // If that was the last element, pop from stack
                    if (diffs->len == 0)
                    {
                        free(diffs->pairs);
                        stack.pop();
                    }
                    break;
                }
            case 3:
                {
                    // d3 = propagate_sigma0(a)
                    struct diff_prob pair = propagate_sigma0(state.a);
                    state.d3              = pair.diff;
                    state.current_step   += 1;
                    break;
                }
            case 4:
                {
                    // d4 = propagate_maj(a, b, c)
                    // Fresh run
                    if (memcmp(&state, &stack.top().first, sizeof(struct prop_state)))
                    {
                        struct diff_vec diffs = propagate_maj(state.a, state.b, state.c, pthresh);
                        if (diffs.len == 0) goto BAILOUT; // Dead end
                        stack.push(std::make_pair(state, diffs));
                    }

                    // Pop the back element from the vec on the stack
                    struct diff_vec *diffs = &stack.top().second;
                    struct diff_prob pair  = diffs->pairs[--diffs->len];
                    state.l2prob          += pair.l2prob;
                    state.d4               = pair.diff;
                    state.current_step    += 1;

                    // If that was the last element, pop from stack
                    if (diffs->len == 0)
                    {
                        free(diffs->pairs);
                        stack.pop();
                    }
                    break;
                }
            case 5:
                {
                    // d5 = propagate_add(d3, d4)
                    // Fresh run
                    if (memcmp(&state, &stack.top().first, sizeof(struct prop_state)))
                    {
                        struct diff_vec diffs = propagate_add(state.d3, state.d4, pthresh);
                        if (diffs.len == 0) goto BAILOUT; // Dead end
                        stack.push(std::make_pair(state, diffs));
                    }

                    // Pop the back element from the vec on the stack
                    struct diff_vec *diffs = &stack.top().second;
                    struct diff_prob pair  = diffs->pairs[--diffs->len];
                    state.l2prob          += pair.l2prob;
                    state.d5               = pair.diff;
                    state.current_step    += 1;

                    // If that was the last element, pop from stack
                    if (diffs->len == 0)
                    {
                        free(diffs->pairs);
                        stack.pop();
                    }
                    break;
                }
            case 6:
                {
                    // d := c; c := propagate_add(b, d2)
                    // Fresh run
                    if (memcmp(&state, &stack.top().first, sizeof(struct prop_state)))
                    {
                        struct diff_vec diffs = propagate_add(state.d, state.d1, pthresh);
                        if (diffs.len == 0) goto BAILOUT; // Dead end
                        stack.push(std::make_pair(state, diffs));
                    }

                    // Pop the back element from the vec on the stack
                    struct diff_vec *diffs = &stack.top().second;
                    struct diff_prob pair  = diffs->pairs[--diffs->len];
                    state.l2prob          += pair.l2prob;
                    state.d                = state.c;
                    state.c                = pair.diff;
                    state.current_step    += 1;

                    // If that was the last element, pop from stack
                    if (diffs->len == 0)
                    {
                        free(diffs->pairs);
                        stack.pop();
                    }
                    break;
                }
            case 7:
                {
                    // b := a; a := propagate_add(d2, d5)
                    // Fresh run
                    if (memcmp(&state, &stack.top().first, sizeof(struct prop_state)))
                    {
                        struct diff_vec diffs = propagate_add(state.d2, state.d5, pthresh);
                        if (diffs.len == 0) goto BAILOUT; // Dead end
                        stack.push(std::make_pair(state, diffs));
                    }

                    // Pop the back element frmo the vec on the stack
                    struct diff_vec *diffs = &stack.top().second;
                    struct diff_prob pair  = diffs->pairs[--diffs->len];
                    state.l2prob          += pair.l2prob;
                    state.b                = state.a;
                    state.a                = pair.diff;
                    state.current_step     = 0;
                    state.current_round   += 1;
                    state.d0 = state.d1 = state.d2 = state.d3 = state.d4 = state.d5 = 0;

                    // If that was the last element, pop from stack
                    if (diffs->len == 0)
                    {
                        free(diffs->pairs);
                        stack.pop();
                    }
                    break;
                }
            default: // Impossible
                break;
        }
        
        // Finished our search
        if (state.current_round == n && state.differential == 0)
        {
            printf("Differential (0x%02x,0x%02x,0x%02x,0x%02x) ~ 2^%f\n", state.a, state.b, state.c, state.d, state.l2prob);
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
    // Consts
    unsigned int seed = time(NULL);
    float pthresh = -4.0f;

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
    fprintf(stdout, "Threshold probability: 2^-%f\n\n", pthresh);
    fflush(stdout);
    propagate(diff, rounds, pthresh);
    free(diff);
    return 0;
}
