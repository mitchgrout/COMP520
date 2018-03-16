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
    uint8_t d0, d1, d2, d3, d4, d5; // Differences in temporary variables
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
void propagate(const char *msg_diff, const size_t n)
{
    // Threshhold probability
    const float pthresh = -5.0f;

    // Stack is used to track propagations which had multiple output differences
    // that we need to check. If a prop_state is on the stack, then it should
    // have a set of diffs to check. These diffs will be read off in reverse
    // order
    // Stack is allowed to increase in size only; should grow by a factor of 2
    /*
    static ssize_t stack_size       = 1024;
    static ssize_t stack_pos        = -1;
    static struct prop_state *stack = (struct prop_state *) calloc(stack_size, sizeof(struct prop_state));
    #define STACK_INCREASE() do { stack_size *= 2; stack = (struct prop_state *) realloc(stack, stack_size * sizeof(struct prop_state)); } while (0)
    #define STACK_IS_EMPTY() (stack_pos < 0)
    #define STACK_POP()      (stack_pos--)
    #define STACK_PEEK()     (stack[stack_pos])
    #define STACK_PUSH(x)    do { if (stack_pos >= stack_size) { STACK_INCREASE(); } stack[++stack_pos] = x; } while (0)
    */
    // Info used for backtracking: a base state, and a set of diff/prob values
    // which could have been used for propagation
    std::stack<std::pair<struct prop_state, struct diff_vec>> stack;

    // Build the default prop_state object
    struct prop_state sstate;
    sstate.current_round = 0;
    sstate.current_step  = 0;
    sstate.l2prob        = 0.0f;
    sstate.differential  = 0;
    sstate.d0 = 0; sstate.d1 = 0; sstate.d2 = 0;
    sstate.d3 = 0; sstate.d4 = 0; sstate.d5 = 0;
    struct diff_vec svec = { NULL, 0 };
    stack.push(std::make_pair(sstate, svec));
    
    // Stack contains points where we can restart a propagation
    while (!stack.empty())
    {
START:
        // Grab the propagation state on top
        struct prop_state state = stack.top().first;

        // Run this propagation through to completion
        while (state.current_round < n)
        {
#if 0
            // Give an indicator of stack size [display what round/step we're on]
            for (int i = 0; i < 8 * state.current_round + state.current_step; i++)
                putchar('#');
            putchar('\n');
#endif
            switch (state.current_step)
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
                        // d1 = propagate_add(d0, m)
                        // Fresh run
                        if (memcmp(&state, &stack.top().first, sizeof(struct prop_state)))
                        {
                            // Mix in message difference
                            uint8_t m = 0x00;
                            for (int i = 0; i < 8; i++) 
                                m |= (msg_diff[8 * state.current_round + i] == 'x'? 1 : 0) << (8 - i);

                            struct diff_vec diffs = propagate_add(state.d0, m, pthresh);
                            if (diffs.len == 0) goto BAILOUT; // Dead end
                            stack.push(std::make_pair(state, diffs));
                            goto START; // Restart
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
                            goto START; // Restart
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
                            goto START; // Restart
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
                            goto START; // Restart
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
                            goto START; // Restart
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
                            goto START; // Restart
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
        }
BAILOUT:
        // Either hit a bad probability and bailed out, or finished. In either
        // case we will backtrack to last item on stack
        if (state.current_round == n && state.differential == 0)
        {
            // We hit the end and found a differential
            printf("Differential (0x%02x,0x%02x,0x%02x,0x%02x) ~ 2^%f\n", state.a, state.b, state.c, state.d, state.l2prob);
        }
        stack.pop();
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
    // Set up the PRNG
    unsigned int seed = time(NULL);
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
    fprintf(stdout, "Rounds: %zu\n\n", rounds);
    fflush(stdout);
    propagate(diff, rounds);
    // NOTE: DONT TRY TO FREE(DIFF) FOR SOME REASON
    return 0;
}
