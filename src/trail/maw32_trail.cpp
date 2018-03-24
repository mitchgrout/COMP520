// libc
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// STL
#include <stack>
#include <vector>
#include <map>
using namespace std;

// Consts
#define BLOCK_SIZE 64
#define DIGEST_SIZE 32

// MAW internals
static inline uint8_t rotr(const uint8_t x, const uint8_t n) { return (x >> n) | (x << (8 - n)); }
static inline uint8_t maj(const uint8_t x, const uint8_t y, const uint8_t z) { return (x & y) ^ (x & z) ^ (y & z); }
static inline uint8_t sigma0(const uint8_t x) { return rotr(x, 2) ^ rotr(x, 3) ^ rotr(x, 5); }
static inline uint8_t sigma1(const uint8_t x) { return rotr(x, 1) ^ rotr(x, 4) ^ (x >> 3); }
static inline uint8_t add(const uint8_t x, const uint8_t y) { return x + y; }

// Propagation state
struct prop_state
{
    size_t round;               // Which round we are currently propagating
    size_t step;                // Which step in particular
    uint8_t sched[16];          // Differences in message schedule
    union 
    {
        struct { uint8_t a, b, c, d; }; // Current registers
        uint32_t diff;          // Registers as a 32-bit int
    };
    uint8_t t1, t2, maj;        // Differences in temporary variables
    union
    {
        uint8_t trail[4*16];    // Trail of differences
        uint32_t trail32[16];   //
    };
};

// Util 
static inline int prop_state_equal(const struct prop_state left, const struct prop_state right)
{
    return left.round == right.round &&
           left.step  == right.step; 
}

// Linear
static inline uint8_t propagate_sigma0(const uint8_t d_m)
{
    return sigma0(0 ^ d_m) ^ sigma0(0);
}

// Linear
static inline uint8_t propagate_sigma1(const uint8_t d_m)
{
    return sigma1(0 ^ d_m) ^ sigma1(0);
}

// Util
static inline vector<uint8_t> filter(const map<uint8_t, size_t>& samples, const size_t sample_size, const float l2pthresh)
{
    vector<uint8_t> results;
    for (const auto& elem : samples)
    {
        float prob = log2f(elem.second) - log2f(sample_size);
        if (prob >= l2pthresh) results.push_back(elem.first);
    }
    return results;
}

// Nonlinear
static map<uint16_t, vector<uint8_t>> key_memo;
vector<uint8_t> propagate_keymix(const uint8_t d_x, const size_t round, const float l2pthresh)
{
    // Memoization
    uint16_t key = ((d_x << 8) | round) & 0xffff;
    if (key_memo.find(key) != key_memo.end()) return key_memo[key];

    const uint8_t K[16] = 
    {
        0xb7, 0xe1, 0x51, 0x62, 0x8a, 0xed, 0x2a, 0x6a,
        0xbf, 0x71, 0x58, 0x80, 0x9c, 0xf4, 0xf3, 0xc7
    };

    const size_t sample_size = 256;
    map<uint8_t, size_t> counts;

    // Sample
    for (int n = 0; n < sample_size; n++)
    {
        uint8_t x = n;
        counts[add(K[round], x ^ d_x) ^ add(K[round], x)]++;
    }
    vector<uint8_t> result = filter(counts, sample_size, l2pthresh);
    key_memo[key] = result;
    return result;
}

// Nonlinear
static map<uint16_t, vector<uint8_t>> add_memo;
vector<uint8_t> propagate_add(const uint8_t d_x, const uint8_t d_y, const float l2pthresh)
{
    // Memoization
    uint16_t key = ((d_x << 8) | d_y) & 0xffff;
    if (add_memo.find(key) != add_memo.end()) return add_memo[key];

    const size_t sample_size = 256*256;
    map<uint8_t, size_t> counts;

    // Sample
    for (int n = 0; n < sample_size; n++)
    {
        uint8_t x = (n >> 8) & 0xff,
                y = (n >> 0) & 0xff;
        counts[add(x ^ d_x, y ^ d_y) ^ add(x, y)]++;
    }
    vector<uint8_t> result = filter(counts, sample_size, l2pthresh);
    add_memo[key] = result;
    return result;
}

// Nonlinear
static map<uint32_t, vector<uint8_t>> maj_memo;
vector<uint8_t> propagate_maj(const uint8_t d_x, const uint8_t d_y, const uint8_t d_z, const float l2pthresh)
{
    // Memoization
    uint32_t key = ((d_x << 16) | (d_y << 8) | d_z);
    if (maj_memo.find(key) != maj_memo.end()) return maj_memo[key];

    const size_t sample_size = 256*256;
    map<uint8_t, size_t> counts;

    // Sample
    // for (int x = 0; x < 256; x++) for (int y = 0; y < 256; y++) for (int z = 0; z < 256; z++)
    for (int n = 0; n < sample_size; n++)
    {
        uint8_t x = rand() & 0xff,
                y = rand() & 0xff,
                z = rand() & 0xff;
        counts[maj(x ^ d_x, y ^ d_y, z ^ d_z) ^ maj(x, y, z)]++; 
    }
    vector<uint8_t> result = filter(counts, sample_size, l2pthresh);
    maj_memo[key] = result;
    return result;
}

// Take a differential, and propagate through to round n.
void propagate(const char *msg_diff, size_t n, const float pthresh)
{
    // 16 rounds max
    if (n > 16)
    {
        n = 16;
        fprintf(stderr, "Warning: Cannot propagate over more than 16 rounds maximum.\n");
        return;
    }

    // Statistics
    size_t total_trails = 0,
           zero_trails  = 0;
    time_t start_time = time(NULL);

    // Backtracking value
    stack<pair<struct prop_state, vector<uint8_t>>> stack;

    // Build the default state objects
    struct prop_state sstate;
    memset(&sstate, 0, sizeof(sstate));
    vector<uint8_t> svec;

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
    stack.push(make_pair(sstate, svec));
   
    bool b = false;
    // Stack contains points where we can restart a propagation
    while (!stack.empty())
    {
        // Grab the propagation state on top
        struct prop_state state = stack.top().first;
        // !!! Somehow need for force sstate off the stack
        if (b && prop_state_equal(sstate, state)) break;
        b = true;
        // !!!

        // Run this propagation through to completion
        while (state.round < n)
        {
            // We can actually use a quick heuristic here: 
            // registers a, c, are not changed, only moved to new registers. So
            // to get an output difference of 0 in the final round, a=c=0.
            if (state.round == n - 1 && (state.a != 0 || state.c != 0)) goto BAILOUT;
            switch (state.step)
            {
                #define PROP_START(...)                         \
                if (!prop_state_equal(state, stack.top().first))\
                {                                               \
                    vector<uint8_t> vec = __VA_ARGS__;          \
                    if (vec.size() == 0) goto BAILOUT;          \
                    stack.push(make_pair(state, vec));          \
                }

                #define PROP_INTROS                         \
                vector<uint8_t>& vec = stack.top().second;  \
                uint8_t diff         = vec.back();          \
                vec.pop_back();

                #define PROP_END if (vec.size() == 0) stack.pop();

                #define t (state.round)

                case 0:
                    {
                        // t1 = sigma1(b)
                        state.t1    = propagate_sigma1(state.b);
                        state.step += 1;
                        break;
                    }

                case 1:
                    {
                        // t1 = t1 + d
                        PROP_START(propagate_add(state.t1, state.d, pthresh));
                        PROP_INTROS
                            state.t1      = diff;
                        state.step   += 1;
                        PROP_END
                            break;
                    }

                case 2:
                    {
                        // t1 = t1 + K[t]
                        PROP_START(propagate_keymix(state.t1, state.round, pthresh));
                        PROP_INTROS
                            state.t1      = diff;
                        state.step   += t < 8? 3 : 1;
                        PROP_END
                            break;
                    }

                case 3:
                    {
                        // t >= 8: W[t] = sigma0(W[t-3]) + W[t-4]
                        PROP_START(propagate_add(propagate_sigma0(state.sched[t-3]), state.sched[t-4], pthresh));
                        PROP_INTROS
                            state.sched[t] = diff;
                        state.step    += 1;
                        PROP_END
                            break;
                    }

                case 4:
                    {
                        // t >= 8: W[t] = sigma1(W[t-8]) + W[t]
                        PROP_START(propagate_add(propagate_sigma1(state.sched[t-8]), state.sched[t], pthresh));
                        PROP_INTROS
                            state.sched[t] = diff;
                        state.step    += 1;
                        PROP_END
                            break;
                    }

                case 5:
                    {
                        // t1 = t1 + W[t]
                        PROP_START(propagate_add(state.t1, state.sched[t], pthresh));
                        PROP_INTROS
                            state.t1      = diff;
                        state.step   += 1;
                        PROP_END
                            break;
                    }
                case 6:
                    {
                        // t2 = sigma0(a)
                        state.t2    = propagate_sigma0(state.a);
                        state.step += 1;
                        break;
                    }
                case 7:
                    {
                        // maj = maj(a, b, c)
                        PROP_START(propagate_maj(state.a, state.b, state.c, pthresh));
                        PROP_INTROS
                            state.maj     = diff;
                        state.step   += 1;
                        PROP_END
                            break;
                    }

                case 8:
                    {
                        // t2 = t2 + maj
                        PROP_START(propagate_add(state.t2, state.maj, pthresh));
                        PROP_INTROS
                            state.t2      = diff;
                        state.step   += 1;
                        PROP_END
                            break;
                    }

                case 9:
                    {
                        // d = c; c = b + t1
                        PROP_START(propagate_add(state.b, state.t1, pthresh));
                        PROP_INTROS
                            state.d       = state.c;
                        state.c       = diff;
                        state.step   += 1;
                        PROP_END
                            break;
                    }

                case 10:
                    {
                        // b = a; a = t1 + t2
                        PROP_START(propagate_add(state.t1, state.t2, pthresh));
                        PROP_INTROS
                            state.b          = state.a;
                        state.a          = diff;
                        state.trail32[t] = state.diff;
                        state.step       = 0;
                        state.round     += 1;
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
        }
        // Finished our search
        if (state.round == n)
        {
            total_trails++;
            zero_trails += !state.diff;
        }
        if (0)
        {
BAILOUT:
            total_trails++;
        }
    }
    time_t end_time = time(NULL);
    printf("Elapsed time: %ld sec\n", end_time - start_time);
    printf("Zero trails: %zu/%zu\n\n", zero_trails, total_trails);
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
    const float pthresh = -4.0f;

    do
    {
        seed = rand();
        
        // Set up the PRNG
        fprintf(stdout, "SEED: %u\n", seed);
        srand(seed);

        // Build an input differential
        char *diff = make_input_diff();
        size_t rounds = 16;

        // Print the differential + round count
        for (int i = 1; i <= BLOCK_SIZE; i++)
        {
            putc(diff[i-1], stdout);
            if (!(i % 8)) putc('\n', stdout);
        }
        fprintf(stdout, "Rounds: %zu\n", rounds);
        fprintf(stdout, "Threshold probability: 2^%f\n", pthresh);
        fprintf(stdout, "==================================\n");
        propagate(diff, rounds, pthresh);
        fflush(stdout);
        free(diff);
    }
    while (1);
    return 0;
}
