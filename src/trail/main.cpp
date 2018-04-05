// Entry point to maw_trail

#ifndef __MAIN
#define __MAIN

#include <stdlib.h>
#include <stdint.h>
#include "utils.cpp"
#include "maw32_utils.cpp"
#include "maw32_trail.cpp"
#include <utility>

// Determine if an input differential is viable
static inline bool is_viable(uint8_t *W, const size_t rounds, const float l2pthresh,
                             size_t t, size_t ctr)
{
    if (t >= 16) return ctr >= 0;   // 2/8
    uint8_t w0 = sigma0(W[t-3]),
            w1 = sigma1(W[t-8]);
    for (uint8_t t1: propagate_add(w0, w1, l2pthresh))
    for (uint8_t t2: propagate_add(W[t-4], t1, l2pthresh))
    {
        W[t] = t2;
        if (is_viable(W, rounds, l2pthresh, t+1, ctr + !t2)) return true;
    }
    return false;
}

// Successor of an array of bytes
static inline bool succ(uint8_t *ptr, size_t len)
{
    // Increment elements right-to-left, stopping if there is no carry
    while (len) if (++ptr[--len]) return true;
    return false;
}

// Create an input differential randomly
uint8_t *make_input_diff(size_t rounds, float l2pthresh)
{
    // The entire message schedule
    static uint8_t sched[16];

    // Randomly assign differences for last three words
    for (int idx = 0; idx < 8; idx++)
    {
        sched[idx] = (idx < 4)? 0 : rand();
    }
    // If not viable, mutate slightly
    while (!is_viable(sched, rounds, l2pthresh, 8, 0)) succ(sched, 8) ;
    return sched;
}

// Entry point
int main(int argc, char **argv)
{
    const float pthresh = -4.000000f;
    const size_t rounds = 16;
    const char *key_fname = "/Scratch/key-file-4.000000.bin",
               *add_fname = "/Scratch/add-file-4.000000.bin",
               *maj_fname = "/Scratch/maj-file-4.000000.bin";

    log(stdout, "Initializing...");
    if (load_key_memo(key_fname)) log(stdout, "Loaded key memos from %s", key_fname);
    else log(stdout, "Failed to load key memos from %s", key_fname);
    if (load_add_memo(add_fname)) log(stdout, "Loaded add memos from %s", add_fname);
    else log(stdout, "Failed to load add memos from %s", add_fname);
    if (load_maj_memo(maj_fname)) log(stdout, "Loaded maj memos from %s", maj_fname);
    else log(stdout, "Failed to load maj memos from %s", maj_fname);
    log(stdout, "Done!");

    // Init PRNG
    unsigned int seed = time(NULL);
    srand(seed);
    do
    {
        // Reseed so we can deterministically perform a propagation
        seed = rand();
        srand(seed);

        // Build an input differential, print it out
        uint8_t *diff = make_input_diff(rounds, pthresh);
        log(stdout, "Fingerprint:");
        for (int i = 0; i < 8; i++)
        {
            #define C2B(c) ((c)? 'x' : '-')
            log(stdout, "%c%c%c%c%c%c%c%c", 
                    C2B(diff[i] & 0x80),
                    C2B(diff[i] & 0x40),
                    C2B(diff[i] & 0x20),
                    C2B(diff[i] & 0x10),
                    C2B(diff[i] & 0x08),
                    C2B(diff[i] & 0x04),
                    C2B(diff[i] & 0x02),
                    C2B(diff[i] & 0x01));
            #undef C2B
        }
        log(stdout, "Seed: %u", seed);
        log(stdout, "Rounds: %zu/16", rounds);
        log(stdout, "Threshold probability: 2^%f", pthresh);
        std::pair<size_t, size_t> result = propagate(diff, rounds, pthresh);
        log(stdout, "Trails: %zu/%zu\n", result.first, result.second);
    }
    while (1);
    return 0;
}
#endif // __MAIN
