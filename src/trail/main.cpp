// Entry point to maw_trail

#ifndef __MAIN
#define __MAIN

#include <stdlib.h>
#include <stdint.h>
#include "utils.cpp"
#include "maw32_utils.cpp"
#include "maw32_trail.cpp"
#include <utility>

// Read a byte from /dev/urandom
static inline uint8_t csprng_rand()
{
    static FILE *handle = fopen("/dev/urandom", "r");
    return getc(handle);
}

// Determine if an input differential is viable
static inline bool is_viable(uint8_t *W, const size_t rounds, const float l2pthresh,
                             size_t t, size_t ctr)
{
    if (t >= rounds) return true;
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

// Create an input differential randomly
uint8_t *make_input_diff(size_t rounds, float l2pthresh)
{
    // The entire message schedule
    static uint8_t sched[16];

    // Randomly assign differences for last three words
    do for (int idx = 0; idx < 8; idx++) sched[idx] = (idx < 4)? 0 : csprng_rand();
    while (!is_viable(sched, rounds, l2pthresh, 8, 0)) ;
    return sched;
}

// Entry point
int main(int argc, char **argv)
{
    const float pthresh = -3.000000f;
    const size_t rounds = 12;
    const char *key_fname = "/Scratch/key-file-3.000000.bin",
               *add_fname = "/Scratch/add-file-3.000000.bin",
               *maj_fname = "/Scratch/maj-file-3.000000.bin";

    log(stdout, "Initializing...");
    if (load_key_memo(key_fname)) log(stdout, "Loaded key memos from %s", key_fname);
    else log(stdout, "Failed to load key memos from %s", key_fname);
    if (load_add_memo(add_fname)) log(stdout, "Loaded add memos from %s", add_fname);
    else log(stdout, "Failed to load add memos from %s", add_fname);
    if (load_maj_memo(maj_fname)) log(stdout, "Loaded maj memos from %s", maj_fname);
    else log(stdout, "Failed to load maj memos from %s", maj_fname);
    log(stdout, "Done!\n");
    log(stdout, "Rounds: %zu/16", rounds);
    log(stdout, "Threshold probability: 2^%f\n", pthresh);
 
    // Init PRNG
    unsigned int seed = time(NULL);
    srand(seed);
    do
    {
        // Reseed so we can deterministically perform a propagation
        seed = rand(); // csprng_rand();
        srand(seed);

        // Build an input differential, print it out
        uint8_t *diff = make_input_diff(rounds, pthresh);
        std::pair<size_t, size_t> result = propagate(diff, rounds, pthresh);
        if (result.first)
        {
            log(stdout, "Fingerprint: %02x%02x%02x%02x%02x%02x%02x%02x",
                    diff[0], diff[1], diff[2], diff[3],
                    diff[4], diff[5], diff[6], diff[7]);
            log(stdout, "Seed: %u", seed);
            log(stdout, "Trails: %zu/%zu\n", result.first, result.second);
        }
    }
    while (1);
    return 0;
}
#endif // __MAIN
