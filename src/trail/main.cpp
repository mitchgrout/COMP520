// Entry point to maw_trail

#ifndef __MAIN
#define __MAIN

#include <stdlib.h>
#include <stdint.h>
#include "utils.cpp"
#include "maw32_utils.cpp"
#include "maw32_trail.cpp"
#include <vector>
#include <tuple>
using namespace std;

// A representation of a gene for our genetic algorithm
typedef struct
{
    uint8_t diff[16];       // Input difference used
    size_t zero_trails;     // How many output differences of zero observed
    size_t total_trails;    // How many total outputs observed
} gene_t;

// Determine if a gene is 'alive', i.e. usable for breeding
static inline int is_alive(gene_t gene)
{
    return gene.total_trails;
}

// Kill a gene, i.e. make it unusable for breeding
static inline void kill_gene(gene_t *gene)
{
    memset(gene, 0, sizeof(gene_t));
}

// Determine how fit a gene is
static inline double get_fitness(gene_t gene)
{
    // Fitness of a dead gene is zero, otherwise it is the fraction of of zero
    // trails to all observed trails
    if (!is_alive(gene)) return 0.0f;
    else return (gene.zero_trails * 1.0f) / (gene.total_trails * 1.0f);
}

// Pretty-print a gene
static inline void print_gene(gene_t gene)
{
    log(stdout, "(Fingerprint: %02x%02x%02x%02x%02x%02x%02x%02x, Fitness: %lf)", 
            gene.diff[0], gene.diff[1], gene.diff[2], gene.diff[3],
            gene.diff[4], gene.diff[5], gene.diff[6], gene.diff[7],
            get_fitness(gene));
}

// Choose a random gene from the pool, weighted by their fitness
static size_t dice(gene_t *pool, size_t len)
{
    // Start by computing the total fitness
    double total_weight = 0.0f;
    for (size_t idx = 0; idx < len; idx++) total_weight += get_fitness(pool[idx]);
    
    // Roll a random number in [0, total_weight]
    double result = (rand() * total_weight) / RAND_MAX;
    
    // Subtract the fitness of each gene from the total weight
    // When we get <=0, select that gene
    for (size_t idx = 0; idx < len; idx++) if (is_alive(pool[idx]))
    {
        result -= get_fitness(pool[idx]);
        if (result <= 0) return idx;
    }
    return -1;
}

// Compute a 32-bit true-random number
static inline uint32_t csprng_rand32()
{
    static FILE *handle = fopen("/dev/urandom", "r");
    return (getc(handle) << 24) |
           (getc(handle) << 16) |
           (getc(handle) <<  8) |
           (getc(handle) <<  0) ;
}

// min(a,b)
static inline size_t max(size_t a, size_t b) { return a < b? b : a; }

// Determine if an input differential is viable
static inline bool is_viable(uint8_t *W, const size_t rounds, const float l2pthresh,
                             size_t t, size_t ctr)
{
    const size_t x = 4, y = 1;  // (y/x) expansions must be zero
    if (t >= rounds) return x*ctr >= y*(max(8, rounds) - 8);
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

// Create an input differential randomly.
void make_input_diff(uint8_t *sched, size_t rounds, float l2pthresh)
{
    // Zero out the first 4 words
    memset(sched, 0, 4);
    // Randomly assign differences for last four words, and check viability
    do for (int idx = 4; idx < 8; idx++) sched[idx] = rand() & 0xff;
    while (!is_viable(sched, rounds, l2pthresh, 8, 0)) ;
}

// Cross over two message differences
void cross(uint8_t *dest, const uint8_t *left, const uint8_t *right, const size_t mid)
{
    // masks[n] gives the mask for top n bits
    // ~masks[n] gives mask for lower (8-n) bits
    static const uint8_t masks[] = 
    { 0x00, 0x80, 0xc0, 0xe0, 0xf0, 0xf8, 0xfc, 0xfe, 0xff };

    // Byte midpoint
    uint8_t byte_mid = mid / 8;
    // Bit split
    uint8_t bit_split = mid % 8;

    // Copy a chunk of left
    memcpy(dest, left, byte_mid);
    if (bit_split)
    {
        // There is a little bit of left remaining
        dest[byte_mid] = 
            (left[byte_mid]  & masks[bit_split]) | 
            (right[byte_mid] & ~masks[bit_split]);
        memcpy(dest+byte_mid+1, right+byte_mid+1, 8-byte_mid-1);
    }
    else
    {
        // All of left has been copied
        memcpy(dest+byte_mid, right+byte_mid, 8-byte_mid);
    }
}

// Entry point
int main(int argc, char **argv)
{
    // Runtime consts
    const float pthresh = -3.000000f;
    const size_t rounds = 8;
    const char *key_fname = "/Scratch/key-file-3.000000.bin",
               *add_fname = "/Scratch/add-file-3.000000.bin",
               *maj_fname = "/Scratch/maj-file-3.000000.bin";

    // Setup
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
 
    // Reseed PRNG
    srand(csprng_rand32());

    // Gather enough differentials to use for genetic algorithms
    // Does not have to be on the stack hence static
    #define POOL_SIZE 32
    static gene_t pool[POOL_SIZE];
    for (int idx = 0; idx < POOL_SIZE; )
    {
        // Build a differential, print it out
        make_input_diff(pool[idx].diff, rounds, pthresh);
        std::pair<size_t, size_t> result = propagate(pool[idx].diff, rounds, pthresh);
        if (result.first)
        {
            pool[idx].zero_trails  = result.first;
            pool[idx].total_trails = result.second;
            print_gene(pool[idx]);
            idx++;
        }
    }
    log(stdout, "");

    // We now have a full gene pool, begin breeding
    for (size_t pool_num = 1; ; pool_num++)
    {
        // Create a copy; we will destroy this to rebuild the real pool
        static gene_t pool_copy[POOL_SIZE];
        memset(pool_copy, 0, sizeof(pool_copy));

        // Allow half the pool to live; this is done by using `dice` to
        // statistically pick the best solutions
        size_t idx = 0;
        for ( ; idx < POOL_SIZE/2; )
        {
            // Pick an index of a living gene
            size_t survivor_idx = dice(pool, POOL_SIZE);
            // Copy it to the pool
            pool_copy[idx++] = pool[survivor_idx];
            // and remove the gene from the original pool
            kill_gene(&pool[survivor_idx]);
        }
        // Copy the survivors into the real pool
        memcpy(pool, pool_copy, sizeof(pool_copy));
        // and begin breeding new genes
        for ( ; idx < POOL_SIZE; idx++)
        {
            // Note that we have this loop to prevent getting stuck with 'bad'
            // choices for parent genes
            while (1)
            {
                // 1/4 to mutate
                if (!(rand() % 4))
                {
                    // Pick a random gene from our survivors
                    size_t parent_idx = dice(pool, POOL_SIZE/2);
                    // Copy the diff into new_gene
                    memcpy(pool[idx].diff, pool_copy[parent_idx].diff, 8*sizeof(uint8_t));
                    // Flip a random bit in the dense section
                    int bit_idx = 32 + (rand() % 32);
                    pool[idx].diff[bit_idx/8] ^= 1 << (8 - (bit_idx % 8));
                }
                else
                {
                    // Pick two random disinct genes
                    size_t parent1_idx = dice(pool, POOL_SIZE/2);
                    size_t parent2_idx;
                    do parent2_idx = dice(pool,POOL_SIZE/2); while(parent1_idx == parent2_idx);
                    // Pick a midpoint
                    int mid = 32 + (rand() % 32);
                    // and cross
                    cross(pool[idx].diff, 
                          pool_copy[parent1_idx].diff,
                          pool_copy[parent2_idx].diff, 
                          mid);
                }
                // Try to propagate it
                std::pair<size_t, size_t> result = propagate(pool[idx].diff, rounds, pthresh);
                if (result.first)
                {
                    pool[idx].zero_trails  = result.first;
                    pool[idx].total_trails = result.second;
                    break;
                }
                // otherwise continue
            }
        }
        // Everything has been repopulated.
        log(stdout, "Population %zu bred.", pool_num);
        gene_t *best = &pool[0];
        for (int i = 0; i < POOL_SIZE; i++)
        {
            print_gene(pool[i]);
            if (get_fitness(pool[i]) > get_fitness(*best)) best = &pool[i];
        }
        log(stdout, "Current best:");
        print_gene(*best);
        log(stdout, "");
    }
    return 0;
}
#endif // __MAIN
