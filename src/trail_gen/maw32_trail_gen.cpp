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
#include <utility>
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

// Derivatives
static inline uint8_t maj_diff(const uint8_t x,   const uint8_t y,   const uint8_t z,
                               const uint8_t d_x, const uint8_t d_y, const uint8_t d_z)
{
    return maj(x, y, z) ^ maj(x^d_x, y^d_y, z^d_z);
}

static inline uint8_t add_diff(const uint8_t x,   const uint8_t y, 
                               const uint8_t d_x, const uint8_t d_y)
{
    return add(x, y) ^ add(x ^ d_x, y ^ d_y);
}

static inline uint8_t keymix_diff(const uint8_t x, const uint8_t d_x, const size_t round)
{
    const uint8_t K[16] = 
    {
        0xb7, 0xe1, 0x51, 0x62, 0x8a, 0xed, 0x2a, 0x6a,
        0xbf, 0x71, 0x58, 0x80, 0x9c, 0xf4, 0xf3, 0xc7
    };
    return add(x, K[round]) ^ add(x ^ d_x, K[round]);
}

// Util; takes a collection of inputs and frequencies, and returns a 
// collection of inputs and log2 probabilities, all greater than l2pthresh
static inline vector<pair<uint8_t,int8_t>> filter(const map<uint8_t, size_t>& samples, const size_t sample_size, const float l2pthresh)
{
    vector<pair<uint8_t,int8_t>> results;
    for (const auto& elem : samples)
    {
        float prob = log2f(elem.second) - log2f(sample_size);
        if (prob >= l2pthresh) results.push_back(make_pair(elem.first, (uint8_t) floor(prob)));
    }
    return results;
}

// Nonlinear
vector<pair<uint8_t,int8_t>> propagate_keymix(const uint8_t d_x, const size_t round, const float l2pthresh)
{
    const size_t sample_size = 256;
    map<uint8_t, size_t> counts;

    // Sample
    for (int n = 0; n < sample_size; n++)
    {
        uint8_t x = n;
        counts[keymix_diff(x, d_x, round)]++;
    }
    return filter(counts, sample_size, l2pthresh);
}

// Nonlinear
vector<pair<uint8_t,int8_t>> propagate_add(const uint8_t d_x, const uint8_t d_y, const float l2pthresh)
{
    const size_t sample_size = 256*256;
    map<uint8_t, size_t> counts;

    // Sample
    for (int n = 0; n < sample_size; n++)
    {
        uint8_t x = (n >> 8) & 0xff,
                y = (n >> 0) & 0xff;
        counts[add_diff(x, y, d_x, d_y)]++;
    }
    return filter(counts, sample_size, l2pthresh);
}

// Nonlinear
vector<pair<uint8_t,int8_t>> propagate_maj(const uint8_t d_x, const uint8_t d_y, const uint8_t d_z, const float l2pthresh)
{
    const size_t sample_size = 256*256;
    map<uint8_t, size_t> counts;

    // Sample
    for (int n = 0; n < sample_size; n++)
    {
        uint8_t x = rand() & 0xff,
                y = rand() & 0xff,
                z = rand() & 0xff;
        counts[maj_diff(x, y, z, d_x, d_y, d_z)]++;
    }
    return filter(counts, sample_size, l2pthresh);
}

// Entry point
int main(int argc, char **argv)
{
    if (argc != 2)
    {
        puts("Usage: ./gen [probability]");
        return 1;
    }
    float pthresh = atof(argv[1]);
    unsigned int seed = time(NULL);
    srand(seed);
    printf("SEED: %u\n", seed);
    printf("PTHRESH: %f\n", pthresh);

    char key_fname[256],
         add_fname[256],
         maj_fname[256];
    sprintf(key_fname, "/Scratch/key-file%f.bin", pthresh);
    sprintf(add_fname, "/Scratch/add-file%f.bin", pthresh);
    sprintf(maj_fname, "/Scratch/maj-file%f.bin", pthresh);
    FILE *key_file = fopen(key_fname, "w+"),
         *add_file = fopen(add_fname, "w+"),
         *maj_file = fopen(maj_fname, "w+");
 
    // NOTE: File format is as follows
    // [args...],[len], [output_diff],[output_diff_log2prob], ...
    // where each term is a single uint8_t
    // For example, for maj_memo, the output may look like:
    // [x=0],[y=1], [z=2], [len=4], [output=0x0], [prob=-2], [output=0x1], [prob=-2], ...
    puts("Creating keymix memo...");
    for (int round = 0; round < 16; round++) for (int d_x = 0; d_x < 256; d_x++)
    {
        vector<pair<uint8_t,int8_t>> result = propagate_keymix(d_x, round, pthresh);
        putc(d_x & 0xff, key_file);
        putc(round & 0xff, key_file);
        putc(result.size() & 0xff, key_file);
        for (auto out_pair : result) 
        {
            putc(out_pair.first, key_file); 
            putc(out_pair.second, key_file); 
        }
    }
    fclose(key_file);

    puts("Creating add memo...");
    for (int d_x = 0; d_x < 256; d_x++) for (int d_y = 0; d_y < 256; d_y++)
    {
        vector<pair<uint8_t,int8_t>> result = propagate_add(d_x, d_y, pthresh);
        putc(d_x & 0xff, add_file);
        putc(d_y & 0xff, add_file);
        putc(result.size() & 0xff, add_file);
        for (auto out_pair : result) 
        {
            putc(out_pair.first, add_file); 
            putc(out_pair.second, add_file); 
        }
    }
    fclose(add_file);

    puts("Creating maj memo...");
    for (int d_x = 0; d_x < 256; d_x++)
    for (int d_y = 0; d_y < 256; d_y++)
    for (int d_z = 0; d_z < 256; d_z++)
    {
        vector<pair<uint8_t,int8_t>> result = propagate_maj(d_x, d_y, d_z, pthresh);
        putc(d_x & 0xff, maj_file);
        putc(d_y & 0xff, maj_file);
        putc(d_z & 0xff, maj_file);
        putc(result.size() & 0xff, maj_file);
        for (auto out_pair : result) 
        {
            putc(out_pair.first, maj_file);
            putc(out_pair.second, maj_file); 
        }
    }
    fclose(maj_file);
    puts("Done!");
    return 0;
}
