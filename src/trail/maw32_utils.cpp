// Utility functions associated with the MAW32 hash function

#ifndef __MAW_UTILS
#define __MAW_UTILS

#include <stddef.h>
#include <stdint.h>

// Consts:

// Input block size
#define MAW32_BLOCK_SIZE  64
// Output digest size
#define MAW32_DIGEST_SIZE 32


// MAW32 internal functions:

// Right-rotate
static inline uint8_t rotr(const uint8_t x, const uint8_t n)
{
    return (x >> n) | (x << (8 - n)); 
}

// Major (major-parity)
static inline uint8_t maj(const uint8_t x, const uint8_t y, const uint8_t z)
{
    return (x & y) ^ (x & z) ^ (y & z);
}

// Sigma0 (non-truncating)
static inline uint8_t sigma0(const uint8_t x)
{
    return rotr(x, 2) ^ rotr(x, 3) ^ rotr(x, 5);
}

// Sigma1 (truncating)
static inline uint8_t sigma1(const uint8_t x)
{
    return rotr(x, 1) ^ rotr(x, 4) ^ (x >> 3);
}

// Addition (modulo 256)
static inline uint8_t add(const uint8_t x, const uint8_t y)
{
    return x + y;
}

// MAW32 internal function differences:

// Major
static inline uint8_t maj_diff(const uint8_t x,   const uint8_t y,   const uint8_t z,
                               const uint8_t d_x, const uint8_t d_y, const uint8_t d_z)
{
    return maj(x, y, z) ^ maj(x^d_x, y^d_y, z^d_z);
}

// Addition
static inline uint8_t add_diff(const uint8_t x,   const uint8_t y, 
                               const uint8_t d_x, const uint8_t d_y)
{
    return add(x, y) ^ add(x ^ d_x, y ^ d_y);
}

// Keymix (const-addition)
static inline uint8_t keymix_diff(const uint8_t x, const uint8_t d_x, const size_t round)
{
    const uint8_t K[16] = 
    {
        0xb7, 0xe1, 0x51, 0x62, 0x8a, 0xed, 0x2a, 0x6a,
        0xbf, 0x71, 0x58, 0x80, 0x9c, 0xf4, 0xf3, 0xc7
    };
    return add(x, K[round]) ^ add(x ^ d_x, K[round]);
}
#endif // __MAW_UTILS
