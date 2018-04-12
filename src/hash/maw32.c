// Implementation of the MAW32 algorithm.

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include "util.h"

#define LOG

// Utilities for MAW32
#define MAW32_BLOCK_SIZE  8 // 64 bits => 8 bytes for every block
#define MAW32_DIGEST_SIZE 4 // 32 bits => 4 bytes for the digest
static inline uint8_t rotr(uint8_t x, uint8_t n) { return (x >> n) | (x << (8 - n)); }
static inline uint8_t maj(uint8_t x, uint8_t y, uint8_t z) { return (x & y) ^ (x & z) ^ (y & z); }
static inline uint8_t sigma0(uint8_t x) { return rotr(x, 2) ^ rotr(x, 3) ^ rotr(x, 5); }
static inline uint8_t sigma1(uint8_t x) { return rotr(x, 1) ^ rotr(x, 4) ^ (x >> 3); }

// Fetch the next block from the input; if padding is required, this will be
// emplaced in a local buffer. Repeat calls to next_block with the same pointer
// will return sequential blocks until the final block has been consumed; at
// this point, NULL will be returned. If a different pointer is passed, then
// the first block of the new pointer is fetched, and the process restarts
static uint8_t *next_block(const uint8_t *ptr, size_t len)
{
    static const uint8_t *lastptr = NULL;       // Start of the first block of the input we are chunking
    static size_t idx = 0;                      // Index of the next block
    static uint8_t local[MAW32_BLOCK_SIZE];     // Local buffer for storing partial blocks
    static int needs_empty_block = 0;           // Whether or not we have a case where the final block is empty

    // New input
    if (lastptr != ptr)
    {
        lastptr = ptr;
        idx = 0;
        needs_empty_block = 0;
    }

    // Edge case: we hit case #3 last time, and need to produce a block with all zeroes + length
    if (needs_empty_block)
    {
        memset(local, 0, MAW32_BLOCK_SIZE - 4);
        *(uint32_t *)(local + MAW32_BLOCK_SIZE - 4) = h2b32(len * 8);
        idx += MAW32_BLOCK_SIZE;
        needs_empty_block = 0;
        return local;
    }

    // Marker that we have exhausted all the data in ptr
    if (idx > len)
    {
        lastptr = NULL;
        return NULL;
    }
    
    // Pointer to the next [potentially partial] block
    const uint8_t *block = ptr + idx;
    // How many bytes of actual data remain
    size_t partial = len - idx;

    // Do we have a full block to send
    if (partial >= MAW32_BLOCK_SIZE)
    {
        // Mark that this block has been consumed
        idx += MAW32_BLOCK_SIZE;
        return block;
    }
    // A partial block which can fit at least 1 byte padding + length
    else if (partial <= MAW32_BLOCK_SIZE - 1 - 4)
    {
        // Copy the partial block into local, then append bits as per the standard
        memcpy(local, block, partial);
        local[partial++] = 0x80;
        while (partial < MAW32_BLOCK_SIZE - 4) { local[partial++] = 0x00; }
        // Append the length (in bits) as a big-endian 64-bit int
        *(uint32_t *)(local + partial) = h2b32(len * 8);
        idx += MAW32_BLOCK_SIZE;
        return local;
    }
    // A partial block which requires the creation of an extra block afterwards (56 to 63 bytes)
    else
    {
        // Copy the partial block into local, then append bits as per the standard
        memcpy(local, block, partial);
        local[partial++] = 0x80;
        while (partial < MAW32_BLOCK_SIZE) { local[partial++] = 0x00; }
        idx += MAW32_BLOCK_SIZE;
        needs_empty_block = 1;
        return local;
    }
}

// Compute the MAW32 hash of the input.
// Params:
// - ptr: A non-null pointer to an array of data to hash
// - len: The length of ptr in bytes
// - buf: A block of at least 16+1 bytes in which to store the hash; if null, 
//        a local buffer will be used
// Returns: The string representation of the 32-bit digest. If buf is NULL,
//          this will be a local buffer, otherwise buf will be returned. If
//          ptr is NULL, then NULL is returned
char *maw32_hash(const uint8_t *ptr, size_t len, char *buf)
{
    // No pointer
    if (!ptr) { return NULL; }

    // Set up the IV
    uint8_t H[4] = 
    {
        0x24, 0x3f, 0x6a, 0x88
    };

    // Set up the constants
    uint8_t K[16] = 
    {
        0xb7, 0xe1, 0x51, 0x62, 
        0x8a, 0xed, 0x2a, 0x6a,
        0xbf, 0x71, 0x58, 0x80,
        0x9c, 0xf4, 0xf3, 0xc7
    };

    // Current block
    uint8_t *M = NULL;
    int block = 0;
    while ((M = next_block(ptr, len)) != NULL)
    {
        // Registers
        uint8_t a = H[0], b = H[1], c = H[2], d = H[3];
        // Message schedule
        uint8_t W[16];
    
        // Transform this block
        for (int t = 0; t < 16; t++)
        {
            // Set up the message schedule for this round
            if (t < 8) { W[t] = M[t]; }
            else { W[t] = sigma0(W[t-3]) + W[t-4] + sigma1(W[t-8]); }
            uint8_t t1 = d + sigma1(b) + K[t] + W[t];
            uint8_t t2 = sigma0(a) + maj(a, b, c);
            d = c;
            c = b + t1;
            b = a;
            a = t1 + t2;
        }
        
        // Update H
        H[0] += a; H[1] += b; H[2] += c; H[3] += d;
    }

    // Copy H into buf, using a local buffer if necessary
    static char local[2*MAW32_DIGEST_SIZE+1];
    if (!buf) { buf = local; }
    sprintf(buf, "%02x%02x%02x%02x", H[0], H[1], H[2], H[3]);
    return buf;
}
