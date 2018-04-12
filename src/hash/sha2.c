// Implementation of SHA-2 algorithms, specifically:
// - SHA-256

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "util.h"

// SHA-256:
#define SHA256_BLOCK_SIZE  64   // 512 bits => 64 bytes for every block
#define SHA256_DIGEST_SIZE 32   // 256 bits => 32 bytes for the digest
static inline uint32_t ch32(uint32_t x, uint32_t y, uint32_t z)  { return (x & y) ^ ((~x) & z);        }
static inline uint32_t maj32(uint32_t x, uint32_t y, uint32_t z) { return (x & y) ^ (x & z) ^ (y & z); }
static inline uint32_t Sigma0_256(uint32_t x) { return rotr32(x, 2)  ^ rotr32(x, 13) ^ rotr32(x, 22);  }
static inline uint32_t Sigma1_256(uint32_t x) { return rotr32(x, 6)  ^ rotr32(x, 11) ^ rotr32(x, 25);  }
static inline uint32_t sigma0_256(uint32_t x) { return rotr32(x, 7)  ^ rotr32(x, 18) ^ (x >> 3);       }
static inline uint32_t sigma1_256(uint32_t x) { return rotr32(x, 17) ^ rotr32(x, 19) ^ (x >> 10);      }

// SHA-512:
#define SHA512_BLOCK_SIZE  128  // 1024 bits => 128 bytes for every block
#define SHA512_DIGEST_SIZE 64   // 512 bits  => 64  bytes for every block
static inline uint32_t ch64(uint32_t x, uint32_t y, uint32_t z)  { return (x & y) ^ ((~x) & z);        }
static inline uint32_t maj64(uint32_t x, uint32_t y, uint32_t z) { return (x & y) ^ (x & z) ^ (y & z); }
static inline uint32_t Sigma0_512(uint32_t x) { return rotr64(x, 28) ^ rotr64(x, 34) ^ rotr64(x, 39);  }
static inline uint32_t Sigma1_512(uint32_t x) { return rotr64(x, 14) ^ rotr64(x, 18) ^ rotr64(x, 41);  }
static inline uint32_t sigma0_512(uint32_t x) { return rotr64(x, 1)  ^ rotr64(x, 8)  ^ (x >> 7);       }
static inline uint32_t sigma1_512(uint32_t x) { return rotr64(x, 19) ^ rotr64(x, 61) ^ (x >> 6);       }

// Fetch the next block from the input; if padding is required, this will be
// emplaced in a local buffer. Repeat calls to next_block with the same pointer
// will return sequential blocks until the final block has been consumed; at
// this point, NULL will be returned. If a different pointer is passed, then
// the first block of the new pointer is fetched, and the process restarts
static uint8_t *next_block(const uint8_t *ptr, size_t len)
{
    static const uint8_t *lastptr = NULL;       // Start of the first block of the input we are chunking
    static size_t idx = 0;                      // Index of the next block
    static uint8_t local[SHA256_BLOCK_SIZE];    // Local buffer for storing partial blocks
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
        memset(local, 0, SHA256_BLOCK_SIZE - 8);
        *(uint64_t *)(local + SHA256_BLOCK_SIZE - 8) = h2b64(len * 8);
        idx += SHA256_BLOCK_SIZE;
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
    if (partial >= SHA256_BLOCK_SIZE)
    {
        // Mark that this block has been consumed
        idx += SHA256_BLOCK_SIZE;
        return block;
    }
    // A partial block which can fit at least 1 byte padding + length
    else if (partial <= SHA256_BLOCK_SIZE - 1 - 8)
    {
        // Copy the partial block into local, then append bits as per the standard
        memcpy(local, block, partial);
        local[partial++] = 0x80;
        while (partial < SHA256_BLOCK_SIZE - 8) { local[partial++] = 0x00; }
        // Append the length (in bits) as a big-endian 64-bit int
        *(uint64_t *)(local + partial) = h2b64(len * 8);
        idx += SHA256_BLOCK_SIZE;
        return local;
    }
    // A partial block which requires the creation of an extra block afterwards (56 to 63 bytes)
    else
    {
        // Copy the partial block into local, then append bits as per the standard
        memcpy(local, block, partial);
        local[partial++] = 0x80;
        while (partial < SHA256_BLOCK_SIZE) { local[partial++] = 0x00; }
        idx += SHA256_BLOCK_SIZE;
        needs_empty_block = 1;
        return local;
    }
}

// Compute the SHA256 hash of the input.
// Params:
// - ptr: A non-null pointer to an array of data to hash
// - len: The length of ptr in bytes
// - buf: A block of at least 32+1 bytes in which to store the hash; if null,
//        a local buffer will be used
// Returns: The string representation of the 256-bit digest. If buf is NULL,
//          this will be a local buffer, otherwise buf will be returned. If
//          ptr is NULL, then NULL is returned
char *sha256_hash(const uint8_t *ptr, size_t len, char *buf)
{
    // No pointer
    if (!ptr) { return NULL; }

    // Set up the IV
    uint32_t H[8] = 
    {
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
        0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19 
    };

    // Set up the constants
    uint32_t K[64] = 
    {
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
    };

    // Current block
    uint8_t *M = NULL;
    int block = 0;
    while ((M = next_block(ptr, len)) != NULL)
    {
        // Registers
        uint32_t a = H[0], b = H[1], c = H[2], d = H[3],
                 e = H[4], f = H[5], g = H[6], h = H[7];
        // Message schedule
        uint32_t W[64];

        // Transform this block
        for (int t = 0; t < 64; t++)
        {
            // Set up the message schedule for this round
            if (t < 16) { W[t] = b2h32(((uint32_t *)M)[t]); }
            else { W[t] = sigma1_256(W[t-2]) + W[t-7] + sigma0_256(W[t-15]) + W[t-16]; }
            
            // Transform the block
            uint32_t T1 = h + Sigma1_256(e) + ch32(e, f, g) + K[t] + W[t];
            uint32_t T2 = Sigma0_256(a) + maj32(a, b, c);
            h = g; g = f; f = e; e = d + T1;
            d = c; c = b; b = a; a = T1 + T2;          
        }
 
        // Update H
        H[0] += a; H[1] += b; H[2] += c; H[3] += d;
        H[4] += e; H[5] += f; H[6] += g; H[7] += h; 
    }

    // Copy H into buf, using a local buffer if necessary
    char local[2*SHA256_DIGEST_SIZE+1];
    if (!buf) { buf = local; }
    sprintf(buf, "%08x%08x%08x%08x%08x%08x%08x%08x",
            H[0], H[1], H[2], H[3], H[4], H[5], H[6], H[7]);
    return buf;
}
