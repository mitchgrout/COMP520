// Implementation of the MAW32 algorithm.
// Description:
// MAW32 is an experimental  hash function designed by Mitchell Grout while
// studying at the University of Waikato. Due to the small block size, digest
// size, and rounds, this is NOT recommended for use in security; instead, this
// is recommended for simple cryptanalysis of hash functions.
// Specification of domain/codomain of the compression function:
// - Block size:       64 bits
// - Digest size:      32 bits
// - Number of rounds: 16
//
// MAW32 uses the Merkle-Damgard construction to allow variable-length block
// sizes. The padding is performed as follows:
// - Append a 1-bit to the input
// - Append 0-bits until the input length (in bits) is congruent to 32 (mod 64)
// - Append the length of the initial input in bits as a 32-bit int (big endian) to the input
//
// We define the following functions:
// - rotr(x, n) = (x >> n) | (x << (8 - n))
// - maj(x,y,z) = (x & y) ^ (x & z) ^ (y & z)
// - sigma0(x)  = rotr(x, 2) ^ rotr(x, 3) ^ rotr(x, 5)
// - sigma1(x)  = rotr(x, 1) ^ rotr(x, 4) ^ (x >> 3)
//
// The compression function is as follows:
// - Set up the key schedule as 16 uint8_t values, using the fractional expansion of e (A170873)
//      uint8_t K[16] = 
//      {
//          0xb7, 0xe1, 0x51, 0x62, 
//          0x8a, 0xed, 0x2a, 0x6a,
//          0xbf, 0x71, 0x58, 0x80,
//          0x9c, 0xf4, 0xf3, 0xc7
//      };
// - Set up the IV as 4 uint8_t values, using the fractioal expansion of pi (A062964)
//      uint8_t H[4] = 
//      {
//          0x24, 0x3f, 0x6a, 0x88
//      };
// - Pad the input as described above
// - For blocks M_i, from i = 1 .. N:
//   {
//      Set up the message schedule, W_t, from t = 1 .. 16:
//      {
//          W_t := M[t]                                     for t = 1 .. 8
//          W_t := sigma0(W[t-3]) + W[t-4] + sigma1(W[t-8]) for t = 9 .. 16
//      }
//
//      let a = H[0], b = H[1], c = H[2], d = H[3]
//      for t = 1 .. 16:
//      {
//          t1 := d + sigma1(b) + K[t] + W[t]
//          t2 := sigma0(a) + maj(a, b, c)
//          d = c
//          c = b + t1
//          b = a
//          a = t1 + t2
//      }
//      let H[0] += a, H[1] += b, H[2] += c, H[3] += d
//   }
//   H is now the MAW32 hash of the input

#include <stdio.h>  // size_t
#include <stdint.h> // uint8_t

// Compute the MAW32 hash of the input.
// Params:
// - ptr: A non-null pointer to an array of data to hash
// - len: The length of ptr in bytes
// - buf: A block of at least 16+1 bytes in which to store the hash; if null, 
//        a local buffer will be used
// Returns: The string representation of the 32-bit digest. If buf is NULL,
//          this will be a local buffer, otherwise buf will be returned. If
//          ptr is NULL, then NULL is returned
char *maw32_hash(const uint8_t *ptr, size_t len, char *buf);
