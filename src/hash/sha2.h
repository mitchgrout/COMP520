// Headers for SHA-2 algorithms, specifically:
// - SHA-256

#include <stdio.h>  // size_t
#include <stdint.h> // uint8_t

// Compute the SHA256 hash of the input.
// Params:
// - ptr: A non-null pointer to an array of data to hash
// - len: The length of ptr in bytes
// - buf: A block of at least 32+1 bytes in which to store the hash; if null,
//        a local buffer will be used
// Returns: The string representation of the 256-bit digest. If buf is NULL,
//          this will be a local buffer, otherwise buf will be returned. If
//          ptr is NULL, then NULL is returned
char *sha256_hash(const uint8_t *ptr, size_t len, char *buf);
