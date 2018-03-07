// General utilities

#include <stdint.h>     // uint32_t, uint64_t
#include <byteswap.h>   // __bswap_32, __bswap_64

// Determine if the system is little-endian
static int is_little_endian()
{
    const union
    {
        long one;
        char little;
    } is_endian = { 1 };
    return is_endian.little;
}

// Left-rotate a 32-bit word
static inline uint32_t rotl32(uint32_t x, uint8_t n)
{
    return (x << n) | (x >> (32 - n));
}

// Right-rotate a 32-bit word
static inline uint32_t rotr32(uint32_t x, uint8_t n)
{
    return (x >> n) | (x << (32 - n));
}

// Left-rotate a 64-bit word
static inline uint32_t rotl64(uint32_t x, uint8_t n) 
{
    return (x << n) | (x >> (64 - n)); 
}

// Right-rotate a 64-bit word
static inline uint32_t rotr64(uint32_t x, uint8_t n)
{
    return (x >> n) | (x << (64 - n));
}

// Convert a 32-bit word to big endian, from host endian
static inline uint32_t h2b32(uint32_t x)
{
    if (is_little_endian())
    {
        return __bswap_32(x);
    }
    else
    {
        return x;
    }
}

// Convert a 32-bit word to host endian, from big endian
static inline uint32_t b2h32(uint32_t x)
{
    return h2b32(x);
}

// Convert a 64-bit word to big endian, from host endian
static inline uint64_t h2b64(uint64_t x)
{
    if (is_little_endian())
    {
        return __bswap_64(x);
    }
    else
    {
        return x;
    }
}

// Convert a 64-bit word to host endian, from big endian
static inline uint64_t b2h64(uint64_t x)
{
    return h2b64(x);
}
