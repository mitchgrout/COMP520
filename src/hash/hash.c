#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#include "sha2.h"
#include "maw32.h"

// Description of a hash algorithm
struct hash_algo
{
    char *name;         // Name (lowercase)
    size_t blk_size;    // Block size in bytes
    size_t dig_size;    // Digest size in bytes
    char *(*hash)(const uint8_t *, size_t, size_t, char *); // Pointer to the function
};

// Const array of supported hash functions
static const struct hash_algo known_algos[] = 
{
    { "maw32", 8, 4, maw32_hash },
    { "sha256", 64, 32, sha256_hash },
};
#define LEN(_arr) (sizeof(_arr)/sizeof((_arr)[0]))

// Mutate a string in-place, making it lowercase
inline char *to_lower_str(char *str)
{
    for (char *sent = str; *sent; sent++) if (*sent >= 'A' && *sent <= 'Z') *sent += 'a' - 'A';
    return str;
}

// Determine if the given string is a decimal string
inline bool is_digit_str(char *str)
{
    for (char *sent = str; *sent; sent++) if (*sent < '0' || *sent > '9') return false;
    return true;
}

// Parse a hex byte
inline bool parse_uint8(char *str, uint8_t *result)
{
    return sscanf(str, "%02hhX", result);
}

// Parse an integer, first checking that it indeed is an integer string
inline bool parse_uint32(char *str, uint32_t *result)
{
    if (!is_digit_str(str)) return false;
    *result = atoi(str);
    return true;
}

// Show the usage for the program
void show_usage()
{
    puts(
        "USAGE: hash [OPT] [ALGO] [...]\n" 
        "\n" 
        "OPT:\n" 
        "  sample (ALGO, min, max, n):\n"  
        "    Randomly sample n bytestrings of length [min, max] inclusive\n"
        "    and compute the ALGO hash of these strings. Statistical\n"
        "    information about ALGO is collected from the results.\n"
        "  iterate (ALGO, n):\n"
        "    Iterate through every bytestring of length n, and compute\n"
        "    the ALGO hash of these strings.\n"
        "  test (ALGO, inp...):\n"
        "    Compute the ALGO hash of every input given in inp..., in order.\n"
        "    The hashes are presented on different lines along with\n"
        "    their inputs.\n"
        "  diff (ALGO, rounds, diff):\n"
        "    Randomly sample inputs, and determine if the input XOR the diff\n"
        "    results in a collision after a given number of rounds.\n"
        "    diff should be expressed as a single hexadecimal value\n"
        "  help ():\n"
        "    Show this message\n" 
        "\n"
        "ALGO:");
    printf("  ");
    for (int i = 0; i < LEN(known_algos); i++)
    {
        printf("%s", known_algos[i].name);
        if (i != LEN(known_algos)-1) printf(", ");
    }
    putchar('\n');
}

// Simple application to run various types of hash algorithms
int main(int argc, char **argv)
{
    // Used to ensure user-supplied data is valid
    #define ASSERT(expr, ...)   \
    if (!(expr))                \
    {                           \
        __VA_ARGS__;            \
        show_usage();           \
        return 1;               \
    }

    // Need at least 2 arguments: { "hash", "[OPT]", ... }
    ASSERT(argc >= 2);
    // Seed the PRNG
    srand(time(NULL));

    // Lowercase argv[1] and argv[2]
    argv[1] = to_lower_str(argv[1]);
    argv[2] = to_lower_str(argv[2]);

    // Get the option
    char *opt  = to_lower_str(argv[1]);
    // This is the only opt which does not have an algorithm specified; deal with it early
    ASSERT(strcmp(opt, "help")); 

    // Decide what hash function is being used
    ASSERT(argc >= 3, puts("Missing hash function"));
    char hashbuf[256];
    struct hash_algo algo = { NULL, -1, -1, NULL };
    to_lower_str(argv[2]);
    for (int i = 0; i < LEN(known_algos); i++) if (!strcmp(argv[2], known_algos[i].name))
    {
        algo = known_algos[i];
        break;
    }
    ASSERT(algo.hash, printf("Unknown hash function '%s'\n", argv[2]));

    // We no longer need to access argv[0] (fname), argv[1] (opt), or argv[2] (algo)
    // Shuffle along to make later logic simpler
    argv += 3;
    argc -= 3;

    // Switch on opt
    if (!strcmp(opt, "sample"))
    {
        ASSERT(argc == 3, printf("Option 'sample' requires 3 arguments, got %d\n", argc));
        uint32_t min, max, n;
        ASSERT(parse_uint32(argv[0], &min), puts("Argument 'min' was not an integer"));
        ASSERT(parse_uint32(argv[1], &max), puts("Argument 'max' was not an integer"));
        ASSERT(parse_uint32(argv[2], &n), puts("Argument 'n' was not an integer"));
        ASSERT(max >= min, puts("Max must be greater than or equal to min"));
        const uint32_t mod = 1 + max - min;
        
        // Allocate a buffer large enough for any of our strings; we will be reusing this
        uint8_t *buf = (uint8_t *) calloc(max, sizeof(uint8_t));

        // Start sampling
        while (n-- > 0)
        {
            // Pick a length in [min, max], and populate it with random data in [0, 256]
            uint32_t len = min + (rand() % mod);
            for (size_t idx = 0; idx < len; idx++) { buf[idx] = rand(); }
            
            // Compute and print the hash
            printf("%s - 0x", algo.hash(buf, len, -1, hashbuf));
            for (size_t idx = 0; idx < len; idx++) { printf("%02x", buf[idx]); }
            printf("\n");
        }
        free(buf);
        return 0;
    }
    else if (!strcmp(opt, "iterate"))
    {
        ASSERT(argc == 1, printf("Option 'iterate' requires 1 argument, got %d\n", argc));
        uint32_t n;
        ASSERT(parse_uint32(argv[0], &n), puts("Argument 'n' was not an integer"));
        // Allocate a buffer large enough to store any of our strings; we will be reusing this
        uint8_t *buf = (uint8_t *) calloc(n, sizeof(uint8_t)); 
        // and fill it with zeroes
        memset(buf, 0, n);

        // Start iterating
        while (1)
        {
            // Hash whatever is in buf
            printf("%s - 0x", algo.hash(buf, n, -1, hashbuf));
            for (size_t idx = 0; idx < n; idx++) { printf("%02x", buf[idx]); }
            printf("\n");

            // Take the successor of the input
            size_t offset;
            for (offset = 0; offset < n; offset++)
            {
                if (buf[n - offset - 1] == 255)
                {
                    buf[n - offset - 1] = 0;
                    continue;
                }
                else
                {
                    buf[n - offset - 1]++;
                    break;
                }
            }
            // buf just wrapped back round to all zeroes; finished
            if (offset == n)
            {
                break;
            }
        }
        free(buf);
        return 0;
    }
    else if (!strcmp(opt, "test"))
    {
        for (size_t idx = 0; idx < argc; idx++)
        {
            char *s = argv[idx];
            printf("%s - \"%s\"\n", algo.hash(s, strlen(s), -1, hashbuf), s);
        }
        return 0;
    }
    else if (!strcmp(opt, "diff"))
    {
        ASSERT(argc == 2, printf("Operation 'diff' requires 2 arguments, got %d\n", argc));
        int rounds;
        ASSERT(parse_uint32(argv[0], &rounds), puts("Argument 'rounds' was not an integer"));
        // No need to free these; program will only be terminated with ^C.
        uint8_t *diff    = malloc(algo.blk_size),
                *input1  = malloc(algo.blk_size),
                *input2  = malloc(algo.blk_size),
                *digest1 = malloc(2*algo.dig_size+1),
                *digest2 = malloc(2*algo.dig_size+1);
        for (int i = 0; i < algo.blk_size; i++)
            ASSERT(parse_uint8(argv[1] + 2*(i+1), diff+i), printf("Could not parse hex byte at index %d\n", i));

        while (1)
        {
            // Produce two inputs with a difference of diff
            for (int i = 0; i < algo.blk_size; i++) input2[i] = diff[i] ^ (input1[i] = rand());

            // Compute and compare their hashes
            algo.hash(input1, algo.blk_size, rounds, digest1);
            algo.hash(input2, algo.blk_size, rounds, digest2);
            // If a collision is detected, print the inputs + their colliding hash
            if (!strcmp(digest1, digest2))
            {
                for (int i = 0; i < algo.blk_size; i++) printf("%02x", input1[i]);
                printf(" - ");
                for (int i = 0; i < algo.blk_size; i++) printf("%02x", input2[i]);
                printf(" => %s\n", digest1);
            }
        }
    }
    else
    {
        printf("Unknown option \"%s\"\n", opt);
        show_usage();
        return 1;
    }
}
