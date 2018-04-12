#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#include "sha2.h"
#include "maw32.h"

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
        "ALGO:\n"
        "  sha256, maw32");
}

// Simple application to run various types of hash algorithms
int main(int argc, char **argv)
{
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
  
    // Get the option
    char *opt  = to_lower_str(argv[1]);
    // This is the only opt which does not have an algorithm specified; deal with it early
    ASSERT(strcmp(opt, "help")); 

    // Get the algorithm
    ASSERT(argc >= 3, puts("Missing hash function"));
    char *algo = to_lower_str(argv[2]), hashbuf[256];
    size_t algo_blocksize, algo_digestsize;

    // Decide which hash function is being used
    char *(*hash)(const uint8_t *, size_t, size_t, char *) = NULL;
    if (!strcmp(algo, "sha256"))
    {
        hash            = sha256_hash;
        algo_blocksize  = 64; 
        algo_digestsize = 32;
    }
    else if (!strcmp(algo, "maw32"))
    {
        hash            = maw32_hash;
        algo_blocksize  = 8;
        algo_digestsize = 4;
    }
    else ASSERT(0, printf("Unknown hash function '%s'\n", algo));

    // We no longer need to access argv[0] (fname), argv[1] (opt), or argv[2] (algo)
    // Shuffle along to make later logic simpler
    argv += 3;
    argc -= 3;

    // Switch on opt
    if (!strcmp(opt, "sample"))
    {
        ASSERT(argc == 3, printf("Option 'sample' requires 3 arguments, got %d\n", argc));
        int min, max, n;
        ASSERT(parse_uint32(argv[0], &min), puts("Argument 'min' was not an integer"));
        ASSERT(parse_uint32(argv[1], &max), puts("Argument 'max' was not an integer"));
        ASSERT(parse_uint32(argv[2], &n), puts("Arugment 'n' was not an integer"));
        const int mod = 1 + max - min;
        
        // Allocate a buffer large enough for any of our strings; we will be reusing this
        uint8_t *buf = (uint8_t *) calloc(max, sizeof(uint8_t));

        // Start sampling
        while (n-- > 0)
        {
            // Pick a length in [min, max], and populate it with random data in [0, 256]
            int len = min + (rand() % mod);
            for (size_t idx = 0; idx < len; idx++) { buf[idx] = rand(); }
            
            // Compute and print the hash
            printf("%s - 0x", hash(buf, len, -1, hashbuf));
            for (size_t idx = 0; idx < len; idx++) { printf("%02x", buf[idx]); }
            printf("\n");
        }
        free(buf);
        return 0;
    }
    else if (!strcmp(opt, "iterate"))
    {
        ASSERT(argc == 1, printf("Option 'iterate' requires 1 argument, got %d\n", argc));
        int n;
        ASSERT(parse_uint32(argv[0], &n), puts("Argument 'n' was not an integer"));
        // Allocate a buffer large enough to store any of our strings; we will be reusing this
        uint8_t *buf = (uint8_t *) calloc(n, sizeof(uint8_t)); 
        // and fill it with zeroes
        memset(buf, 0, n);

        // Start iterating
        while (1)
        {
            // Hash whatever is in buf
            printf("%s - 0x", hash(buf, n, -1, hashbuf));
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
        j
        for (size_t idx = 0; idx < argc; idx++)
        {
            char *s = argv[idx];
            printf("%s - \"%s\"\n", hash(s, strlen(s), -1, hashbuf), s);
        }
        return 0;
    }
    else if (!strcmp(opt, "diff"))
    {
        ASSERT(argc == 2, printf("Operation 'diff' requires 2 arguments, got %d\n", argc));
        int rounds;
        ASSERT(parse_uint32(argv[0], &rounds), puts("Argument 'rounds' was not an integer"));
        // No need to free these; program will only be terminated with ^C.
        uint8_t *diff    = malloc(algo_blocksize),
                *input1  = malloc(algo_blocksize),
                *input2  = malloc(algo_blocksize),
                *digest1 = malloc(2*algo_digestsize+1),
                *digest2 = malloc(2*algo_digestsize+1);
        for (int i = 0; i < algo_blocksize; i++)
        {
            ASSERT(parse_uint8(argv[1] + 2*(i+1), diff+i), printf("Could not parse hex byte at index %d\n", i));
        }

        while (1)
        {
            // Produce two inputs with a difference of diff
            for (int i = 0; i < algo_blocksize; i++) input2[i] = diff[i] ^ (input1[i] = rand());

            // Compute and compare their hashes
            hash(input1, algo_blocksize, rounds, digest1);
            hash(input2, algo_blocksize, rounds, digest2);
            // If a collision is detected, print the inputs + their colliding hash
            if (!strcmp(digest1, digest2))
            {
                for (int i = 0; i < algo_blocksize; i++) printf("%02x", input1[i]);
                printf(" - ");
                for (int i = 0; i < algo_blocksize; i++) printf("%02x", input2[i]);
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
