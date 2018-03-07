#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#include "sha2.h"
#include "maw32.h"

// Mutate a string in-place, making it lowercase
char *tolowerstr(char *str)
{
    // Don't try if the input is invalid
    if (!str)
    {
        return str;
    }

    char *sent = str;
    while (*sent)
    {
        if (*sent >= 'A' && *sent <= 'Z') *sent += 'a' - 'A';
        sent++;
    }
    return str;
}

// Determine if the given string is a decimal string
bool isdigitstr(char *str)
{
    // Empty strings are not digit strings
    if (!str || !*str)
    {
        return false;
    }

    char *sent = str;
    while (*sent)
    {
        if (*sent < '0' || *sent > '9') return false;
        sent++;
    }
    return true;
}

// Parse an integer, first checking that it indeed is an integer string
bool parse_int32(char *str, uint32_t *result)
{
    if (!isdigitstr(str)) return false;
    *result = atoi(str);
    return true;
}

// Show the usage for the program
void show_usage()
{
    puts(
        "USAGE: hash [OPT] [...]\n" 
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
        "  help - Show this message\n" 
        "\n"
        "ALGO:\n"
        "  sha256, maw32");
}

// Simple application to run various types of hash algorithms
int main(int argc, char **argv)
{
    // Need at least 2 arguments: { "hash", "[OPT]", ... }
    if (argc < 2)
    {
        show_usage();
        return 1;
    }
   
    // Get the option
    char *opt  = tolowerstr(argv[1]);
    // This is the only opt which does not have an algorithm specified; deal with it early
    if (!strcmp(opt, "help"))
    {
        show_usage();
        return 0;
    }

    // Get the algorithm
    if (argc < 3)
    {
        printf("Missing algorithm");
        show_usage();
        return 1;
    }
    char *algo = tolowerstr(argv[2]);
    char hashbuf[2048];

    // Decide which hash function is being used
    char *(*hash)(const uint8_t *, size_t, char *) = NULL;
    if (!strcmp(algo, "sha256"))     hash = sha256_hash;
    else if (!strcmp(algo, "maw32")) hash = maw32_hash;
    else
    {
        printf("Unknown hash function \"%s\"\n", algo);
        show_usage();
        return 1;
    }

    // Switch on opt
    if (!strcmp(opt, "sample"))
    {
        // sample has 4 arguments
        if (argc != 6)
        {
            printf("Incorrect number of arguments for option \"sample\"; requires 4 arguments\n");
            show_usage();
            return 1;
        }

        // Get our arguments
        int min, max, n;
        if (!parse_int32(argv[3], &min))
        {
            printf("Argument \"min\" was not an integer\n");
            show_usage();
            return 1;
        }
        if (!parse_int32(argv[4], &max))
        {
            printf("Argument \"max\" was not an integer\n");
            show_usage();
            return 1;
        }
        if (!parse_int32(argv[5], &n))
        {
            printf("Argument \"n\" was not an integer\n");
            show_usage();
            return 1;
        }
        const int mod = 1 + max - min;
        
        // Allocate a buffer large enough for any of our strings; we will be reusing this
        uint8_t *buf = (uint8_t *) calloc(max, sizeof(uint8_t));

        // Seed the prng
        srand(time(NULL));

        // Start sampling
        while (n-- > 0)
        {
            // Pick a length in [min, max], and populate it with random data in [0, 256]
            int len = min + (rand() % mod);
            for (size_t idx = 0; idx < len; idx++) { buf[idx] = rand(); }
            
            // Compute and print the hash
            printf("%s - 0x", hash(buf, len, hashbuf));
            for (size_t idx = 0; idx < len; idx++) { printf("%02x", buf[idx]); }
            printf("\n");
        }
        free(buf);
        return 0;
    }
    else if (!strcmp(opt, "iterate"))
    {
        // iterate has 1 argument1
        if (argc != 4)
        {
            printf("Incorrect number of arguments for option \"iterate\"; requires 3 arguments\n");
            show_usage();
            return 1;
        }

        // Get our arguments
        int n;
        if (!parse_int32(argv[3], &n))
        {
            printf("Argument \"n\" was not an integer\n");
            show_usage();
            return 1;
        }
   
        // Allocate a buffer large enough to store any of our strings; we will be reusing this
        uint8_t *buf = (uint8_t *) calloc(n, sizeof(uint8_t)); 
        // and fill it with zeroes
        memset(buf, 0, n);

        // Start iterating
        while (1)
        {
            // Hash whatever is in buf
            printf("%s - 0x", hash(buf, n, hashbuf));
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
        for (size_t idx = 3; idx < argc; idx++)
        {
            char *s = argv[idx];
            printf("%s - \"%s\"\n", hash(s, strlen(s), hashbuf), s);
        }
        return 0;
    }
    else
    {
        printf("Unknown option \"%s\"\n", opt);
        show_usage();
        return 1;
    }
}
