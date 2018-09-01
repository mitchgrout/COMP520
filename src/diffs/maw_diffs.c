#include <stdio.h>
#include <stdint.h>
#include <string.h>

// Utilities for MAW32
uint8_t K[16] = 
{
    0xb7, 0xe1, 0x51, 0x62, 
    0x8a, 0xed, 0x2a, 0x6a,
    0xbf, 0x71, 0x58, 0x80,
    0x9c, 0xf4, 0xf3, 0xc7
};
static inline uint8_t rotr(uint8_t x, uint8_t n) { return (x >> n) | (x << (8 - n)); }
static inline uint8_t maj(uint8_t x, uint8_t y, uint8_t z) { return (x & y) ^ (x & z) ^ (y & z); }
static inline uint8_t sigma0(uint8_t x) { return rotr(x, 2) ^ rotr(x, 3) ^ rotr(x, 5); }
static inline uint8_t sigma1(uint8_t x) { return rotr(x, 1) ^ rotr(x, 4) ^ (x >> 3); }
static inline uint8_t key_mix(uint8_t m, size_t k) { return m + K[k]; }
static inline uint8_t add(uint8_t x, uint8_t y) { return x + y; }

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

// Show the usage for the program
void show_usage()
{
    puts(
        "USAGE: maw_diff [FUNC] [...]\n" 
        "\n" 
        "FUNC:\n"
        "  sigma0 (d_m): Iterate with differential d_m\n"
        "  sigma1 (d_m): Iterate with differential d_m\n"
        "  keymix (k, d_m): Iterate adding round k const with differential d_m\n"
        "  maj (d_x, d_y, d_z): Iterate with three differentials\n"
        "  add (d_x, d_y): Add differentials\n");
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        show_usage();
        return 1;
    }

    char *func = tolowerstr(argv[1]);
    if (!strcmp(func, "sigma0"))
    {
        puts("Differences for sigma0:");
        int d_m;
        sscanf(argv[2], "0x%02x", &d_m);
        int counts[256];
        for (int i = 0; i < 256; i++) counts[i] = 0;
        for (int m = 0; m < 256; m++)
        {
            counts[sigma0(m) ^ sigma0(m ^ d_m)]++;
        }
        for (int d_h = 0; d_h < 256; d_h++) if (counts[d_h])
        {
            printf("0x%02x : %d/256\n", d_h, counts[d_h]);
        }
        return 0;
    }
    else if (!strcmp(func, "sigma1"))
    {
        puts("Differences for sigma1:");
        int d_m;
        sscanf(argv[2], "0x%02x", &d_m);
        int counts[256];
        for (int i = 0; i < 256; i++) counts[i] = 0;
        for (int m = 0; m < 256; m++)
        {
            counts[sigma1(m) ^ sigma1(m ^ d_m)]++;
        }
        for (int d_h = 0; d_h < 256; d_h++) if (counts[d_h])
        {
            printf("0x%02x : %d/256\n", d_h, counts[d_h]);
        }
        return 0;
    }
    else if (!strcmp(func, "keymix"))
    {
        int k;
        sscanf(argv[2], "%d", &k);
        printf("Differences for keymix-%d:\n", k);
        int d_m;
        sscanf(argv[3], "0x%02x", &d_m);
        int counts[256];
        for (int i = 0; i < 256; i++) counts[i] = 0;
        for (int m = 0; m < 256; m++)
        {
            counts[key_mix(m, k) ^ key_mix(m ^ d_m, k)]++;
        }
        for (int d_h = 0; d_h < 256; d_h++) if (counts[d_h])
        {
            printf("0x%02x : %d/256\n", d_h, counts[d_h]);
        }
        return 0;
    }
    else if (!strcmp(func, "maj"))
    {
        puts("Differences for maj:");
        int d_x, d_y, d_z;
        sscanf(argv[2], "0x%02x", &d_x);
        sscanf(argv[3], "0x%02x", &d_y);
        sscanf(argv[4], "0x%02x", &d_z);
        int counts[256];
        for (int i = 0; i < 256; i++) counts[i] = 0;
        for (int x = 0; x < 256; x++)
        for (int y = 0; y < 256; y++)
        for (int z = 0; z < 256; z++)
        {
            counts[maj(x, y, z) ^ maj(x ^ d_x, y ^ d_y, z ^ d_z)]++;
        }
        for (int d_h = 0; d_h < 256; d_h++) if (counts[d_h])
        {
            printf("0x%02x : %d/16777216\n", d_h, counts[d_h]);
        }
        return 0;
    }
    else if (!strcmp(func, "add"))
    {
        puts("Differences for +:");
        int d_x, d_y;
        sscanf(argv[2], "0x%02x", &d_x);
        sscanf(argv[3], "0x%02x", &d_y);
        int counts[256];
        for (int i = 0; i < 256; i++) counts[i] = 0;
        for (int x = 0; x < 256; x++)
        for (int y = 0; y < 256; y++)
        {
            counts[add(x, y) ^ add(x ^ d_x, y ^ d_y)]++;
        }
        for (int d_h = 0; d_h < 256; d_h++) if (counts[d_h])
        {
            printf("0x%02x : %d/65536\n", d_h, counts[d_h]);
        }
        return 0;
    }
    else
    {
        printf("Unknown function '%s'\n", func);
        show_usage();
        return 1;
    }
}
