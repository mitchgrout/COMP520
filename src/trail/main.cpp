// Entry point to maw_trail

#ifndef __MAIN
#define __MAIN

#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <semaphore.h>

#include "utils.cpp"
#include "maw32_utils.cpp"
#include "maw32_trail.cpp"

#include <queue>
#include <vector>
#include <tuple>
#include <thread>
#include <random>
using namespace std;

// Running configuration
typedef struct
{
    float pthresh;
    size_t rounds;
    size_t nthreads;
    uint32_t seed;
} conf_t;

// A representation of a gene for our genetic algorithm
typedef struct
{
    uint8_t diff[16];       // Input difference used
    size_t zero_trails;     // How many output differences of zero observed
    size_t total_trails;    // How many total outputs observed
} gene_t;

// Determine if a gene is 'alive', i.e. usable for breeding
static inline int is_alive(gene_t gene)
{
    return gene.total_trails;
}

// Kill a gene, i.e. make it unusable for breeding
static inline void kill_gene(gene_t *gene)
{
    memset(gene, 0, sizeof(gene_t));
}

// Determine how fit a gene is
static inline double get_fitness(gene_t gene)
{
    // Fitness of a dead gene is zero, otherwise it is the fraction of of zero
    // trails to all observed trails
    if (!is_alive(gene)) return 0.0f;
    else return (gene.zero_trails * 1.0f) / (gene.total_trails * 1.0f);
}

// Pretty-print a gene
static inline void print_gene(gene_t gene, char *annotation)
{
    log(stdout, "(Fingerprint: %02x%02x%02x%02x%02x%02x%02x%02x, Fitness: %lf)%s", 
            gene.diff[0], gene.diff[1], gene.diff[2], gene.diff[3],
            gene.diff[4], gene.diff[5], gene.diff[6], gene.diff[7],
            get_fitness(gene), annotation);
}

// Choose a random gene from the pool, weighted by their fitness
static size_t dice(mt19937& gen, gene_t *pool, size_t len)
{
    // Start by computing the total fitness
    double total_weight = 0.0f;
    for (size_t idx = 0; idx < len; idx++) total_weight += get_fitness(pool[idx]);
    
    // Roll a random number in [0, total_weight]
    double result = (gen() * total_weight) / RAND_MAX;
    
    // Subtract the fitness of each gene from the total weight
    // When we get <=0, select that gene
    for (size_t idx = 0; idx < len; idx++) if (is_alive(pool[idx]))
    {
        result -= get_fitness(pool[idx]);
        if (result <= 0) return idx;
    }
    return -1;
}

// min(a,b)
static inline size_t max(size_t a, size_t b) { return a < b? b : a; }

// Determine if an input differential is viable
static inline bool is_viable(uint8_t *W, const size_t rounds, const float l2pthresh,
                             size_t t, size_t ctr)
{
    const size_t x = 4, y = 1;  // (y/x) expansions must be zero
    if (t >= rounds) return x*ctr >= y*(max(8, rounds) - 8);
    uint8_t w0 = sigma0(W[t-3]),
            w1 = sigma1(W[t-8]);
    for (uint8_t t1: propagate_add(w0, w1, l2pthresh))
    for (uint8_t t2: propagate_add(W[t-4], t1, l2pthresh))
    {
        W[t] = t2;
        if (is_viable(W, rounds, l2pthresh, t+1, ctr + !t2)) return true;
    }
    return false;
}

// Create an input differential randomly.
void make_input_diff(mt19937& gen, uint8_t *sched, size_t rounds, float l2pthresh)
{
    // Zero out the first 4 words
    memset(sched, 0, 4);
    // Randomly assign differences for last four words, and check viability
    do for (int idx = 4; idx < 8; idx++) sched[idx] = gen() & 0xff;
    while (!is_viable(sched, rounds, l2pthresh, 8, 0)) ;
}

// Cross over two message differences
void cross(uint8_t *dest, const uint8_t *left, const uint8_t *right, const size_t mid)
{
    // masks[n] gives the mask for top n bits
    // ~masks[n] gives mask for lower (8-n) bits
    static const uint8_t masks[] = 
    { 0x00, 0x80, 0xc0, 0xe0, 0xf0, 0xf8, 0xfc, 0xfe, 0xff };

    // Byte midpoint
    uint8_t byte_mid = mid / 8;
    // Bit split
    uint8_t bit_split = mid % 8;

    // Copy a chunk of left
    memcpy(dest, left, byte_mid);
    if (bit_split)
    {
        // There is a little bit of left remaining
        dest[byte_mid] = 
            (left[byte_mid]  & masks[bit_split]) | 
            (right[byte_mid] & ~masks[bit_split]);
        memcpy(dest+byte_mid+1, right+byte_mid+1, 8-byte_mid-1);
    }
    else
    {
        // All of left has been copied
        memcpy(dest+byte_mid, right+byte_mid, 8-byte_mid);
    }
}

// Lock for the following queue
pthread_mutex_t queue_lock = PTHREAD_MUTEX_INITIALIZER;
// Queue for storing genes created by slave_make_trails
queue<gene_t> slave_pool;
// Semaphore used to notify get_next_gene when it can run
// The value represents how many values are available in the queue
sem_t queue_notif;

// Push a gene to the slave pool
static inline void put_next_gene(gene_t gene)
{
    pthread_mutex_lock(&queue_lock);
    slave_pool.push(gene);
    sem_post(&queue_notif);
    pthread_mutex_unlock(&queue_lock);
}

// Fetch the next gene from slave_pool
static inline gene_t get_next_gene()
{
    while (1)
    {
        sem_wait(&queue_notif);
        pthread_mutex_lock(&queue_lock);
        gene_t result = slave_pool.front();
        slave_pool.pop();
        pthread_mutex_unlock(&queue_lock);
        return result;
    }
}

// Used as a worker thread to generate potential trails
void *slave_make_trails(void *arg)
{
    conf_t config = *(conf_t *) arg;
    mt19937 gen(config.seed);
    while (1)
    {
        gene_t gene;
        // Build a differential, print it out
        make_input_diff(gen, gene.diff, config.rounds, config.pthresh);
        std::pair<size_t, size_t> result = propagate(gene.diff, config.rounds, config.pthresh);
        if (result.first)
        {
            gene.zero_trails  = result.first;
            gene.total_trails = result.second;
            put_next_gene(gene);
       }
    }
    pthread_exit(NULL);
}

void show_usage()
{
    puts(
        "USAGE: maw_trail [...]\n"
        "\n"
        "Arguments:\n"
        "  -d           Dry run. Runs all setup but does not\n"
        "               generate any trails.\n"
        "  -i           Random only. Does not apply the\n"
        "               genetic algorithm to generate results."
        "  -n count     Specify the number of threads to use.\n"
        "               Defaults to half of the threads on the CPU.\n"
        "  -p prob      Specify the threshhold probability as a\n"
        "               log2 value. Defaults to -3.000000.\n"
        "  -r rounds    Specify the number of rounds to propagate.\n"
        "               Defaults to 8.\n"
        "  -s size      Specify the gene pool size.\n"
        "               Defaults to 32.\n"
        "  -m rate      Specify the immigration rate.\n"
        "               Defaults to 0.05 (5%)");
}

// Entry point
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

    // Default args
    bool dry_run     = false;
    bool random_only = false;
    size_t nthreads  = (size_t) ceil(sysconf(_SC_NPROCESSORS_CONF) / 2.0f);
    float pthresh    = -3.000000f;
    size_t rounds    = 8;
    size_t pool_size = 32;
    float immigration_rate = 0.05;
    char key_fname[] = "/Scratch/key-file-*.******.bin",
         add_fname[] = "/Scratch/add-file-*.******.bin",
         maj_fname[] = "/Scratch/maj-file-*.******.bin";

    // Read args
    int option = -1;
    while ((option = getopt(argc, argv, "hidn:p:r:s:m:")) != -1) switch (option)
    {
        case 'd':
            dry_run = true;
            break;

        case 'i':
            random_only = true;;
            break;

        case 'n':
            nthreads = (unsigned) atoll(optarg);
            ASSERT(nthreads, log(stdout, "Error: Cannot use zero threads"));
            if (nthreads >= sysconf(_SC_NPROCESSORS_CONF))
            {
                log(stdout, "Warning: Requesting to use more threads than available on CPU");
            }
            break;

        case 'p': 
            pthresh = atof(optarg);
            ASSERT(pthresh < 0, log(stdout, "Error: Cannot have a positive threshhhold probability"));
            break;

        case 'r':
            rounds = (unsigned) atoll(optarg);
            ASSERT(rounds >= 1 && rounds <= 16, log(stdout, "Error: Rounds must be between 1 and 16"));
            break;

        case 's':
            pool_size = (unsigned) atoll(optarg);
            ASSERT(pool_size >= 16, log(stdout, "Error: Pool size must be greater than or equal to 16"));
            break;
    
        case 'm':
            immigration_rate = atof(optarg);
            ASSERT(0.0 <= immigration_rate, log(stdout, "Error: Immigration rate must be >= 0"));
            ASSERT(immigration_rate <= 0.5, log(stdout, "Error: Immigration rate must be <= 0.5"));
            break;

        default:
            ASSERT(0);
    }

    // Alter *_fname appropriately, by setting *.****** to pthresh
    sprintf(key_fname+17, "%f.bin", pthresh);
    sprintf(add_fname+17, "%f.bin", pthresh);
    sprintf(maj_fname+17, "%f.bin", pthresh);

    // Setup
    log(stdout, "Initializing...");
    log(stdout, "Threads: %zu", nthreads);
    log(stdout, "Rounds: %zu/16", rounds);
    log(stdout, "Threshold probability: 2^%f", pthresh);
    log(stdout, "Random only: %s", random_only? "true" : "false");
    log(stdout, "Pool size: %zu", pool_size);
    log(stdout, "Immigration rate: %f\n", immigration_rate);
    
    // Load memos
    if (load_key_memo(key_fname)) log(stdout, "Loaded key memos from %s", key_fname);
    else log(stdout, "Failed to load key memos from %s", key_fname);
    if (load_add_memo(add_fname)) log(stdout, "Loaded add memos from %s", add_fname);
    else log(stdout, "Failed to load add memos from %s", add_fname);
    if (load_maj_memo(maj_fname)) log(stdout, "Loaded maj memos from %s", maj_fname);
    else log(stdout, "Failed to load maj memos from %s", maj_fname);
    log(stdout, "Done!\n");
    if (dry_run) return 0;

    // Set up RNGs
    random_device devrand;
    mt19937 gen(devrand());

    // Set up semaphores
    sem_init(&queue_notif, 0, 0);

    // Spawn nthread many workers
    thread *tids = (thread *) calloc(nthreads, sizeof(thread));
    conf_t *configs = (conf_t *) calloc(nthreads, sizeof(conf_t));
    for (size_t idx = 0; idx < nthreads; idx++)
    {
        configs[idx].pthresh  = pthresh;
        configs[idx].rounds   = rounds;
        configs[idx].nthreads = nthreads;
        configs[idx].seed     = devrand();
        tids[idx] = thread(slave_make_trails, configs+idx); 
    }

    // If we're only after random data
    if (random_only) while (1)
    {
        print_gene(get_next_gene(), (char *)" - Immigration");
    }

    // Gather enough differentials to use for genetic algorithms
    // Does not have to be on the stack hence static
    gene_t *pool      = (gene_t *) calloc(pool_size, sizeof(gene_t)),
           *pool_copy = (gene_t *) calloc(pool_size, sizeof(gene_t));
    for (int idx = 0; idx < pool_size; idx++)
    {
        pool[idx] = get_next_gene();
        print_gene(pool[idx], (char *)" - Immigration");
    }
   
    log(stdout, "Beginning optimization");

    // We now have a full gene pool, begin breeding
    for (size_t pool_num = 1; ; pool_num++)
    {
        // Create a copy; we will destroy this to rebuild the real pool
        memset(pool_copy, 0, pool_size * sizeof(gene_t));

        // Allow half the pool to live; this is done by using `dice` to
        // statistically pick the best solutions
        size_t idx = 0;
        for ( ; idx < pool_size/2; )
        {
            // Pick an index of a living gene
            size_t survivor_idx = dice(gen, pool, pool_size);
            // Copy it to the pool
            pool_copy[idx++] = pool[survivor_idx];
            print_gene(pool[survivor_idx], (char *)" - Survivor");
            // and remove the gene from the original pool
            kill_gene(&pool[survivor_idx]);
        }
        // Copy the survivors into the real pool
        memcpy(pool, pool_copy, pool_size * sizeof(gene_t));
        // add immigrants
        for ( ; idx < (int)ceil((pool_size/2) * (1.0+immigration_rate)); idx++)
        {
            pool[idx] = get_next_gene();
            print_gene(pool[idx], (char *)" - Immigration");
        }
        // and begin breeding new genes
        for ( ; idx < pool_size; idx++)
        {
            // Roll a number in [0,16)
            int result = gen() & 0xf;
            // 1/16 to introduce immigrate a gene in
            if (0 && rand == 0) 
            {
                // Blocking operation
                pool[idx] = get_next_gene();
                print_gene(pool[idx], (char *)" - Immigration");
            }
            // Note that we have this loop to prevent getting stuck with 'bad'
            // choices for parent genes
            else while (1)
            {
                // 1/4 to mutate
                if (result < 4)
                {
                    // Pick a random gene from our survivors
                    size_t parent_idx = dice(gen, pool, pool_size/2);
                    // Copy the diff into new_gene
                    memcpy(pool[idx].diff, pool_copy[parent_idx].diff, 8*sizeof(uint8_t));
                    // Flip a random bit in the dense section
                    int bit_idx = 32 + (gen() % 32);
                    pool[idx].diff[bit_idx/8] ^= 1 << (8 - (bit_idx % 8));
                }
                // Rest is crossing over
                else
                {
                    // Pick two random disinct genes
                    size_t parent1_idx = dice(gen, pool, pool_size/2);
                    size_t parent2_idx;
                    do parent2_idx = dice(gen, pool,pool_size/2); while(parent1_idx == parent2_idx);
                    // Pick a midpoint
                    int mid = 32 + (gen() % 32);
                    // and cross
                    cross(pool[idx].diff, 
                            pool_copy[parent1_idx].diff,
                            pool_copy[parent2_idx].diff, 
                            mid);
                }
                // Try to propagate it
                std::pair<size_t, size_t> result = propagate(pool[idx].diff, rounds, pthresh);
                if (result.first)
                {
                    pool[idx].zero_trails  = result.first;
                    pool[idx].total_trails = result.second;
                    print_gene(pool[idx], (char *)" - Generated");
                    break;
                }
                // otherwise continue
            }
        }
        gene_t *best = &pool[0];
        for (int i = 1; i < pool_size; i++)
        {
            if (get_fitness(pool[i]) > get_fitness(*best)) best = &pool[i];
        }
        print_gene(*best, (char *)" - Best");
 
        // Everything has been repopulated.
        log(stdout, "Population %zu bred.", pool_num);
   }
    return 0;
}
#endif // __MAIN
