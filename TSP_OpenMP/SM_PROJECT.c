#define _CRT_SECURE_NO_WARNINGS 
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <time.h> 
#include <omp.h> 

#define N 200
#define T_START 1000.0
#define T_MIN 0.001
#define ALPHA 0.999
#define ITERATIONS 5000000

// Data Type
typedef struct {
    double x, y;
} City;

City cities[N];
int global_best_tour[N];   // Stores the globally best tour among all threads

//Functions

// Calculate Euclidean distance between two cities
double distance_between(City a, City b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return sqrt(dx * dx + dy * dy);
}

// Compute total tour length (closed loop)
double total_distance(const int tour[]) {
    double s = 0.0;
    for (int i = 0; i < N - 1; ++i)
        s += distance_between(cities[tour[i]], cities[tour[i + 1]]);
    s += distance_between(cities[tour[N - 1]], cities[tour[0]]);
    return s;
}

// Swap two city indices
void swap_idx(int tour[], int i, int j) {
    int t = tour[i];
    tour[i] = tour[j];
    tour[j] = t;
}

// Linear Congruential Generator (LCG) , it is a simple per-thread random generator
static inline unsigned int lcg_next(unsigned int* state) {
    *state = (*state * 1103515245u + 12345u) & 0x7fffffffu;
    return *state;
}

// Return random number between 0 and 1
static inline double rand01(unsigned int* state) {
    return (double)lcg_next(state) / (double)0x7fffffffu;
}

// Shuffle tour using the thread's random seed
void shuffle_tour(int tour[], unsigned int* seed) {
    for (int i = 0; i < N; ++i)
        tour[i] = i;
    for (int i = N - 1; i > 0; --i) {
        int j = (int)(lcg_next(seed) % (unsigned int)(i + 1));
        swap_idx(tour, i, j);
    }
}

//Simulated Annealing

// Single run of Simulated Annealing (no change made here)
double sa_single_run(unsigned int seed_init, int best_tour_out[]) {
    unsigned int seed = seed_init;
    int tour[N];
    int best_tour[N];

    shuffle_tour(tour, &seed);

    double current_cost = total_distance(tour);
    double best_cost = current_cost;
    double T = T_START;

    for (long long step = 0; step < ITERATIONS && T > T_MIN; ++step) {
        int i = (int)(lcg_next(&seed) % N);
        int j = (int)(lcg_next(&seed) % N);
        while (j == i)
            j = (int)(lcg_next(&seed) % N);

        swap_idx(tour, i, j);
        double new_cost = total_distance(tour);
        double delta = new_cost - current_cost;

        // Accept new solution if better, or probabilistically if worse
        if (delta < 0 || rand01(&seed) < exp(-delta / T)) {
            current_cost = new_cost;
        }
        else {
            swap_idx(tour, i, j); // Undo the swap
        }

        // Update best solution found so far
        if (current_cost < best_cost) {
            best_cost = current_cost;
            for (int k = 0; k < N; ++k)
                best_tour[k] = tour[k];
        }

        // Cooling
        T *= ALPHA;
    }

    // Copy best tour to output
    for (int k = 0; k < N; ++k)
        best_tour_out[k] = best_tour[k];

    return best_cost;
}

// Parallel Multi-Start SA

// Run multiple SA in parallel threads and find the global best
double run_parallel_SA(int num_threads, double* best_out, int* best_thread_out) {
    double global_best_distance = 1e300;
    int global_best_thread = -1;
    unsigned int base_seed = (unsigned int)time(NULL);

    omp_set_num_threads(num_threads);
    double start = omp_get_wtime();

#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        unsigned int seed = (unsigned int)(clock() ^ (time(NULL) + tid * 2654435761u));

        int local_best_tour[N];
        double local_best_distance = sa_single_run(seed, local_best_tour);

        // Update global best safely (critical section)
#pragma omp critical
        {
            if (local_best_distance < global_best_distance) {
                global_best_distance = local_best_distance;
                global_best_thread = tid;
                for (int k = 0; k < N; ++k)
                    global_best_tour[k] = local_best_tour[k];
            }
        }
    }

    double end = omp_get_wtime();
    double elapsed_ms = (end - start) * 1000.0;

    *best_out = global_best_distance;
    *best_thread_out = global_best_thread;

    return elapsed_ms;
}

//Main

int main() {
    // 1) Read city coordinates
    FILE* f = fopen("C:\\Users\\HP\\Downloads\\cities200.txt", "r");
    if (!f) {
        printf("Could not open the file!\n");
        return 1;
    }

    char header[256];
    fgets(header, sizeof(header), f);

    int id;
    for (int i = 0; i < N; ++i)
        fscanf(f, "%d,%lf,%lf", &id, &cities[i].x, &cities[i].y);
    fclose(f);

    // 2) Define thread counts to test
    int thread_counts[] = { 1, 4, 8, 16 };
    int num_tests = sizeof(thread_counts) / sizeof(thread_counts[0]);

    double times[16], speedup[16], efficiency[16], bests[16];
    int best_threads[16];

    double seq_time_ms = 0.0;
    double best_dist;
    int best_thread;

    // For tracking overall best result
    double global_best_dist = 1e300;
    int global_best_idx = -1;

    // 3) Run tests for each number of threads
    for (int i = 0; i < num_tests; ++i) {
        int threads = thread_counts[i];
        double t_ms = run_parallel_SA(threads, &best_dist, &best_thread);

        if (threads == 1)
            seq_time_ms = t_ms;  // baseline for speedup

        double s = seq_time_ms / t_ms;
        double e = s / threads;

        times[i] = t_ms;
        speedup[i] = s;
        efficiency[i] = e;
        bests[i] = best_dist;
        best_threads[i] = best_thread;

        // Keep track of the best global result
        if (best_dist < global_best_dist) {
            global_best_dist = best_dist;
            global_best_idx = i;
        }
    }

    // 4) Print results table
    printf("\nThreads    Time(ms)    Speedup    Eff(%%)    BestDist     BestThread\n");
    for (int i = 0; i < num_tests; ++i) {
        printf("%-10d %-10.2f %-10.2f %-9.2f %-12.4f %-10d\n",
            thread_counts[i],
            times[i],
            speedup[i],
            efficiency[i] * 100.0,
            bests[i],
            best_threads[i]);
    }

    //Print the overall best result 
    printf("\n Overall Best Result \n");
    printf("Threads: %d | Best Distance: %.4f\n",
        thread_counts[global_best_idx],
        global_best_dist);

    //Print the best tour path 
    printf("Best Tour Order:\n");
    for (int i = 0; i < N; ++i)
        printf("%d ", global_best_tour[i]);
    printf("\n");
    return 0;
}
