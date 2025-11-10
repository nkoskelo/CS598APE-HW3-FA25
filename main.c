#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>
#include <stdbool.h>
#include <omp.h>
#include <time.h>

// #define MEASURE_CONVERGENCE

#ifdef MEASURE_CONVERGENCE
char* visited;
int visited_count = 0;
void init_visited(int L) {
    visited = (char*)malloc(sizeof(char) * L * L);
    memset(visited, 0, sizeof(char) * L * L);
    visited_count = 0;
}

void mark_visited(int L, int i, int j) {
    int idx = i * L + j;
    if (!visited[idx]) {
        visited[idx] = 1;
        visited_count++;
    }
}
void free_visited() {
    free(visited);
}
void log_convergence(int L, int step, double energy) {
  double prop_visited = visited_count / (double)(L * L);
  double avg_energy_per_spin = energy / (L * L);
  printf("Convergence for step %d: visited = %d, energy = %f, prop_visited = %f, avg_energy_per_spin = %f\n", step, visited_count, energy, prop_visited, avg_energy_per_spin);
}
#endif

int L;          // Lattice size (L x L)
double T;       // Temperature
double J = 1.0; // Coupling constant
char *lattice;
double *exponential_lookup_table; // Only a small number of options so just precompute them.

const int bits = 3;
const int mask = (1 << bits) - 1;
inline int latticeIndex(int i, int j) {
  
  // int lineIndex = ((i & mask) << bits) + (j & mask);
  // int tileIndex = (i >> bits) * (L >> bits) + (j >> bits);
  // return tileIndex * (1 << (2 * bits)) + lineIndex;
  return i * L + j;
}

// Waste space so that cache coherence doesn't burn us
#define CACHE_LINE_SIZE 64
struct padded_drand48_data
{
  struct drand48_data data;
  char padding[CACHE_LINE_SIZE - sizeof(struct drand48_data)];
};

struct padded_drand48_data *rand_buffers;

void initialize_random_state(int num_threads)
{
  rand_buffers = aligned_alloc(CACHE_LINE_SIZE, num_threads * sizeof(struct padded_drand48_data));
  for (int i = 0; i < num_threads; ++i)
  {
    // Pull seeds from the original RNG in some sufficiently
    // uncorrelated way.
    long seed = rand() + time(NULL) + 100 * i;
    srand48_r(seed, &rand_buffers[i].data);
  }
}

void cleanup_random_state()
{
  free(rand_buffers);
}

void initialize_exponential_table(int num_threads){
  exponential_lookup_table = aligned_alloc(CACHE_LINE_SIZE,(8)*sizeof(double));

  // Only need 0, 2, 4 for indices.
  for (int i = 0; i < 8; i++){
    exponential_lookup_table[i] = 0;
  }
  exponential_lookup_table[0] = 1;
  exponential_lookup_table[2] = exp(-2*J*2 / T);
  exponential_lookup_table[4] = exp(-2*J*4 / T);
}
void free_exponential_table(){
  free(exponential_lookup_table);
}

double gen_rand_double_r()
{
  int tid = omp_get_thread_num();
  double result;
  drand48_r(&rand_buffers[tid].data, &result);
  return result;
}

float tdiff(struct timeval *start, struct timeval *end) {
  return (end->tv_sec - start->tv_sec) + 1e-6 * (end->tv_usec - start->tv_usec);
}

unsigned long long seed = 100;

unsigned long long randomU64() {
  seed ^= (seed << 21);
  seed ^= (seed >> 35);
  seed ^= (seed << 4);
  return seed;
}

double randomDouble() {
  unsigned long long next = randomU64();
  next >>= (64 - 26);
  unsigned long long next2 = randomU64();
  next2 >>= (64 - 26);
  return ((next << 27) + next2) / (double)(1LL << 53);
}



void initializeLattice() {
  lattice = (char *)aligned_alloc(64, sizeof(char) * L * L);
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      lattice[latticeIndex(i, j)] = (randomDouble() < 0.5) ? -1 : 1;
    }
  }
}


double calculateTotalEnergy() {
  double energy = 0.0;

  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      char spin = lattice[latticeIndex(i, j)];

      // Implemented with a wrapping boundary condition.
      char up = lattice[latticeIndex((i - 1 + L) % L, j)];
      char down = lattice[latticeIndex((i + 1) % L, j)];
      char left = lattice[latticeIndex(i, (j - 1 + L) % L)];
      char right = lattice[latticeIndex(i, (j + 1) % L)];

      energy += -J * spin * (up + down + left + right);
    }
  }
  return 0.5 * energy;
}

char calculateEnergyDifferenceIfFlipped(int i, int j){
    // Energy = sum(s_i * s_j) / 2 (since each pair counted twice).
    // Flipping subtracts the contributed term twice, once to 
    // get rid of the current contribution, and once to add the negative contribution.
    char spin = lattice[latticeIndex(i, j)];
    char up = lattice[latticeIndex((i - 1 + L) % L, j)];
    char down = lattice[latticeIndex((i + 1) % L, j)];
    char left = lattice[latticeIndex(i, (j - 1 + L) % L)];
    char right = lattice[latticeIndex(i, (j + 1) % L)];
    return spin * (up + down + left + right);
}

double calculateMagnetization() {
  double mag = 0.0;
  // This part is can be naively parallel
  // but it is not in the timed loop.
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      mag += lattice[latticeIndex(i, j)];
    }
  }
  return mag / (L * L);
}



// row_start : inclusive, row_stop : exclusive
void metropolisHastingsStep(int row_start, int row_stop) {
  // Just needs to be a random value between [0, L)
  int i = (int) (gen_rand_double_r() * (row_stop - row_start)) + row_start;
  int j = (int) (gen_rand_double_r() * L);

  #ifdef MEASURE_CONVERGENCE
    mark_visited(L, i, j);
  #endif

  char dE = calculateEnergyDifferenceIfFlipped(i, j);
  bool accept = dE <= 0.0 || gen_rand_double_r() < exponential_lookup_table[dE];
  if (accept) {
    lattice[latticeIndex(i, j)] *= -1;
  }
}

void saveLatticeImage(const char *png_filename) {
  char ppm_filename[256];
  snprintf(ppm_filename, sizeof(ppm_filename), "temp_%s.ppm", png_filename);

  FILE *f = fopen(ppm_filename, "wb");
  if (!f) {
    printf("Error: Could not create temporary file %s\n", ppm_filename);
    return;
  }

  fprintf(f, "P6\n");
  fprintf(f, "%d %d\n", L, L);
  fprintf(f, "255\n");

  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      unsigned char r, g, b;
      if (lattice[latticeIndex(i, j)] == 1) {
        r = 255;
        g = 255;
        b = 255;
      } else {
        r = 0;
        g = 50;
        b = 200;
      }
      fwrite(&r, 1, 1, f);
      fwrite(&g, 1, 1, f);
      fwrite(&b, 1, 1, f);
    }
  }

  fclose(f);

  char cmd[512];
  snprintf(cmd, sizeof(cmd), "convert %s %s 2>/dev/null", ppm_filename,
           png_filename);
  int result = system(cmd);

  if (result == 0) {
    printf("Saved visualization to %s\n", png_filename);
    remove(ppm_filename);
  } else {
    rename(ppm_filename, png_filename);
    printf("Saved visualization to %s (install ImageMagick for PNG)\n",
           png_filename);
  }
}

void sanityCheck(double energy, double mag_per_spin, const char *stage) {
  double energy_per_spin = energy / (L * L);
  double Tc = 2.0 * J / log(1.0 + sqrt(2.0));

  printf("Sanity check [%s]:\n", stage);

  // 1. Energy per spin
  if (energy_per_spin < -2.0 * J - 0.01 || energy_per_spin > 2.0 * J + 0.01) {
    printf("  [ERROR] Energy per spin (%.4f) outside expected bounds "
           "[%.2f, %.2f]\n",
           energy_per_spin, -2.0 * J, 2.0 * J);
  } else {
    printf("  [OK] Energy per spin = %.4f (within bounds [%.2f, %.2f])\n",
           energy_per_spin, -2.0 * J, 2.0 * J);
  }

  // 2. Magnetization per spin
  if (fabs(mag_per_spin) > 1.01) {
    printf("  [ERROR] Magnetization per spin (%.4f) outside physical bounds "
           "[-1, 1]\n",
           mag_per_spin);
  } else {
    printf("  [OK] Magnetization per spin = %.4f (within bounds [-1, 1])\n",
           mag_per_spin);
  }

  printf("\n");
}

void freeLattice() {
  free(lattice);
}

int main(int argc, const char **argv) {
  if (argc < 4) {
    printf("Usage: %s <lattice_size> <temperature> <steps>\n", argv[0]);
    printf("Example: %s 100 2.269 10000000\n", argv[0]);
    printf("\n2D Ising Model\n");
    printf("Critical temperature: Tc = 2J/ln(1+√2) ≈ 2.26918531421\n");
    return 1;
  }

  L = atoi(argv[1]);
  T = atof(argv[2]);
  int steps = atoi(argv[3]);

  int n_threads = 4;
  int rows_per_thread = (L + n_threads - 1) / n_threads;
  omp_set_num_threads(n_threads);
  initialize_random_state(n_threads);

  initialize_exponential_table(n_threads);

  printf("2D Ising Model\n");
  printf("=================================================\n");
  printf("Lattice size: %d x %d (%d spins)\n", L, L, L * L);
  printf("Temperature: T = %.4f (Tc ≈ 2.269)\n", T);
  printf("Coupling constant: J = %.2f\n", J);
  printf("Number of Metropolis-Hastings steps: %d\n", steps);
  printf("=================================================\n\n");

  initializeLattice();

  double initial_energy = calculateTotalEnergy();
  double initial_mag = calculateMagnetization();
  printf("Initial energy: %.4f\n", initial_energy);
  printf("Initial magnetization: %.4f\n\n", initial_mag);

  sanityCheck(initial_energy, initial_mag, "Initial state");

  saveLatticeImage("initial_state.png");

  #ifdef MEASURE_CONVERGENCE
    init_visited(L);
  #endif

  struct timeval start, end;
  gettimeofday(&start, NULL);


  #pragma omp parallel for schedule(dynamic, 10000)
  // #pragma omp parallel for
  for (int step = 0; step < steps; step++) {
    #ifdef MEASURE_CONVERGENCE
    if (step % 100 == 0) {
      log_convergence(L, step, calculateTotalEnergy());
    }
    #endif
    // One can save off the energy between iterations
    // so that only the update needs to be computed.
    int thread_id = omp_get_thread_num();
    int row_start = thread_id * rows_per_thread;
    int row_stop = (thread_id + 1) * rows_per_thread;
    row_stop = (row_stop > L) ? L : row_stop;
    metropolisHastingsStep(row_start, row_stop);
  }

  gettimeofday(&end, NULL);

  double final_energy = calculateTotalEnergy();
  double final_mag = calculateMagnetization();

  printf("\nFinal energy: %.4f\n", final_energy);
  printf("Final magnetization: %.4f\n\n", final_mag);

  sanityCheck(final_energy, final_mag, "Final state");

  printf("Total time: %0.6f seconds\n\n", tdiff(&start, &end));

  saveLatticeImage("final_state.png");

  freeLattice();
  #ifdef MEASURE_CONVERGENCE
    free_visited();
  #endif
  cleanup_random_state();
  free_exponential_table();
  return 0;
}
