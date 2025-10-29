#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>
#include <stdbool.h>

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

int L;          // Lattice size (L x L)
double T;       // Temperature
double J = 1.0; // Coupling constant
int **lattice;

void initializeLattice() {
  lattice = (int **)malloc(sizeof(int *) * L);
  for (int i = 0; i < L; i++) {
    lattice[i] = (int *)malloc(sizeof(int) * L);
    for (int j = 0; j < L; j++) {
      lattice[i][j] = (randomDouble() < 0.5) ? -1 : 1;
    }
  }
}

double calculateTotalEnergy() {
  double energy = 0.0;

  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      int spin = lattice[i][j];

      // Implemented with a wrapping boundary condition.
      int up = lattice[(i - 1 + L) % L][j];
      int down = lattice[(i + 1) % L][j];
      int left = lattice[i][(j - 1 + L) % L];
      int right = lattice[i][(j + 1) % L];

      energy += -J * spin * (up + down + left + right);
    }
  }
  return 0.5 * energy;
}

double calculateEnergyDifferenceIfFlipped(int i, int j){
    double delta = 0.0;
    // Energy = sum(s_i * s_j) / 2 (since each pair counted twice).
    // Flipping subtracts the contributed term twice, once to 
    // get rid of the current contribution, and once to add the negative contribution.
    int spin = lattice[i][j];
    int up = lattice[(i - 1 + L) % L][j];
    int down = lattice[(i + 1) % L][j];
    int left = lattice[i][(j - 1 + L) % L];
    int right = lattice[i][(j + 1) % L];
    delta +=  2 * J * spin * (up + down + left + right);
    return delta;
}

double calculateMagnetization() {
  double mag = 0.0;
  // This part is can be naively parallel
  // but it is not in the timed loop.
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      mag += lattice[i][j];
    }
  }
  return mag / (L * L);
}

void metropolisHastingsStep() {
  // Just needs to be a random value between [0, L)
  int i = (int)(randomDouble() * L);
  int j = (int)(randomDouble() * L);

  double dE = calculateEnergyDifferenceIfFlipped(i, j);
  bool accept = dE <= 0.0 || randomDouble() < exp(-dE / T);
  if (accept) {
    lattice[i][j] *= -1;
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
      if (lattice[i][j] == 1) {
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
  for (int i = 0; i < L; i++) {
    free(lattice[i]);
  }
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

  struct timeval start, end;
  gettimeofday(&start, NULL);

  for (int step = 0; step < steps; step++) {
    // One can save off the energy between iterations
    // so that only the update needs to be computed.
    metropolisHastingsStep();
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
  return 0;
}
