#include "simulation_v2.h"
#include <mpi.h>
#include <omp.h>
#include <cstddef>
#include <random>

std::mt19937 thread_rng() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::mt19937 rng(1234 + rank * 10000 + omp_get_thread_num());
  return rng;
}

void Simulation::generate_particles(seed_density seed, int particles_per_cell) {
  int num_cells = lNx * Ny;
  std::vector<double> rho_values(num_cells, 0.0);

  double sum = 0.0;
#pragma omp parallel for reduction(+ : sum) collapse(2)
  for (int i = 0; i < lNx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      int gi = lx0 + i;
      rho_values[i * lNx + j] = seed(gi, j);
      sum += rho_values[i * lNx + j];
    }
  }

  double global_sum = 0.0;
  MPI_Allreduce(&global_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
  double mass_per_density = params.MASS / global_sum;

  std::vector<std::vector<Particle>> thread_particles(omp_get_max_threads());

#pragma omp parallel
  {
    std::mt19937 rng(4321 + rank * 10000 + omp_get_thread_num());
    std::uniform_real_distribution<double> u(0.0, 1.0);
    int tid = omp_get_thread_num();

#pragma omp for collapse(2) nowait
    for (int i = 0; i < lNx; ++i) {
      for (int j = 0; j < Ny; ++j) {
        int gi = lx0 + i;
        double x_min = gi * dx;
        double y_min = j * dy;

        double rho = rho_values[i * lNx + j];
        double cell_mass = rho * mass_per_density;

        std::poisson_distribution<int> poisson(particles_per_cell);
        int n_particles = poisson(rng);

        if (n_particles == 0)
          continue;

        double particle_mass = cell_mass / n_particles;
        for (int p = 0; p < n_particles; ++p) {
          double x = x_min + dx * u(rng);
          double y = y_min + dy * u(rng);

          thread_particles[tid].emplace_back(Particle{vec2(x, y), vec2(0.0), vec2(0.0), particle_mass});
        }
      }
    }
  }

  for (auto& vec : thread_particles)
    particles.insert(particles.end(), vec.begin(), vec.end());
}
