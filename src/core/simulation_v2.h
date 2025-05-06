#pragma once

#include <fftw3-mpi.h>
#include <mpi.h>
#include <cstddef>
#include <functional>
#include <vector>
#include "common.h"

class Simulation {
 public:
  using seed_density = std::function<double(double, double)>;
  inline static const seed_density seed_uniform = [](double x, double y) {
    return 1.0;
  };
  inline static const seed_density seed_gaussian = [](double x, double y) {
    constexpr double sigma = 0.5;
    constexpr double norm = 1.0 / (2.0 * M_PI * sigma * sigma);
    return norm * std::exp(-(x * x + y * y) / (2.0 * sigma * sigma));
  };

  struct Params {
    float GRAVITY = 0.01;          // Gravitational constant
    float SOFTENING = 0.001;       // Softening length
    float TIMESTEP = 0.01;         // Integration timestep
    float RADIUS = 10.0;           // Periodic boundary condition radius
    float MASS = 50.0;             // Total mass of the universe
    float PARTICLES_PER_CELL = 4;  // Total number of particles
    ptrdiff_t RESOLUTION = 1024;   // Resolution for density texture
  };

  Params params;

  MPI_Comm comm = MPI_COMM_NULL;  // worker MPI_Comm
  int rank, size = -1;            // worker rank/size
  int wrank, wsize = -1;          // global rank/size (including main process)

  ptrdiff_t Nx, Ny;    // mesh dimensions (cells)
  double Lx, Ly;       // mesh dimensions (length)
  double dx, dy;       // cell dimensions (length)
  ptrdiff_t lNx, lx0;  // rank strip length, strip start

  std::vector<double> rho;      // density field
  std::vector<double> rho_ext;  // density field (with halo)

  fftw_complex* rho_k = nullptr;  // complex density field (frequency space)
  fftw_plan rho_fft = nullptr;

  std::vector<Particle> particles;

 public:
  Simulation(Params params) : params(params) {
    Nx = params.RESOLUTION;
    Ny = params.RESOLUTION;
    Lx = params.RADIUS;
    Ly = params.RADIUS;
    dx = (double)Lx / Nx;
    dy = (double)Ly / Ny;

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);

    MPI_Comm_split(MPI_COMM_WORLD, wrank != 0, wrank, &comm);
    if (wrank == 0)
      return;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    ptrdiff_t lalloc = fftw_mpi_local_size_2d(Nx, Ny, comm, &lNx, &lx0);
    rho_k = fftw_alloc_complex(lalloc);
    rho_fft = fftw_plan_dft_2d(Nx, Ny, rho_k, rho_k, FFTW_FORWARD, FFTW_MEASURE);

    rho = std::vector<double>(lNx * Ny, 0.0);
    rho_ext = std::vector<double>((lNx + 2) * Ny, 0.0);

    generate_particles(seed_uniform);
    assign_masses();
  }

  ~Simulation() {
    fftw_destroy_plan(rho_fft);
    fftw_free(rho_k);

    if (comm != MPI_COMM_NULL)
      MPI_Comm_free(&comm);
  }

  std::vector<Particle> gather_particles() const;
  std::vector<double> gather_rho() const;

 private:
  void generate_particles(seed_density seed);
  void assign_masses();
};
