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
    x -= 0.5;
    y -= 0.5;
    constexpr double sigma = 0.05;
    constexpr double norm = 1.0 / (2.0 * M_PI * sigma * sigma);
    return norm * std::exp(-(x * x + y * y) / (2.0 * sigma * sigma));
  };

  struct Params {
    float GRAVITY = 1.0;           // Gravitational constant
    float SOFTENING = 0.2;         // Softening length
    float TIMESTEP = 0.014;        // Integration timestep
    float RADIUS = 1.0;            // Periodic boundary condition radius
    float MASS = 1.0;              // Total mass of the universe
    float PARTICLES_PER_CELL = 8;  // Total number of particles
    bool USE_SCALE_FACTOR = false;
    ptrdiff_t RESOLUTION = 1024;  // Resolution for density texture
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

  fftw_complex* scratch_k = nullptr;  // Fourier scratch space
  fftw_complex* xgrad = nullptr;      // gradient of potential (x-axis)
  fftw_complex* ygrad = nullptr;      // gradient of potential (y-axis)
  fftw_plan fplan = nullptr;
  fftw_plan bxplan = nullptr;
  fftw_plan byplan = nullptr;

  std::vector<vec2> ff;  // force field

  double a = 1.0;     // scale factor
  double H = 0.0;     // scale factor
  double adot = 0.0;  // scale factor time derivative
  double rho0 = 0.0;  // mean density at a=1

  double t = 0.0;
  double dt = 0.0;

  std::vector<Particle> particles;

 public:
  Simulation(Params params) : params(params) {
    Nx = params.RESOLUTION;
    Ny = params.RESOLUTION;
    Lx = params.RADIUS;
    Ly = params.RADIUS;
    dx = (double)Lx / Nx;
    dy = (double)Ly / Ny;
    dt = params.TIMESTEP;

    rho0 = params.MASS / (Lx * Ly);
    adot = std::sqrt(8.0 * M_PI * params.GRAVITY * rho0);
    H = adot / a;

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);

    MPI_Comm_split(MPI_COMM_WORLD, wrank != 0, wrank, &comm);
    if (wrank == 0)
      return;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    ptrdiff_t lalloc = fftw_mpi_local_size_2d(Nx, Ny, comm, &lNx, &lx0);
    scratch_k = fftw_alloc_complex(lalloc);
    xgrad = fftw_alloc_complex(lalloc);
    ygrad = fftw_alloc_complex(lalloc);
    fplan = fftw_mpi_plan_dft_2d(Nx, Ny, scratch_k, scratch_k, comm, FFTW_FORWARD, FFTW_MEASURE);
    bxplan = fftw_mpi_plan_dft_2d(Nx, Ny, xgrad, xgrad, comm, FFTW_BACKWARD, FFTW_MEASURE);
    byplan = fftw_mpi_plan_dft_2d(Nx, Ny, ygrad, ygrad, comm, FFTW_BACKWARD, FFTW_MEASURE);

    rho = std::vector<double>(lNx * Ny, 0.0);
    rho_ext = std::vector<double>((lNx + 2) * Ny, 0.0);
    ff = std::vector<vec2>(lNx * Ny, vec2(0.0));

    generate_particles(seed_uniform);
    assign_masses();
    compute_forces();
  }

  void timestep() {
    if (wrank != 0) {
      update_positions();
      assign_masses();
      compute_forces();
    }

    t += dt;
  }

  ~Simulation() {
    fftw_destroy_plan(fplan);
    fftw_destroy_plan(bxplan);
    fftw_destroy_plan(byplan);
    fftw_free(scratch_k);
    fftw_free(xgrad);
    fftw_free(ygrad);

    if (comm != MPI_COMM_NULL)
      MPI_Comm_free(&comm);
  }

  std::vector<Particle> gather_particles() const;
  std::vector<double> gather_rho() const;
  std::vector<vec2> gather_ff() const;

 private:
  void generate_particles(seed_density seed);
  void assign_masses();
  void compute_forces();
  void update_positions();

  vec2 cic_force(vec2 p);
};
