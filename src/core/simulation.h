#pragma once

#include "common.h"
#include <fftw3-mpi.h>
#include <functional>
#include <mpi.h>
#include <omp.h>
#include <random>
#include <vector>

class Simulation {
  MPI_Datatype get_mpi_particle_type() const {
    static MPI_Datatype mpi_particle_type;
    static bool initialized = false;

    if (!initialized) {
      const int n_blocks = 7;
      int block_lengths[n_blocks] = {1, 1, 1, 1, 1, 1, 1};
      MPI_Aint offsets[n_blocks];
      for (int i = 0; i < n_blocks; i++) {
        offsets[i] = sizeof(double) * i;
      }
      MPI_Datatype types[n_blocks] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                                      MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                                      MPI_DOUBLE};
      MPI_Type_create_struct(n_blocks, block_lengths, offsets, types,
                             &mpi_particle_type);
      MPI_Type_commit(&mpi_particle_type);
      initialized = true;
    }

    return mpi_particle_type;
  }

  std::vector<Particle>
  generate_local_particles(int lxstart, int lxres,
                           std::function<double(double, double)> frho,
                           int rank) {
    int num_cells = lxres * params.RESOLUTION;
    std::vector<double> rho_values(num_cells, 0.0);
    double local_density_sum = 0.0;

#pragma omp parallel
    {
      std::mt19937 rng(1234 + rank * 10000 + omp_get_thread_num());
      std::uniform_real_distribution<double> u(0.0, 1.0);

      double thread_sum = 0.0;

#pragma omp for collapse(2) nowait
      for (int i = 0; i < lxres; ++i) {
        for (int j = 0; j < params.RESOLUTION; ++j) {
          int global_i = lxstart + i;
          double x_center = (global_i + 0.5) * dx;
          double y_center = (j + 0.5) * dy;

          double rho = frho(x_center, y_center);
          rho_values[i * params.RESOLUTION + j] = rho;
          thread_sum += rho;
        }
      }

#pragma omp atomic
      local_density_sum += thread_sum;
    }

    double global_density_sum = 0.0;
    MPI_Allreduce(&local_density_sum, &global_density_sum, 1, MPI_DOUBLE,
                  MPI_SUM, comm);
    double mass_per_density = params.MASS / global_density_sum;

    std::vector<std::vector<Particle>> thread_particles(omp_get_max_threads());

#pragma omp parallel
    {
      std::mt19937 rng(4321 + rank * 10000 + omp_get_thread_num());
      std::uniform_real_distribution<double> u(0.0, 1.0);
      int tid = omp_get_thread_num();

#pragma omp for collapse(2) nowait
      for (int i = 0; i < lxres; ++i) {
        for (int j = 0; j < params.RESOLUTION; ++j) {
          int global_i = lxstart + i;
          double x_min = global_i * dx;
          double y_min = j * dy;

          double rho = rho_values[i * params.RESOLUTION + j];
          double cell_mass = rho * mass_per_density;

          double expected_particles = params.PARTICLES_PER_CELL;
          std::poisson_distribution<int> poisson(expected_particles);

          int n_particles = poisson(rng);

          if (n_particles == 0)
            continue;

          double particle_mass = cell_mass / n_particles;

          for (int p = 0; p < n_particles; ++p) {
            double x = x_min + dx * u(rng);
            double y = y_min + dy * u(rng);

            thread_particles[tid].emplace_back(
                Particle{vec2(x - params.RADIUS / 2, y - params.RADIUS / 2),
                         vec2(0.0), vec2(0.0), particle_mass});
          }
        }
      }
    }

    std::vector<Particle> particles;
    for (auto &vec : thread_particles)
      particles.insert(particles.end(), vec.begin(), vec.end());

    return particles;
  }

public:
  struct Params {
    float GRAVITY = 0.01;         // Gravitational constant
    float SOFTENING = 0.001;      // Softening length
    float TIMESTEP = 0.01;        // Integration timestep
    float RADIUS = 10.0;          // Periodic boundary condition radius
    float MASS = 50.0;            // Total mass of the universe
    float PARTICLES_PER_CELL = 4; // Total number of particles
    ptrdiff_t RESOLUTION = 1024;  // Resolution for density texture
  };

  Params params;

  std::vector<Particle> particles;
  std::vector<double> density;

  // MPI variables
  int rank;
  int size;

  float dx;
  float dy;

  ptrdiff_t lxres, lxstart = 0;
  ptrdiff_t lalloc = 0;

  double *lrho = nullptr;
  double *lphi = nullptr;
  fftw_complex *rhok = nullptr;

  std::vector<double> lrho_ext;

  fftw_plan fplan = nullptr;
  fftw_plan bplan = nullptr;

  MPI_Comm comm = MPI_COMM_NULL;

  Simulation(Params params) : params(params) {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    MPI_Comm_split(MPI_COMM_WORLD, world_rank != 0, world_rank, &comm);
    if (world_rank == 0) {
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      size -= 1;
      return;
    }

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    lalloc = fftw_mpi_local_size_2d(
        params.RESOLUTION, params.RESOLUTION / 2 + 1, comm, &lxres, &lxstart);
    lrho = fftw_alloc_real(2 * lalloc);
    lrho_ext = std::vector<double>((lxres + 2) * params.RESOLUTION, 0.0);
    lphi = fftw_alloc_real(2 * lalloc);
    rhok = fftw_alloc_complex(lalloc);
    fplan = fftw_mpi_plan_dft_r2c_2d(params.RESOLUTION, params.RESOLUTION, lrho,
                                     rhok, comm, FFTW_MEASURE);
    bplan = fftw_mpi_plan_dft_c2r_2d(params.RESOLUTION, params.RESOLUTION, rhok,
                                     lphi, comm, FFTW_MEASURE);

    dx = params.RADIUS / params.RESOLUTION;
    dy = params.RADIUS / params.RESOLUTION;

    // Initialize particles
    particles = generate_local_particles(
        lxstart, lxres, [&](double x, double y) { return 1.0; }, rank);

    // Initialize density texture
    mass_assignment();
  }

  ~Simulation() {
    if (rank != 0) {
      fftw_destroy_plan(fplan);
      fftw_destroy_plan(bplan);
      fftw_free(lrho);
      fftw_free(lphi);
      fftw_free(rhok);

      if (comm != MPI_COMM_NULL) {
        MPI_Comm_free(&comm);
      }
    }
  }

  void timestep() {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == 0) {
      return;
    }
    mass_assignment();
  }

  void mass_assignment() {
    std::fill(lrho_ext.begin(), lrho_ext.end(), 0.0);

    std::vector<double> send_left(params.RESOLUTION, 0.0);
    std::vector<double> send_right(params.RESOLUTION, 0.0);

    for (const auto &p : particles) {
      double fx = std::fmod((p.p.x + params.RADIUS * 0.5) / dx,
                            static_cast<double>(params.RESOLUTION));
      double fy = std::fmod((p.p.y + params.RADIUS * 0.5) / dy,
                            static_cast<double>(params.RESOLUTION));
      if (fx < 0)
        fx += params.RESOLUTION;
      if (fy < 0)
        fy += params.RESOLUTION;

      int gx = static_cast<int>(std::floor(fx));
      int gy = static_cast<int>(std::floor(fy));

      double dx1 = fx - gx, dx0 = 1.0 - dx1;
      double dy1 = fy - gy, dy0 = 1.0 - dy1;

      for (int di = 0; di <= 1; ++di) {
        int i_glob = (gx + di) % params.RESOLUTION;
        double wx = (di == 0 ? dx0 : dx1);

        bool owned = (i_glob >= lxstart) && (i_glob < lxstart + lxres);

        for (int dj = 0; dj <= 1; ++dj) {
          int j_glob = (gy + dj) % params.RESOLUTION;
          double wy = (dj == 0 ? dy0 : dy1);

          double mass = (wx * wy) * (p.mass / (dx * dy));

          if (owned) {
            int lx = i_glob - lxstart + 1;
            lrho_ext[lx * params.RESOLUTION + j_glob] += mass;
          } else {
            int left_edge =
                (lxstart - 1 + params.RESOLUTION) % params.RESOLUTION;
            int right_edge = (lxstart + lxres) % params.RESOLUTION;

            if (i_glob == left_edge)
              send_left[j_glob] += mass;
            else if (i_glob == right_edge)
              send_right[j_glob] += mass;
          }
        }
      }
    }

    std::vector<double> recv_left(params.RESOLUTION, 0.0);
    std::vector<double> recv_right(params.RESOLUTION, 0.0);

    int left = (rank == 0 ? size - 1 : rank - 1);
    int right = (rank == size - 1 ? 0 : rank + 1);

    MPI_Sendrecv(send_left.data(), params.RESOLUTION, MPI_DOUBLE, left, 0,
                 recv_right.data(), params.RESOLUTION, MPI_DOUBLE, right, 0,
                 comm, MPI_STATUS_IGNORE);

    MPI_Sendrecv(send_right.data(), params.RESOLUTION, MPI_DOUBLE, right, 1,
                 recv_left.data(), params.RESOLUTION, MPI_DOUBLE, left, 1, comm,
                 MPI_STATUS_IGNORE);

    for (int j = 0; j < params.RESOLUTION; ++j) {
      lrho_ext[1 * params.RESOLUTION + j] += recv_left[j];
      lrho_ext[(lxres)*params.RESOLUTION + j] += recv_right[j];
    }

    MPI_Sendrecv(&lrho_ext[lxres * params.RESOLUTION], params.RESOLUTION,
                 MPI_DOUBLE, right, 2, &lrho_ext[0], params.RESOLUTION,
                 MPI_DOUBLE, left, 2, comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&lrho_ext[1 * params.RESOLUTION], params.RESOLUTION,
                 MPI_DOUBLE, left, 3,
                 &lrho_ext[(lxres + 1) * params.RESOLUTION], params.RESOLUTION,
                 MPI_DOUBLE, right, 3, comm, MPI_STATUS_IGNORE);

    for (int i = 0; i < lxres; ++i)
      std::copy_n(&lrho_ext[(i + 1) * params.RESOLUTION], params.RESOLUTION,
                  &lrho[i * params.RESOLUTION]);

    MPI_Barrier(comm);
  }

  // Gather methods

  std::vector<Particle> gather_particles() const {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == 0) {
      int total_particles = 0;
      MPI_Recv(&total_particles, 1, MPI_INT, 1, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);

      std::vector<Particle> particles_global(total_particles);
      MPI_Recv(particles_global.data(), total_particles * sizeof(Particle),
               MPI_BYTE, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      return particles_global;
    }

    if (comm == MPI_COMM_NULL)
      return {};

    int local_count = static_cast<int>(particles.size());

    std::vector<int> counts(size), displs(size);
    MPI_Gather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);

    std::vector<Particle> global_particles;
    if (rank == 0) {
      int total = 0;
      displs[0] = 0;
      for (int i = 0; i < size; ++i) {
        if (i > 0)
          displs[i] = displs[i - 1] + counts[i - 1];
        total += counts[i];
      }
      global_particles.resize(total);
    }

    MPI_Datatype mpi_particle_type = get_mpi_particle_type();
    MPI_Gatherv(particles.data(), local_count, mpi_particle_type,
                global_particles.data(), counts.data(), displs.data(),
                mpi_particle_type, 0, comm);

    if (rank == 0) {
      int total = static_cast<int>(global_particles.size());
      MPI_Send(&total, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(global_particles.data(), total * sizeof(Particle), MPI_BYTE, 0,
               1, MPI_COMM_WORLD);
    }

    return {};
  }

  std::vector<double> gather_rho() const {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == 0) {
      std::vector<double> rho_global(params.RESOLUTION * params.RESOLUTION);
      MPI_Recv(rho_global.data(), rho_global.size(), MPI_DOUBLE, 1, 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      return rho_global;
    }

    if (comm == MPI_COMM_NULL)
      return {};

    ptrdiff_t local_size = lxres * params.RESOLUTION;
    int local_count = static_cast<int>(local_size);

    std::vector<int> recvcounts(size), displs(size);
    MPI_Gather(&local_count, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0,
               comm);

    if (rank == 0) {
      displs[0] = 0;
      for (int i = 1; i < size; ++i)
        displs[i] = displs[i - 1] + recvcounts[i - 1];
    }

    std::vector<double> rho_global;
    if (rank == 0)
      rho_global.resize(params.RESOLUTION * params.RESOLUTION);

    MPI_Gatherv(lrho, local_count, MPI_DOUBLE, rho_global.data(),
                recvcounts.data(), displs.data(), MPI_DOUBLE, 0, comm);

    if (rank == 0) {
      MPI_Send(rho_global.data(), rho_global.size(), MPI_DOUBLE, 0, 0,
               MPI_COMM_WORLD);
    }

    return {};
  }
};
