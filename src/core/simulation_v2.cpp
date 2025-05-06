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

void Simulation::generate_particles(seed_density seed) {
  particles.clear();

  int num_cells = lNx * Ny;
  std::vector<double> rho_values(num_cells, 0.0);

  double local_sum = 0.0;
#pragma omp parallel for reduction(+ : local_sum) collapse(2)
  for (int i = 0; i < lNx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      int gi = lx0 + i;
      double xc = (gi + 0.5) * dx;
      double yc = (j + 0.5) * dy;
      double rho = seed(xc, yc);

      rho_values[i * Ny + j] = rho;
      local_sum += rho;
    }
  }

  double global_sum = 0.0;
  MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
  double mass_per_density = params.MASS / global_sum;

  std::vector<std::vector<Particle>> thread_particles(omp_get_max_threads());

#pragma omp parallel
  {
    std::mt19937 rng = thread_rng();
    std::uniform_real_distribution<double> u(0.0, 1.0);
    int tid = omp_get_thread_num();

#pragma omp for collapse(2) nowait
    for (int i = 0; i < lNx; ++i) {
      for (int j = 0; j < Ny; ++j) {
        int gi = lx0 + i;
        double x_min = gi * dx;
        double y_min = j * dy;

        double rho = rho_values[i * Ny + j];
        double cell_mass = rho * mass_per_density;

        std::poisson_distribution<int> poisson(params.PARTICLES_PER_CELL);
        int n_particles = poisson(rng);

        if (n_particles == 0)
          continue;

        double particle_mass = cell_mass / n_particles;
        for (int p = 0; p < n_particles; ++p) {
          double x = x_min + dx * u(rng);
          double y = y_min + dy * u(rng);

          thread_particles[tid].emplace_back(
              Particle{vec2(x, y), vec2(0.0), vec2(0.0), particle_mass});
        }
      }
    }
  }

  for (auto& vec : thread_particles)
    particles.insert(particles.end(), vec.begin(), vec.end());
}

void Simulation::assign_masses() {
  std::fill(rho_ext.begin(), rho_ext.end(), 0.0);

  std::vector<double> send_left(Ny, 0.0), send_right(Ny, 0.0);

  // precompute once per particle
  int left_edge = (int(lx0) - 1 + int(Nx)) % int(Nx);
  int right_edge = (int(lx0) + int(lNx)) % int(Nx);

  for (auto& p : particles) {
    double fx = std::fmod(p.p.x / dx, Nx);
    double fy = std::fmod(p.p.y / dy, Ny);
    if (fx < 0)
      fx += Nx;
    if (fy < 0)
      fy += Ny;

    int gx = int(std::floor(fx));
    int gy = int(std::floor(fy));

    double dx1 = fx - gx, dx0 = 1 - dx1;
    double dy1 = fy - gy, dy0 = 1 - dy1;

    for (int di = 0; di <= 1; ++di) {
      int i_glob = (gx + di + Nx) % Nx;
      double wx = (di == 0 ? dx0 : dx1);

      bool owned = (i_glob >= lx0 && i_glob < lx0 + lNx);

      for (int dj = 0; dj <= 1; ++dj) {
        int j_glob = (gy + dj + Ny) % Ny;
        double wy = (dj == 0 ? dy0 : dy1);

        double m = p.mass * wx * wy / (dx * dy);

        if (owned) {
          int iext = (i_glob - lx0) + 1;
          rho_ext[size_t(iext) * size_t(Ny) + size_t(j_glob)] += m;
        } else {
          if (i_glob == left_edge)
            send_left[j_glob] += m;
          if (i_glob == right_edge)
            send_right[j_glob] += m;
        }
      }
    }
  }

  std::vector<double> recv_left(Ny, 0.0);
  std::vector<double> recv_right(Ny, 0.0);

  int left = (rank == 0 ? size - 1 : rank - 1);
  int right = (rank == size - 1 ? 0 : rank + 1);

  MPI_Sendrecv(send_left.data(), Ny, MPI_DOUBLE, left, 0, recv_right.data(), Ny, MPI_DOUBLE, right,
               0, comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(send_right.data(), Ny, MPI_DOUBLE, right, 1, recv_left.data(), Ny, MPI_DOUBLE, left,
               1, comm, MPI_STATUS_IGNORE);

  for (int j = 0; j < params.RESOLUTION; ++j) {
    rho_ext[1 * params.RESOLUTION + j] += recv_left[j];
    rho_ext[(lNx)*params.RESOLUTION + j] += recv_right[j];
  }

  MPI_Sendrecv(&rho_ext[lNx * Ny], Ny, MPI_DOUBLE, right, 2, &rho_ext[0], Ny, MPI_DOUBLE, left, 2,
               comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(&rho_ext[1 * Ny], Ny, MPI_DOUBLE, left, 3, &rho_ext[(lNx + 1) * Ny], Ny, MPI_DOUBLE,
               right, 3, comm, MPI_STATUS_IGNORE);

  for (int i = 0; i < lNx; ++i)
    std::copy_n(&rho_ext[(i + 1) * Ny], Ny, &rho[i * Ny]);
}

const MPI_Datatype get_mpi_particle_type() {
  static MPI_Datatype mpi_particle_type;
  static bool initialized = false;

  if (!initialized) {
    const int n_blocks = 7;
    int block_lengths[n_blocks] = {1, 1, 1, 1, 1, 1, 1};
    MPI_Aint offsets[n_blocks];
    for (int i = 0; i < n_blocks; i++) {
      offsets[i] = sizeof(double) * i;
    }
    MPI_Datatype types[n_blocks] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                                    MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Type_create_struct(n_blocks, block_lengths, offsets, types, &mpi_particle_type);
    MPI_Type_commit(&mpi_particle_type);
    initialized = true;
  }

  return mpi_particle_type;
}

std::vector<Particle> Simulation::gather_particles() const {
  if (wrank == 0) {
    int N = 0;
    MPI_Recv(&N, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    std::vector<Particle> particles_global(N);
    MPI_Recv(particles_global.data(), N * sizeof(Particle), MPI_BYTE, 1, 1, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    return particles_global;
  }

  int N = particles.size();
  std::vector<int> counts(size), displs(size);
  MPI_Gather(&N, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);

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
  MPI_Gatherv(particles.data(), N, mpi_particle_type, global_particles.data(), counts.data(),
              displs.data(), mpi_particle_type, 0, comm);

  if (rank == 0) {
    int total = static_cast<int>(global_particles.size());
    MPI_Send(&total, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(global_particles.data(), total * sizeof(Particle), MPI_BYTE, 0, 1, MPI_COMM_WORLD);
  }

  return {};
}

std::vector<double> Simulation::gather_rho() const {
  if (wrank == 0) {
    std::vector<double> rho_global(Nx * Ny);
    MPI_Recv(rho_global.data(), rho_global.size(), MPI_DOUBLE, 1, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    return rho_global;
  }

  std::vector<int> recvcounts(size), displs(size);
  int local_count = lNx * Ny;
  MPI_Gather(&local_count, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, comm);

  if (rank == 0) {
    displs[0] = 0;
    for (int i = 1; i < size; ++i)
      displs[i] = displs[i - 1] + recvcounts[i - 1];
  }

  std::vector<double> rho_global;
  if (rank == 0)
    rho_global.resize(Nx * Ny);
  MPI_Gatherv(rho.data(), local_count, MPI_DOUBLE, rho_global.data(), recvcounts.data(),
              displs.data(), MPI_DOUBLE, 0, comm);
  if (rank == 0) {
    MPI_Send(rho_global.data(), rho_global.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }

  return {};
}
