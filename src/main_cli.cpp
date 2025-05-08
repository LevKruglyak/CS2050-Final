#include <mpi.h>
#include "simulation.h"

#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <sstream>
#include <string>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

struct CLIOptions {
  std::string config_file;
  int iterations = 1000;
  int save_frequency = 100;
  std::string out_pref = "density";
};

CLIOptions parse_args(int argc, char** argv) {
  CLIOptions opts;
  for (int i = 1; i < argc; ++i) {
    std::string s = argv[i];
    if (s == "--config" && i + 1 < argc)
      opts.config_file = argv[++i];
    else if (s == "--iterations" && i + 1 < argc)
      opts.iterations = std::atoi(argv[++i]);
    else if (s == "--save-freq" && i + 1 < argc)
      opts.save_frequency = std::atoi(argv[++i]);
    else if (s == "--out-prefix" && i + 1 < argc)
      opts.out_pref = argv[++i];
    else {
      std::cerr << "Unknown flag " << s << "\n";
      std::exit(1);
    }
  }
  return opts;
}

void load_config_json(const std::string& path, Simulation::Params& p) {
  std::ifstream in(path);
  if (!in) {
    std::cerr << "Failed to open config: " << path << "\n";
    std::exit(1);
  }
  nlohmann::json j;
  in >> j;

  if (j.contains("GRAVITY"))
    p.GRAVITY = j["GRAVITY"].get<float>();
  if (j.contains("SOFTENING"))
    p.SOFTENING = j["SOFTENING"].get<float>();
  if (j.contains("TIMESTEP"))
    p.TIMESTEP = j["TIMESTEP"].get<float>();
  if (j.contains("RADIUS"))
    p.RADIUS = j["RADIUS"].get<float>();
  if (j.contains("MASS"))
    p.MASS = j["MASS"].get<float>();
  if (j.contains("PERLIN_NOISE_SCALE"))
    p.PERLIN_NOISE_SCALE = j["PERLIN_NOISE_SCALE"].get<float>();
  if (j.contains("PERLIN_NOISE_PERTURBATION"))
    p.PERLIN_NOISE_PERTURBATION = j["PERLIN_NOISE_PERTURBATION"].get<float>();
  if (j.contains("PERLIN_NOISE_OCTAVES"))
    p.PERLIN_NOISE_OCTAVES = j["PERLIN_NOISE_OCTAVES"].get<int>();
  if (j.contains("PARTICLES_PER_CELL"))
    p.PARTICLES_PER_CELL = j["PARTICLES_PER_CELL"].get<float>();
  if (j.contains("RESOLUTION"))
    p.RESOLUTION = j["RESOLUTION"].get<ptrdiff_t>();
  if (j.contains("USE_SCALE_FACTOR"))
    p.USE_SCALE_FACTOR = j["USE_SCALE_FACTOR"].get<bool>();
}

void write_png(const std::string& fname, const std::vector<float>& img, int W, int H) {
  std::vector<unsigned char> out(W * H);
  for (int i = 0; i < W * H; ++i) {
    float v = img[i];
    if (v < 0)
      v = 0;
    if (v > 1)
      v = 1;
    out[i] = static_cast<unsigned char>(v * 255.0f);
  }

  stbi_write_png(fname.c_str(), W, H, 1, out.data(), W);
}

int main(int argc, char** argv) {
  auto opts = parse_args(argc, argv);

  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  fftw_init_threads();
  fftw_mpi_init();

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  Simulation::Params params;
  if (rank == 0 && !opts.config_file.empty()) {
    load_config_json(opts.config_file, params);
  }

  MPI_Bcast(&params, sizeof(params), MPI_BYTE, 0, MPI_COMM_WORLD);

  auto density_hdr = [&params](float input) {
    input *= (0.3 * params.RADIUS * params.RADIUS / params.MASS);
    return 1.0 - exp(-input * 3.0);
  };

  double average_total = 0;
  double average_mass = 0;
  double average_solve = 0;
  double average_update = 0;
  {
    Simulation sim(params);
    for (int step = 0; step < opts.iterations; ++step) {
      auto profile = sim.timestep();

      if (rank == 0) {
        average_total += profile.total_time / opts.iterations;
        average_mass += (profile.mass_local_accum + profile.mass_halo_exchange) / opts.iterations;
        average_solve +=
            (profile.fft_forward + profile.fft_backward + profile.spectral_solve) / opts.iterations;
        average_update +=
            (profile.update_positions + profile.reassign_particles + profile.force_halo_exchange) /
            opts.iterations;
      }

      if (opts.save_frequency != 0 && step % opts.save_frequency == 0) {
        // gather the global density on rank 0
        auto global_rho = sim.gather_rho();
        if (rank == 0) {
          int W = params.RESOLUTION, H = params.RESOLUTION;
          // normalize to [0,1]
          std::vector<float> img(W * H);
          for (int i = 0; i < W * H; ++i)
            img[i] = density_hdr(global_rho[i]);

          std::ostringstream filename;
          filename << opts.out_pref << "_" << step / opts.save_frequency << ".png";
          write_png(filename.str(), img, W, H);
          // std::cout << "Wrote " << filename.str() << "\n";
        }
      }
    }
  }

  if (rank == 0) {
    std::cout                     //
        << average_total << " "   //
        << average_mass << " "    //
        << average_solve << " "   //
        << average_update << " "  //
        << std::endl;             //
  }

  fftw_cleanup_threads();
  MPI_Finalize();
}
