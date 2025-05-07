// #include "simulation.h"
//
// #include <fftw3-mpi.h>
// #include <mpi.h>
// #include <stb_image_write.h>
// #include <nlohmann/json.hpp>
//
// #include <cstdlib>
// #include <fstream>
// #include <iostream>
// #include <string>
// #include <vector>
//
// // -----------------------------------------------------------------------------
// // Simple flag parser
// struct CLIOptions {
//   std::string config_file;
//   int iterations = 1000;
//   int save_frequency = 100;
//   std::string out_pref = "density";
// };
//
// CLIOptions parse_args(int argc, char** argv) {
//   CLIOptions o;
//   for (int i = 1; i < argc; ++i) {
//     std::string s = argv[i];
//     if (s == "--config" && i + 1 < argc)
//       o.config_file = argv[++i];
//     else if (s == "--iterations" && i + 1 < argc)
//       o.iterations = std::atoi(argv[++i]);
//     else if (s == "--save-freq" && i + 1 < argc)
//       o.save_frequency = std::atoi(argv[++i]);
//     else if (s == "--out-prefix" && i + 1 < argc)
//       o.out_pref = argv[++i];
//     else {
//       std::cerr << "Unknown flag " << s << "\n";
//       std::exit(1);
//     }
//   }
//   return o;
// }
//
// // -----------------------------------------------------------------------------
// // Load JSON config into Simulation::Params
// void load_config_json(const std::string& path, Simulation::Params& p) {
//   std::ifstream in(path);
//   if (!in) {
//     std::cerr << "Failed to open config: " << path << "\n";
//     std::exit(1);
//   }
//   nlohmann::json j;
//   in >> j;
//
//   if (j.contains("GRAVITY"))
//     p.GRAVITY = j["GRAVITY"].get<float>();
//   if (j.contains("SOFTENING"))
//     p.SOFTENING = j["SOFTENING"].get<float>();
//   if (j.contains("TIMESTEP"))
//     p.TIMESTEP = j["TIMESTEP"].get<float>();
//   if (j.contains("RADIUS"))
//     p.RADIUS = j["RADIUS"].get<float>();
//   if (j.contains("MASS"))
//     p.MASS = j["MASS"].get<float>();
//   if (j.contains("PERLIN_NOISE_SCALE"))
//     p.PERLIN_NOISE_SCALE = j["PERLIN_NOISE_SCALE"].get<float>();
//   if (j.contains("PERLIN_NOISE_PERTURBATION"))
//     p.PERLIN_NOISE_PERTURBATION = j["PERLIN_NOISE_PERTURBATION"].get<float>();
//   if (j.contains("PERLIN_NOISE_OCTAVES"))
//     p.PERLIN_NOISE_OCTAVES = j["PERLIN_NOISE_OCTAVES"].get<int>();
//   if (j.contains("PARTICLES_PER_CELL"))
//     p.PARTICLES_PER_CELL = j["PARTICLES_PER_CELL"].get<float>();
//   if (j.contains("RESOLUTION"))
//     p.RESOLUTION = j["RESOLUTION"].get<ptrdiff_t>();
// }
//
// // -----------------------------------------------------------------------------
// // Write a grayscale PNG [0..1] → [0..255]
// void write_png(const std::string& fname, const std::vector<float>& img, int W, int H) {
//   // convert float→uchar
//   std::vector<unsigned char> out(W * H);
//   for (int i = 0; i < W * H; ++i) {
//     float v = img[i];
//     if (v < 0)
//       v = 0;
//     if (v > 1)
//       v = 1;
//     out[i] = static_cast<unsigned char>(v * 255.0f);
//   }
//   // 1 channel, stride = W
//   stbi_write_png(fname.c_str(), W, H, 1, out.data(), W);
// }
//
// // -----------------------------------------------------------------------------
int main(int argc, char** argv) {}
//   // 1) parse CLI
//   auto opts = parse_args(argc, argv);
//
//   // 2) MPI + FFTW init
//   int provided;
//   MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
//   fftw_mpi_init();
//   fftw_init_threads();
//
//   int rank;
//   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//   // 3) load JSON config on rank 0
//   Simulation::Params params;
//   if (rank == 0 && !opts.config_file.empty()) {
//     load_config_json(opts.config_file, params);
//   }
//   // 4) broadcast to all ranks
//   MPI_Bcast(&params, sizeof(params), MPI_BYTE, 0, MPI_COMM_WORLD);
//
//   // 5) construct simulation on every rank
//   Simulation sim(params);
//
//   // 6) main timestep loop
//   for (int step = 1; step <= opts.iterations; ++step) {
//     sim.timestep();
//
//     if (step % opts.save_frequency == 0) {
//       // gather the global density on rank 0
//       auto global_rho = sim.gather_rho();
//       if (rank == 0) {
//         int W = params.RESOLUTION, H = params.RESOLUTION;
//         // normalize to [0,1]
//         float maxr = *std::max_element(global_rho.begin(), global_rho.end());
//         std::vector<float> img(W * H);
//         for (int i = 0; i < W * H; ++i)
//           img[i] = global_rho[i] / maxr;
//
//         std::ostringstream filename;
//         filename << opts.out_pref << "_" << step << ".png";
//         write_png(filename.str(), img, W, H);
//         std::cout << "Wrote " << filename.str() << "\n";
//       }
//     }
//   }
//
//   MPI_Finalize();
//   return 0;
// }
