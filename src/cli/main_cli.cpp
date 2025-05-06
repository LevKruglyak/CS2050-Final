#include "glm/gtc/constants.hpp"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "simulation.h"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

void WriteGrayscalePNG(const std::vector<double> &data, int resolution,
                       const char *filename) {
  std::vector<uint8_t> img(resolution * resolution);

  for (int i = 0; i < resolution * resolution; i++) {
    double val = data[i];
    val = std::max(0.0, std::min(1.0, val));
    img[i] = static_cast<uint8_t>(val * 255.0);
  }

  stbi_write_png(filename, resolution, resolution, 1, img.data(), resolution);
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  // std::unique_ptr<Simulation> sim;
  // if (strncmp(argv[1], "galaxy", 7) == 0) {
  //   std::random_device rd;
  //   std::mt19937 rng(rd());
  //
  //   std::uniform_real_distribution<double> angle_sample(0.0,
  //                                                       glm::tau<double>());
  //   std::uniform_real_distribution<double> radius_sample(0.0, 1.75);
  //
  //   double black_hole_mass = 1.0;
  //
  //   vec2 bcenter = vec2(0.01, 0.01);
  //   std::vector<Particle> particles;
  //   const uint num_particles = 50000;
  //
  //   particles.push_back(Particle{.p = bcenter, .m = black_hole_mass});
  //
  //   double grav_constant = 0.01;
  //   for (uint i = 0; i < num_particles; i++) {
  //     double radius = 1.0 / (2.0 - radius_sample(rng));
  //     double angle = angle_sample(rng);
  //
  //     double x = cos(angle);
  //     double y = sin(angle);
  //     double speed = sqrt(grav_constant * black_hole_mass / radius);
  //
  //     particles.push_back(Particle{.p = radius * vec2(x, y) + bcenter,
  //                                  .v = speed * vec2(-y, x),
  //                                  .m = 0.00001});
  //   }
  //   sim = std::make_unique<Simulation>(particles);
  // } else {
  //   sim = std::make_unique<Simulation>(10.0, 20.0, 10000);
  // }
  //
  // printf("loaded simulation\n");
  // for (int i = 0; i < 10000; i++) {
  //   sim->update();
  //
  //   if (i % 10 == 0) {
  //     std::ostringstream oss;
  //     oss << argv[2] << "_" << std::setw(2) << std::setfill('0') << i / 10
  //         << ".png";
  //     std::vector<double> density(sim->density.size());
  //     for (uint i = 0; i < sim->resolution * sim->resolution; i++) {
  //       density[i] = (sim->density[i] / sim->total_mass) *
  //                    (sim->resolution * sim->resolution) * 0.2;
  //     }
  //     WriteGrayscalePNG(density, sim->resolution, oss.str().c_str());
  //     printf("updated %d\n", i);
  //   }
  // }
  return 0;
}
