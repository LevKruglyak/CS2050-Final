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

  auto sim = std::make_unique<Simulation>(10.0, 20.0, 100000);

  printf("loaded simulation\n");
  for (int i = 0; i < 10000; i++) {
    sim->update();

    if (i % 10 == 0) {
      std::ostringstream oss;
      oss << argv[1] << "_" << std::setw(2) << std::setfill('0') << i / 10
          << ".png";
      std::vector<double> density(sim->density.size());
      for (uint i = 0; i < sim->resolution * sim->resolution; i++) {
        density[i] = (sim->density[i] / sim->total_mass) *
                     (sim->resolution * sim->resolution) * 0.2;
      }
      WriteGrayscalePNG(density, sim->resolution, oss.str().c_str());
      printf("updated %d\n", i);
    }
  }
  return 0;
}
