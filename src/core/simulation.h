#pragma once

#include "glm/common.hpp"
#include <cmath>
#include <fmt/core.h>
#include <fmt/format.h>
#include <glm/glm.hpp>
#include <omp.h>
#include <vector>

using vec2 = glm::dvec2;

struct Particle {
  vec2 pos = vec2();
  vec2 vel = vec2();
  vec2 acc = vec2();
  double mass = 0.0;
};

class Simulation {
public:
  int resolution;
  size_t N;
  double total_mass;

  double grav_constant = 0.01;
  double softening = 0.001;
  double dt = 0.01;
  double radius = 10.0;
  double particles_per_cell = 10;

  std::vector<Particle> particles;
  std::vector<double> density;

  inline vec2 periodic_wrap(vec2 v) {
    return v - glm::round(v / radius) * radius;
  }

  inline vec2 accel(vec2 target, vec2 position) {
    vec2 diff = periodic_wrap(target - position);

    double dist2 = glm::dot(diff, diff) + softening * softening;
    double inv_r3 = 1.0 / (dist2 * std::sqrt(dist2));

    return grav_constant * diff * inv_r3;
  }

  Simulation(std::vector<Particle> particles) : particles(particles) {
    N = particles.size();
    resolution = std::ceil(sqrt((double)N / particles_per_cell));
    density.resize(resolution * resolution, 0);
    total_mass = 0.0;
    for (const auto &particle : particles) {
      total_mass += particle.mass;
    }
    updateDensity();
  }

  void updateDensity() {
    // Update density
    std::fill(density.begin(), density.end(), 0.0f);
    for (size_t i = 0; i < particles.size(); i++) {
      float dX = (particles[i].pos.x / radius + 0.5) * resolution;
      float dY = (particles[i].pos.y / radius + 0.5) * resolution;
      int gridX = (int)dX;
      int gridY = (int)dY;
      dX -= gridX;
      dY -= gridY;

      auto wrapIndex = [&](int i) {
        i = i % resolution;
        if (i < 0)
          i += resolution;
        return i;
      };

      int nX = wrapIndex(dX < 0.5 ? gridX - 1 : gridX + 1);
      int nY = wrapIndex(dY < 0.5 ? gridY - 1 : gridY + 1);

      float oX = abs(0.5f - dX);
      float oY = abs(0.5f - dY);
      float weight00 = (1.0f - oX) * (1.0f - oY);
      float weight10 = oX * (1.0f - oY);
      float weight01 = (1.0f - oX) * oY;
      float weight11 = oX * oY;

      density[gridY * resolution + gridX] += weight00 * particles[i].mass;
      density[gridY * resolution + nX] += weight10 * particles[i].mass;
      density[nY * resolution + gridX] += weight01 * particles[i].mass;
      density[nY * resolution + nX] += weight11 * particles[i].mass;
    }
  }

  void update() {
#pragma omp parallel for schedule(dynamic)
    for (uint i = 0; i < N; i++) {
      vec2 p = particles[i].pos;
      vec2 v = particles[i].vel;
      vec2 a = particles[i].acc;

      // Update position
      vec2 np = p + v * dt + a * (dt * dt * 0.5);

      // Calculate forces
      vec2 na = vec2();
      // #pragma omp parallel for reduction(+ : na)
      for (uint j = 0; j < N; j++) {
        na +=
            (i != j) ? particles[j].mass * accel(particles[j].pos, np) : vec2();
      }

      particles[i] = Particle{
          .pos = periodic_wrap(np),
          .vel = v + (a + na) * (dt * 0.5),
          .acc = na,
          .mass = particles[i].mass,
      };
    }

    updateDensity();
  }
};
