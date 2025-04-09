#pragma once

#include "bvh.h"
#include "common.h"
#include "glm/common.hpp"
#include <cmath>
#include <fmt/core.h>
#include <fmt/format.h>
#include <omp.h>
#include <random>
#include <vector>

class Simulation {
public:
  int resolution;
  size_t N;
  double total_mass;

  double G = 0.01;
  double softening = 0.001;
  double dt = 0.01;
  double radius = 10.0;
  double particles_per_cell = 10;
  double bh_theta = 0.25;

  std::vector<Particle> particles;
  std::vector<double> density;
  std::unique_ptr<BHTree> bh;

  inline vec2 wrap(vec2 v) const { return v - glm::round(v / radius) * radius; }

  inline vec2 accel(vec2 target, vec2 position) const {
    vec2 diff = wrap(target - position);

    double dist2 = glm::dot(diff, diff) + softening * softening;
    double inv_r3 = 1.0 / (dist2 * std::sqrt(dist2));

    return G * diff * inv_r3;
  }

  vec2 bh_accel(const BHTree *node, const Particle &p) const {
    if (node->totalMass <= 0.0)
      return vec2(0);

    vec2 a(0);
    if (!node->children[0] && node->bodies.size() > 0) {
      for (auto &leaf : node->bodies) {
        if (&leaf == &p)
          continue;

        vec2 diff = wrap(leaf.pos - p.pos);
        double dist2 = glm::dot(diff, diff) + softening * softening;
        double inv_r3 = 1.0 / (dist2 * std::sqrt(dist2));
        a += G * diff * inv_r3 * leaf.mass;
      }
      return a;
    } else {
      double s = node->boundary.size();
      double d = glm::length(wrap(node->centerOfMass - p.pos));
      if ((s / d) < bh_theta) {
        vec2 diff = wrap(node->centerOfMass - p.pos);
        double dist2 = glm::dot(diff, diff) + softening * softening;
        double inv_r3 = 1.0 / (dist2 * std::sqrt(dist2));
        return G * diff * inv_r3 * node->totalMass;
      } else {
        for (auto &child : node->children) {
          if (child) {
            a += bh_accel(child.get(), p);
          }
        }
        return a;
      }
    }
  }

  Simulation(double radius, double total_mass, uint N)
      : N(N), total_mass(total_mass), radius(radius) {
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<double> sample(-radius / 2, radius / 2);

    for (uint i = 0; i < N; i++) {
      vec2 x = vec2(sample(rng), sample(rng));

      particles.push_back(Particle{.pos = x, .mass = total_mass / N});
    }

    resolution = std::ceil(sqrt((double)N / particles_per_cell));
    density.resize(resolution * resolution, 0);
    total_mass = 0.0;
    for (const auto &particle : particles) {
      total_mass += particle.mass;
    }
    updateDensity();
    updateBHTree();
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
    updateBHTree();
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

  void naiveUpdatePositions() {
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
          .pos = wrap(np),
          .vel = v + (a + na) * (dt * 0.5),
          .acc = na,
          .mass = particles[i].mass,
      };
    }
  }

  void updateBHTree() {
    AABB globalRegion(vec2(-radius / 2, -radius / 2),
                      vec2(radius / 2, radius / 2));

    bh = std::make_unique<BHTree>(globalRegion);
    for (auto &b : particles) {
      bh->insert(b);
    }
    bh->computeMassDistribution();
  }

  void bhUpdatePositions() {
    updateBHTree();

#pragma omp parallel for schedule(dynamic)
    for (uint i = 0; i < N; i++) {
      vec2 p = particles[i].pos;
      vec2 v = particles[i].vel;
      vec2 a = particles[i].acc;

      vec2 np = p + v * dt + a * (dt * dt * 0.5);
      vec2 na = bh_accel(bh.get(), particles[i]);
      particles[i] = Particle{
          .pos = wrap(np),
          .vel = v + (a + na) * (dt * 0.5),
          .acc = na,
          .mass = particles[i].mass,
      };
    }
  }

  void update() {
    bhUpdatePositions();
    updateDensity();
  }
};
