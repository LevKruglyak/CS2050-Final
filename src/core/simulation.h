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
  double particles_per_cell = 3;
  double bh_theta = 0.25;

  std::vector<Particle> ps;
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

        vec2 diff = wrap(leaf.p - p.p);
        double dist2 = glm::dot(diff, diff) + softening * softening;
        double inv_r3 = 1.0 / (dist2 * std::sqrt(dist2));
        a += G * diff * inv_r3 * leaf.m;
      }
      return a;
    } else {
      double s = node->boundary.size();
      double d = glm::length(wrap(node->centerOfMass - p.p));
      if ((s / d) < bh_theta) {
        vec2 diff = wrap(node->centerOfMass - p.p);
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

      ps.push_back(Particle{.p = x, .m = total_mass / N});
    }

    resolution = std::ceil(sqrt((double)N / particles_per_cell));
    density.resize(resolution * resolution, 0);
    total_mass = 0.0;
    for (const auto &particle : ps) {
      total_mass += particle.m;
    }
    updateDensity();
    updateBHTree();
  }

  Simulation(std::vector<Particle> particles) : ps(particles) {
    N = particles.size();
    resolution = std::ceil(sqrt((double)N / particles_per_cell));
    density.resize(resolution * resolution, 0);
    total_mass = 0.0;
    for (const auto &particle : particles) {
      total_mass += particle.m;
    }
    updateDensity();
    updateBHTree();
  }

  void updateDensity() {
    // Update density
    std::fill(density.begin(), density.end(), 0.0f);
    for (size_t i = 0; i < ps.size(); i++) {
      float dX = (ps[i].p.x / radius + 0.5) * resolution;
      float dY = (ps[i].p.y / radius + 0.5) * resolution;
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

      density[gridY * resolution + gridX] += weight00 * ps[i].m;
      density[gridY * resolution + nX] += weight10 * ps[i].m;
      density[nY * resolution + gridX] += weight01 * ps[i].m;
      density[nY * resolution + nX] += weight11 * ps[i].m;
    }
  }

  void naiveUpdatePositions() {
#pragma omp parallel for schedule(dynamic)
    for (uint i = 0; i < N; i++) {
      vec2 p = ps[i].p;
      vec2 v = ps[i].v;
      vec2 a = ps[i].a;

      // Update position
      vec2 np = p + v * dt + a * (dt * dt * 0.5);

      // Calculate forces
      vec2 na = vec2();
      // #pragma omp parallel for reduction(+ : na)
      for (uint j = 0; j < N; j++) {
        na += (i != j) ? ps[j].m * accel(ps[j].p, np) : vec2();
      }

      ps[i] = Particle{
          .p = wrap(np),
          .v = v + (a + na) * (dt * 0.5),
          .a = na,
          .m = ps[i].m,
      };
    }
  }

  void updateBHTree() {
    AABB globalRegion(vec2(-radius / 2, -radius / 2),
                      vec2(radius / 2, radius / 2));

    bh = std::make_unique<BHTree>(globalRegion);
    for (auto &b : ps) {
      bh->insert(b);
    }
    bh->computeMassDistribution();
  }

  void bhUpdatePositions() {
    updateBHTree();

#pragma omp parallel for schedule(dynamic)
    for (uint i = 0; i < N; i++) {
      vec2 p = wrap(ps[i].p + ps[i].v * dt + ps[i].a * (dt * dt * 0.5));
      vec2 a = bh_accel(bh.get(), ps[i]);

      ps[i].p = p;
      ps[i].v = ps[i].v + (ps[i].a + a) * (dt * 0.5);
      ps[i].a = a;
    }
  }

  void update() {
    bhUpdatePositions();
    updateDensity();
  }
};
