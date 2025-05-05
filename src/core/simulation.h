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
  FlatBHTree bh;

  inline vec2 wrap(vec2 v) const { return v - glm::round(v / radius) * radius; }

  inline vec2 accel(vec2 target, vec2 position) const {
    vec2 diff = wrap(target - position);
    double dist2 = glm::dot(diff, diff) + softening * softening;
    double inv_r3 = 1.0 / (dist2 * std::sqrt(dist2));
    return G * diff * inv_r3;
  }

  vec2 bh_accel_id(FlatBHTree::Id nid, const Particle &p) const {
    const auto &arena = bh.nodes_ref();
    const auto &node = arena[nid];
    if (node.mass <= 0.0)
      return vec2(0);

    const auto invalid = FlatBHTree::InvalidId;
    bool isLeaf = (node.child[0] == invalid);

    if (isLeaf) {
      vec2 a(0);
      for (uint32_t j = 0; j < node.bodyCnt; ++j) {
        const Particle &other = bh.bodies_ref()[node.firstBody + j];

        vec2 diff = wrap(other.p - p.p);
        if (diff.x == 0.0 && diff.y == 0.0)
          continue;

        double dist2 = glm::dot(diff, diff) + softening * softening;
        double inv_r3 = 1.0 / (dist2 * std::sqrt(dist2));
        a += G * diff * inv_r3 * other.m;
      }
      return a;
    }

    double s = node.bounds.size();
    double d = glm::length(wrap(node.com - p.p));

    if ((s / d) < bh_theta) {
      vec2 diff = wrap(node.com - p.p);
      double dist2 = glm::dot(diff, diff) + softening * softening;
      double inv_r3 = 1.0 / (dist2 * std::sqrt(dist2));
      return G * diff * inv_r3 * node.mass;
    } else {
      vec2 a(0);
      for (int k = 0; k < 4; ++k)
        if (node.child[k] != invalid)
          a += bh_accel_id(node.child[k], p);
      return a;
    }
  }

  Simulation(double radius, double total_mass, uint N)
      : N(N), total_mass(total_mass), radius(radius) {
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<double> sample(-radius / 2, radius / 2);

    ps.reserve(N);
    for (uint i = 0; i < N; ++i)
      ps.push_back(Particle{vec2(sample(rng), sample(rng)), vec2(0.0),
                            vec2(0.0), total_mass / N});

    finish_init();
  }

  explicit Simulation(std::vector<Particle> particles)
      : ps(std::move(particles)) {
    N = ps.size();
    total_mass = 0.0;
    for (const auto &b : ps)
      total_mass += b.m;
    finish_init();
  }

  void update() {
    bhUpdatePositions();
    updateDensity();
  }

private:
  void finish_init() {
    resolution = static_cast<int>(
        std::ceil(std::sqrt(static_cast<double>(N) / particles_per_cell)));
    density.resize(resolution * resolution, 0);
    updateDensity();
    updateBHTree();
  }

  void updateDensity() {
    std::fill(density.begin(), density.end(), 0.0);

    for (const auto &part : ps) {
      double dX = (part.p.x / radius + 0.5) * resolution;
      double dY = (part.p.y / radius + 0.5) * resolution;
      int gx = static_cast<int>(dX), gy = static_cast<int>(dY);
      dX -= gx;
      dY -= gy;

      auto wrapIdx = [this](int i) {
        i %= resolution;
        return (i < 0) ? i + resolution : i;
      };

      int nx = wrapIdx(dX < 0.5 ? gx - 1 : gx + 1);
      int ny = wrapIdx(dY < 0.5 ? gy - 1 : gy + 1);

      double oX = std::abs(0.5 - dX), oY = std::abs(0.5 - dY);
      double w00 = (1 - oX) * (1 - oY), w10 = oX * (1 - oY),
             w01 = (1 - oX) * oY, w11 = oX * oY;

      density[gy * resolution + gx] += w00 * part.m;
      density[gy * resolution + nx] += w10 * part.m;
      density[ny * resolution + gx] += w01 * part.m;
      density[ny * resolution + nx] += w11 * part.m;
    }
  }

  void updateBHTree() {
    AABB global({-radius / 2, -radius / 2}, {radius / 2, radius / 2});
    bh.build(ps, global);
  }

  void naiveUpdatePositions() {
#pragma omp parallel for schedule(dynamic)
    for (uint i = 0; i < N; ++i) {
      vec2 p = ps[i].p;
      vec2 v = ps[i].v;
      vec2 a = ps[i].a;

      vec2 np = wrap(p + v * dt + a * (dt * dt * 0.5));

      vec2 na(0);
      for (uint j = 0; j < N; ++j)
        if (i != j)
          na += ps[j].m * accel(ps[j].p, np);

      ps[i].p = np;
      ps[i].v = v + (a + na) * (dt * 0.5);
      ps[i].a = na;
    }
  }

  void bhUpdatePositions() {
    updateBHTree();

#pragma omp parallel for schedule(dynamic)
    for (uint i = 0; i < N; ++i) {
      vec2 p_next = wrap(ps[i].p + ps[i].v * dt + ps[i].a * (dt * dt * 0.5));
      vec2 a_next = bh_accel_id(bh.rootId(), ps[i]);

      ps[i].p = p_next;
      ps[i].v = ps[i].v + (ps[i].a + a_next) * (dt * 0.5);
      ps[i].a = a_next;
    }
  }
};
