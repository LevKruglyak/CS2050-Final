#pragma once

#include "glm/common.hpp"
#include <cmath>
#include <fmt/core.h>
#include <fmt/format.h>
#include <glm/glm.hpp>
#include <highfive/highfive.hpp>
#include <iostream>
#include <omp.h>
#include <vector>

using vec3 = glm::dvec3;

struct Particle {
  vec3 pos = vec3();
  vec3 vel = vec3();
  vec3 acc = vec3();
  double mass = 0.0;
};

struct SimulationParams {
  double grav_constant;
  double softening;
  double dt;
  double radius;
};

class Simulation {
public:
  std::vector<Particle> particles;

  std::vector<double> es;
  std::vector<double> kes;
  std::vector<double> pes;

  SimulationParams params;

  inline vec3 periodic_wrap(vec3 v) {
    return v - glm::round(v / params.radius) * params.radius;
  }

  inline vec3 accel(vec3 target, vec3 position) {
    vec3 diff = periodic_wrap(target - position);

    double dist2 = glm::dot(diff, diff) + params.softening * params.softening;
    double inv_r3 = 1.0 / (dist2 * std::sqrt(dist2));

    return params.grav_constant * diff * inv_r3;
  }

  Simulation(std::vector<Particle> particles, SimulationParams params)
      : particles(particles), params(params) {}

  void update() {
#pragma omp parallel for schedule(dynamic)
    for (uint i = 0; i < particles.size(); i++) {
      vec3 p = particles[i].pos;
      vec3 v = particles[i].vel;
      vec3 a = particles[i].acc;

      // Update position
      vec3 np = p + v * params.dt + a * (params.dt * params.dt * 0.5);

      // Calculate forces
      vec3 na = vec3();
#pragma omp parallel for reduction(+ : na)
      for (uint j = 0; j < particles.size(); j++) {
        na +=
            (i != j) ? particles[j].mass * accel(particles[j].pos, np) : vec3();
      }

      particles[i].pos = periodic_wrap(np);
      particles[i].vel = v + (a + na) * (params.dt * 0.5);
      particles[i].acc = na;
    }

    // Calculate total energy
    double ke = 0.0;
    double pe = 0.0;
#pragma omp parallel for reduction(+ : ke, pe)
    for (uint i = 0; i < particles.size(); i++) {
      vec3 v = particles[i].vel;
      vec3 a = particles[i].acc;
      double m = particles[i].mass;

      ke += 0.5 * m * glm::dot(v, v);
      pe -= 0.5 * m * glm::dot(a, a);
    }

    kes.push_back(ke);
    pes.push_back(pe);
    es.push_back(ke + pe);
  }
};
