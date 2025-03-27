#include <cmath>
#include <fmt/core.h>
#include <fmt/format.h>
#include <glm/glm.hpp>
#include <highfive/highfive.hpp>
#include <iostream>
#include <omp.h>
#include <vector>

static const double GRAV = 4278512.6516; // [Mpc] / [Msun] [yr]^2
static const double SOFTENING = 0.01;    // [Mpc]

using vec3 = glm::dvec3;

struct particle {
  vec3 pos = vec3();
  vec3 vel = vec3();
  vec3 acc = vec3();
  double mass = 0.0;
};

class Simulation {
  double dt = 1e-20;
  std::vector<particle> particles;

  vec3 accel(vec3 target, vec3 position) {
    vec3 diff = target - position;
    float dist2 = glm::dot(diff, diff) + SOFTENING * SOFTENING;
    return GRAV * diff / pow(dist2, 1.5);
  }

public:
  Simulation(std::vector<particle> particles) : particles(particles) {}

  void update() {
#pragma omp parallel for
    for (uint i = 0; i < particles.size(); i++) {
      vec3 p = particles[i].pos;
      vec3 v = particles[i].vel;
      vec3 a = particles[i].acc;
      double m = particles[i].mass;

      // Update position
      p = p + v * dt + a * (dt * dt * 0.5);

      // Calculate forces
      vec3 na = vec3();
      for (uint j = 0; j < particles.size(); j++) {
        na += accel(particles[j].pos, p);
      }
      na *= m;

      particles[i].pos = p;
      particles[i].vel = v + (a + na) * (dt * 0.5);
      particles[i].acc = na;
    }

    // Calculate total energy
    double energy = 0.0;
#pragma omp parallel for reduction(+ : energy)
    for (uint i = 0; i < particles.size(); i++) {
      vec3 p = particles[i].pos;
      vec3 v = particles[i].vel;
      vec3 a = particles[i].acc;
      double m = particles[i].mass;

      energy += 0.5 * m * glm::dot(v, v) - 0.5 * m * glm::dot(a, p);
    }

    std::cout << fmt::format("x-pos {} {}\n", particles[0].pos.x,
                             particles[1].pos.x);
  }
};
