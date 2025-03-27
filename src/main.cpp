#include <cmath>
#include <highfive/highfive.hpp>
#include <iostream>
#include <omp.h>
#include <vector>

#include "math.h"

static const double GRAV = 4278512.6516; // [Mpc] / [Msun] [yr]^2
static const double SOFTENING = 0.1;     // [Mpc]

struct particle {
  vec3 pos;
  vec3 vel;
  vec3 acc;
  double mass;
};

vec3 accel(vec3 target, vec3 position) {
  vec3 diff = target - position;
  float dist2 = diff.dot(diff) + SOFTENING * SOFTENING;
  return GRAV * diff / pow(dist2, 1.5);
}

class Simulation {
  double dt = 1.0;
  std::vector<particle> particles;

public:
  void update() {
#pragma omp parallel for
    for (uint i = 0; i < particles.size(); i++) {
      vec3 p = particles[i].pos;
      vec3 v = particles[i].vel;
      vec3 a = particles[i].acc;

      // Update position
      p = p + v * dt + a * (dt * dt * 0.5);

      // Calculate forces
      vec3 na = vec3();
      for (uint j = 0; j < particles.size(); j++) {
        na += particles[i].mass * accel(particles[j].pos, p);
      }
      na *= particles[i].mass;

      particles[i].pos = p;
      particles[i].vel = v + (a + na) * (dt * 0.5);
      particles[i].acc = na;
    }
  }
};

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  std::cout << "Hello, World!" << std::endl;
  return 0;
}
