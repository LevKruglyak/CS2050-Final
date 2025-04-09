#pragma once

#include <glm/glm.hpp>

using vec2 = glm::dvec2;

struct Particle {
  vec2 pos = vec2();
  vec2 vel = vec2();
  vec2 acc = vec2();
  double mass = 0.0;
};
