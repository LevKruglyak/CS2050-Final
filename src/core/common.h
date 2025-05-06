#pragma once

#include <glm/glm.hpp>

using vec2 = glm::dvec2;

struct Particle {
  vec2 p = vec2();
  vec2 v = vec2();
  vec2 a = vec2();
  double mass;
};
