#pragma once

#include <cmath>           // std::sin, std::cosh, …
#include <glm/common.hpp>  // fract, mix, clamp, abs
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>  // glm::pi, etc.

/* ------------------------------------------------------------------------ */
/*  Small complex‑number helpers (z = x + i y, stored in vec2  = (x,y))     */
/* ------------------------------------------------------------------------ */
namespace cx {

using glm::vec2;
using glm::vec3;

/* real & imaginary parts */
inline float re(const vec2& z) {
  return z.x;
}
inline float im(const vec2& z) {
  return z.y;
}

/*  a + ib   +   c + id */
inline vec2 add(const vec2& a, const vec2& b) {
  return a + b;
}
/*  a + ib   −   c + id */
inline vec2 sub(const vec2& a, const vec2& b) {
  return a - b;
}

/*  (a+ib)(c+id) */
inline vec2 mul(const vec2& a, const vec2& b) {
  return vec2(a.x * b.x - a.y * b.y,   // real
              a.x * b.y + a.y * b.x);  // imag
}

/*  (a+ib)/(c+id) */
inline vec2 div(const vec2& a, const vec2& b) {
  const float denom = b.x * b.x + b.y * b.y;
  return vec2((a.x * b.x + a.y * b.y) / denom, (a.y * b.x - a.x * b.y) / denom);
}

/*  sin(a+ib)  = sin a cosh b + i cos a sinh b */
inline vec2 sin(const vec2& z) {
  return vec2(std::sin(z.x) * std::cosh(z.y), std::cos(z.x) * std::sinh(z.y));
}

/*  cos(a+ib)  = cos a cosh b − i sin a sinh b */
inline vec2 cos(const vec2& z) {
  return vec2(std::cos(z.x) * std::cosh(z.y), -std::sin(z.x) * std::sinh(z.y));
}

/*  tan z = sin z / cos z */
inline vec2 tan(const vec2& z) {
  return div(sin(z), cos(z));
}

/*  log|z| + i arg(z),   principal branch (arg ∈ (–π, π]) */
inline vec2 log(const vec2& z) {
  const float r = glm::length(z);      // |z|
  float theta = std::atan2(z.y, z.x);  // (–π, π]
  return vec2(std::log(r), theta);
}

/*  polar form (r, θ) */
inline vec2 as_polar(const vec2& z) {
  return vec2(glm::length(z), std::atan2(z.y, z.x));  // (r,θ)
}

/*  z^p  where  p is real */
inline vec2 pow(const vec2& z, float p) {
  auto pol = as_polar(z);  // (r,θ)
  float r = std::pow(pol.x, p);
  float ang = pol.y * p;
  return r * vec2(std::cos(ang), std::sin(ang));
}

/*  scaled argument on [0,2]  (utility from the original shader) */
inline float scaledArg(const vec2& z) {
  return std::atan2(z.y, z.x) / glm::pi<float>() + 1.0f;
}

}  // namespace cx
/* ------------------------------------------------------------------------ */
/*  HSV → RGB  (same 1‑liner scheme as the GLSL original)                   */
/* ------------------------------------------------------------------------ */
inline glm::vec3 hsv2rgb(const glm::vec3& hsv) {
  const glm::vec4 K(1.0f, 2.0f / 3.0f, 1.0f / 3.0f, 3.0f);
  const glm::vec3 p = glm::abs(glm::fract(glm::vec3(hsv.x) + glm::vec3(K.x, K.y, K.z)) * 6.0f - glm::vec3(K.w));
  return hsv.z * glm::mix(glm::vec3(K.x), glm::clamp(p - glm::vec3(K.x), 0.0f, 1.0f), hsv.y);
}

inline float hdrTone(float r, float exposure = 1.0f) {
  return 1.0f - std::exp(-exposure * r);
}

inline glm::vec3 complexColour(const glm::vec2& z) {
  const float hue = cx::log(z).y / (2.0f * glm::pi<float>());  // 0‥1
  const float val = hdrTone(glm::length(z), 0.01);
  const glm::vec3 hsv(hue, 1.0f, val);
  return hsv2rgb(hsv);
}
