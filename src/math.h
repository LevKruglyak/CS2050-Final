#pragma once

#include <cmath>

struct vec3 {
  double x;
  double y;
  double z;

  inline vec3() : x(0.0), y(0.0), z(0.0) {}
  inline vec3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}

  inline vec3 operator+(const vec3 &rhs) const {
    return vec3(x + rhs.x, y + rhs.y, z + rhs.z);
  }

  inline vec3 operator-(const vec3 &rhs) const {
    return vec3(x - rhs.x, y - rhs.y, z - rhs.z);
  }

  inline vec3 &operator+=(const vec3 &rhs) {
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
    return *this;
  }

  inline vec3 &operator-=(const vec3 &rhs) {
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
    return *this;
  }

  inline vec3 operator*(double s) const { return vec3(x * s, y * s, z * s); }
  inline vec3 operator/(double s) const { return vec3(x / s, y / s, z / s); }

  inline vec3 &operator*=(double s) {
    x *= s;
    y *= s;
    z *= s;
    return *this;
  }

  inline vec3 &operator/=(double s) {
    x /= s;
    y /= s;
    z /= s;
    return *this;
  }

  inline double dot(const vec3 &rhs) const {
    return x * rhs.x + y * rhs.y + z * rhs.z;
  }

  inline vec3 cross(const vec3 &rhs) const {
    return vec3(y * rhs.z - z * rhs.y, z * rhs.x - x * rhs.z,
                x * rhs.y - y * rhs.x);
  }

  inline double norm() const { return std::sqrt(x * x + y * y + z * z); }
  inline double norm2() const { return x * x + y * y + z * z; }

  inline vec3 normalized() const {
    double len = norm();
    return (len > 0.0) ? (*this / len) : vec3();
  }
};

inline vec3 operator*(double s, const vec3 &v) {
  return vec3(v.x * s, v.y * s, v.z * s);
}
