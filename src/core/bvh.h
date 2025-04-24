#pragma once

#include "common.h"
#include <algorithm>
#include <cassert>
#include <limits>
#include <memory>
#include <vector>

struct AABB {
  vec2 min;
  vec2 max;

  AABB() {
    min = vec2(std::numeric_limits<double>::max());
    max = vec2(std::numeric_limits<double>::lowest());
  }

  AABB(const vec2 &mn, const vec2 &mx) : min(mn), max(mx) {}

  bool contains(const vec2 &p) const {
    return (p.x >= min.x && p.x <= max.x && p.y >= min.y && p.y <= max.y);
  }

  double size() const {
    double dx = max.x - min.x;
    double dy = max.y - min.y;
    return std::max(dx, dy);
  }

  AABB quadrant(int index) const {
    vec2 center = 0.5 * (min + max);
    AABB box;

    switch (index) {
    case 0:
      box.min = vec2(min.x, center.y);
      box.max = vec2(center.x, max.y);
      break;
    case 1:
      box.min = vec2(center.x, center.y);
      box.max = vec2(max.x, max.y);
      break;
    case 2:
      box.min = vec2(min.x, min.y);
      box.max = vec2(center.x, center.y);
      break;
    case 3:
      box.min = vec2(center.x, min.y);
      box.max = vec2(max.x, center.y);
      break;
    default:
      box = *this;
      break;
    }
    return box;
  }
};

class BHTree {
public:
  static constexpr int CAPACITY = 4;
  std::unique_ptr<BHTree> children[4] = {nullptr, nullptr, nullptr, nullptr};
  AABB boundary;
  std::vector<Particle> bodies;

  double totalMass = 0.0;
  vec2 centerOfMass = vec2(0);
  int numParticles = 0;

  BHTree(const AABB &region) : boundary(region) {}

  void insert(const Particle &b) {
    numParticles++;
    if (!boundary.contains(b.p))
      return;

    if (bodies.size() < CAPACITY && children[0] == nullptr) {
      bodies.push_back(b);
      return;
    }

    if (!children[0]) {
      subdivide();
      for (const auto &existingParticle : bodies)
        insertIntoChildren(existingParticle);
      bodies.clear();
    }

    insertIntoChildren(b);
  }

  void computeMassDistribution() {
    if (!children[0]) {
      totalMass = 0.0;
      centerOfMass = vec2(0);
      for (const auto &b : bodies) {
        totalMass += b.m;
        centerOfMass += b.m * b.p;
      }
      if (totalMass > 0.0)
        centerOfMass /= totalMass;
    } else {
      // Internal node
      totalMass = 0.0;
      centerOfMass = vec2(0);

      for (auto &child : children) {
        if (child) {
          child->computeMassDistribution();
          totalMass += child->totalMass;
          centerOfMass += child->totalMass * child->centerOfMass;
        }
      }
      if (totalMass > 0.0)
        centerOfMass /= totalMass;
    }
  }

private:
  void subdivide() {
    for (int i = 0; i < 4; i++) {
      AABB childRegion = boundary.quadrant(i);
      children[i] = std::make_unique<BHTree>(childRegion);
    }
  }

  void insertIntoChildren(const Particle &b) {
    for (int i = 0; i < 4; i++) {
      if (children[i]->boundary.contains(b.p)) {
        children[i]->insert(b);
        break;
      }
    }
  }
};
