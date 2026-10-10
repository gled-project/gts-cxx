// Geometry as values: points, segments, triangles and boxes, and the
// measures computed from them. Nothing here belongs to a mesh.
//
// Copyright (C) 1999 Stéphane Popinet (GTS, triangle.c, point.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include <cmath>
#include <limits>

namespace gts {

struct Vec3
{
  double x = 0, y = 0, z = 0;

  friend Vec3 operator+(Vec3 a, Vec3 b)   { return { a.x + b.x, a.y + b.y, a.z + b.z }; }
  friend Vec3 operator-(Vec3 a, Vec3 b)   { return { a.x - b.x, a.y - b.y, a.z - b.z }; }
  friend Vec3 operator*(double s, Vec3 a) { return { s * a.x, s * a.y, s * a.z }; }
  friend bool operator==(Vec3, Vec3) = default;
};

inline double dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

// The components in the order of GTS's gts_triangle_normal and
// gts_vector_cross, so that results agree to the last bit.
inline Vec3 cross(Vec3 a, Vec3 b)
{
  return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
}

inline double norm(Vec3 a) { return std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z); }

// gts_point_distance
inline double distance(Vec3 a, Vec3 b) { return norm(a - b); }
// gts_point_distance2
inline double distance2(Vec3 a, Vec3 b) { return dot(a - b, a - b); }

struct Segment
{
  Vec3 a, b;
};

struct Triangle
{
  Vec3 a, b, c;
};

// The perimeter of an equilateral triangle of unit area, so that quality()
// is 1 for an equilateral triangle. GTS's GOLDEN_PERIMETER.
inline constexpr double kGoldenPerimeter = 4.5590141139;

// (b - a) x (c - a): twice the area, along the normal of a, b, c taken
// counter-clockwise.
inline Vec3   normal(const Triangle& t)    { return cross(t.b - t.a, t.c - t.a); }
inline double area(const Triangle& t)      { return norm(normal(t)) / 2.; }
inline double perimeter(const Triangle& t) { return distance(t.a, t.b) + distance(t.a, t.c) + distance(t.b, t.c); }

// gts_triangle_quality: 1 for an equilateral triangle, 0 for a degenerate one.
inline double quality(const Triangle& t)
{
  const double p = perimeter(t);
  return p > 0.0 ? kGoldenPerimeter * std::sqrt(area(t)) / p : 0.0;
}

// An axis-aligned box. An empty box has min above max.
struct BBox
{
  Vec3 min { std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max() };
  Vec3 max { -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max() };

  bool empty() const { return min.x > max.x; }
  void add(Vec3 p)
  {
    min = { std::fmin(min.x, p.x), std::fmin(min.y, p.y), std::fmin(min.z, p.z) };
    max = { std::fmax(max.x, p.x), std::fmax(max.y, p.y), std::fmax(max.z, p.z) };
  }
};

} // namespace gts
