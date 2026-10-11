// Robust geometric predicates by Jonathan Richard Shewchuk, with adaptive
// precision arithmetic: the sign of the result is always exact.
//
// Copyright (C) 2026 Matevz Tadel
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "gts/Geometry.h"

namespace gts {

// Positive when a, b and c run counterclockwise in the plane, negative when
// clockwise, zero when they are collinear. The value is approximately twice
// the signed area of the triangle.
double orient2d(const double* a, const double* b, const double* c);

// Positive when d lies below the plane through a, b and c, where below is
// the side from which a, b and c appear clockwise; negative above; zero
// when the four points are coplanar. The value is approximately six times
// the signed volume of the tetrahedron. gts_point_orientation_3d.
double orient3d(const double* a, const double* b, const double* c, const double* d);

inline double orient3d(Vec3 a, Vec3 b, Vec3 c, Vec3 d)
{
  const double pa[3] = { a.x, a.y, a.z }, pb[3] = { b.x, b.y, b.z };
  const double pc[3] = { c.x, c.y, c.z }, pd[3] = { d.x, d.y, d.z };
  return orient3d(pa, pb, pc, pd);
}

// Positive when d lies inside the circle through a, b and c, negative
// outside, zero on it; a, b and c must run counterclockwise, or the sign is
// reversed. Only x and y are read.
double incircle(const double* a, const double* b, const double* c, const double* d);

// Positive when e lies inside the sphere through a, b, c and d, negative
// outside, zero on it; orient3d(a, b, c, d) must be positive, or the sign
// is reversed.
double insphere(const double* a, const double* b, const double* c, const double* d, const double* e);

// gts_point_orientation: orient2d() of the projections on the xy plane.
inline double orient2d(Vec3 a, Vec3 b, Vec3 c)
{
  const double pa[2] = { a.x, a.y }, pb[2] = { b.x, b.y }, pc[2] = { c.x, c.y };
  return orient2d(pa, pb, pc);
}

// gts_point_in_circle: incircle() of the projections of a, b, c and p on
// the xy plane: positive when p lies inside.
inline double in_circle(Vec3 p, Vec3 a, Vec3 b, Vec3 c)
{
  const double pa[2] = { a.x, a.y }, pb[2] = { b.x, b.y }, pc[2] = { c.x, c.y }, pp[2] = { p.x, p.y };
  return incircle(pa, pb, pc, pp);
}

// gts_point_in_sphere: insphere() of a, b, c, d and p: positive when p
// lies inside.
inline double in_sphere(Vec3 p, Vec3 a, Vec3 b, Vec3 c, Vec3 d)
{
  const double pa[3] = { a.x, a.y, a.z }, pb[3] = { b.x, b.y, b.z }, pc[3] = { c.x, c.y, c.z };
  const double pd[3] = { d.x, d.y, d.z }, pp[3] = { p.x, p.y, p.z };
  return insphere(pa, pb, pc, pd, pp);
}

// gts_point_orientation_sos: the sign of orient2d(a, b, c) of the
// projections on the xy plane, never zero, with ties broken by Simulation
// of Simplicity in the order given by rank, as orient3d_sos().
int orient2d_sos(const Vec3* p[3], const unsigned rank[3]);

// gts_point_orientation_3d_sos: the sign of orient3d(a, b, c, d), never
// zero. Ties are broken by Simulation of Simplicity, which needs a fixed
// order of the points: rank gives each point's place in it. gts-cxx ranks
// vertices by index; see QUIRKS.md.
int orient3d_sos(const Vec3* p[4], const unsigned rank[4]);

} // namespace gts
