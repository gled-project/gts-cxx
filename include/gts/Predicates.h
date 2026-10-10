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

// gts_point_orientation_3d_sos: the sign of orient3d(a, b, c, d), never
// zero. Ties are broken by Simulation of Simplicity, which needs a fixed
// order of the points: rank gives each point's place in it. gts-cxx ranks
// vertices by index; see QUIRKS.md.
int orient3d_sos(const Vec3* p[4], const unsigned rank[4]);

} // namespace gts
