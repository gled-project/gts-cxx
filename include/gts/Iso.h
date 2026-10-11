// Iso-surfaces: the surface f(x, y, z) = iso of a function sampled on a
// cartesian grid, by marching cubes (isosurface_cartesian) or by splitting
// each cube into tetrahedra (isosurface_tetra and its two variants).
//
// The function is given plane by plane, as in GTS: a PlaneFunc fills the
// values of one plane z = g.z of the grid. sample() makes one from a function
// of a point.
//
// Copyright (C) 1999 Stéphane Popinet (GTS: iso.c)
// Copyright (C) 2002 Gert Wollny (GTS: isotetra.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "gts/Geometrium.h"

#include <functional>
#include <vector>

namespace gts {

// GtsCartesianGrid: nx * ny * nz points, from (x, y, z) in steps of dx, dy,
// dz.
struct CartesianGrid
{
  unsigned nx = 0, ny = 0, nz = 0;
  double   x = 0, dx = 0, y = 0, dy = 0, z = 0, dz = 0;
};

// The values of a function on one plane of a grid, indexed by the point's
// x index i and y index j.
class PlaneValues
{
public:
  PlaneValues(unsigned nx, unsigned ny) : m_nx(nx), m_ny(ny), m_v(std::size_t(nx) * ny) {}

  double&  operator()(unsigned i, unsigned j)       { return m_v[std::size_t(i) * m_ny + j]; }
  double   operator()(unsigned i, unsigned j) const { return m_v[std::size_t(i) * m_ny + j]; }
  unsigned nx() const { return m_nx; }
  unsigned ny() const { return m_ny; }

private:
  unsigned            m_nx, m_ny;
  std::vector<double> m_v;
};

// The point of g at indices i, j, k, which may be fractional: g.x + i * g.dx
// and so on. The methods compute every position so, grid points and
// crossings alike.
inline Vec3 grid_point(const CartesianGrid& g, double i, double j, double k)
{
  return { g.x + g.dx * i, g.y + g.dy * j, g.z + g.dz * k };
}

// Fills values with the function on plane k of g, whose point (i, j) is
// grid_point(g, i, j, k).
using PlaneFunc = std::function<void(PlaneValues& values, const CartesianGrid& g, unsigned k)>;

// A PlaneFunc that evaluates f at each point of the plane.
PlaneFunc sample(std::function<double(Vec3)> f);

// The four methods add to s the faces of the surface f = iso, oriented
// towards f < iso: f > iso is the inside, as for a potential. For a
// function that grows outwards, such as a distance, negate f and iso.

// gts_isosurface_cartesian: marching cubes, slice by slice.
void isosurface_cartesian(Surface& s, CartesianGrid g, const PlaneFunc& f, double iso);

// gts_isosurface_tetra: each cube of the grid split into five tetrahedra.
void isosurface_tetra(Surface& s, CartesianGrid g, const PlaneFunc& f, double iso);

// gts_isosurface_tetra_bounded: isosurface_tetra on the grid extended by a
// layer of points with f - iso = -1, outside, on each side, a grid step
// beyond it, which closes the surface where it leaves the grid.
void isosurface_tetra_bounded(Surface& s, CartesianGrid g, const PlaneFunc& f, double iso);

// gts_isosurface_tetra_bcl: the body-centred variant, with tetrahedra
// joining the centre of each cube to the centres of its neighbours. It
// covers the grid between the centres of the outermost cubes.
void isosurface_tetra_bcl(Surface& s, CartesianGrid g, const PlaneFunc& f, double iso);

} // namespace gts
