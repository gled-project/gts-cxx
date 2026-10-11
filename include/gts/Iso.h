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

#include <array>
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

// GtsGridPlane: the points of one plane of a structured grid, which need not
// be cartesian, indexed as PlaneValues.
class GridPlane
{
public:
  GridPlane(unsigned nx, unsigned ny) : m_nx(nx), m_ny(ny), m_p(std::size_t(nx) * ny) {}

  Vec3&    operator()(unsigned i, unsigned j)       { return m_p[std::size_t(i) * m_ny + j]; }
  Vec3     operator()(unsigned i, unsigned j) const { return m_p[std::size_t(i) * m_ny + j]; }
  unsigned nx() const { return m_nx; }
  unsigned ny() const { return m_ny; }

private:
  unsigned          m_nx, m_ny;
  std::vector<Vec3> m_p;
};

// GtsIsoSlice: where the surface f = iso crosses the edges of one plane of
// a grid and the edges from it to the next plane. A crossing gets its
// vertex when a face first uses it. Marching cubes goes from slice to
// slice: fill one, make the faces between it and the one before with
// isosurface_slice(), and swap.
class IsoSlice
{
public:
  // Throws std::invalid_argument for fewer than two points in x or y.
  IsoSlice(unsigned nx, unsigned ny);

  // gts_iso_slice_fill: the crossings of the plane of points p1 and values
  // f1, and of the edges to the next plane, p2 and f2, which are null for
  // the last plane. The planes must have the slice's size.
  void fill(const GridPlane& p1, const GridPlane* p2, const PlaneValues& f1, const PlaneValues* f2, double iso);
  // gts_iso_slice_fill_cartesian: the same for plane k of g and plane k + 1.
  void fill_cartesian(const CartesianGrid& g, unsigned k, const PlaneValues& f1, const PlaneValues* f2, double iso);

  unsigned nx() const { return m_nx; }
  unsigned ny() const { return m_ny; }

private:
  friend void isosurface_slice(Surface& s, IsoSlice& s1, IsoSlice& s2);

  // A crossing; orientation tells whether f >= iso at the far end.
  struct Crossing
  {
    bool     exists = false;
    int      orientation = 0;
    Vec3     p;
    VertexId v;
  };
  // Arrays 0 and 1 have ny points per row, array 2 has ny - 1. Array 0
  // holds the edges to the next plane, 1 those along x, 2 those along y.
  Crossing& at(unsigned a, unsigned i, unsigned j) { return m_c[a][std::size_t(i) * (a == 2 ? m_ny - 1 : m_ny) + j]; }
  template <class Point>
  void fill_with(const Point& point, bool next, const PlaneValues& f1, const PlaneValues* f2, double iso);

  unsigned m_nx, m_ny;
  std::array<std::vector<Crossing>, 3> m_c;
};

// gts_isosurface_slice: adds to s the faces of the surface between the
// planes of s1 and s2, oriented towards f < iso.
void isosurface_slice(Surface& s, IsoSlice& s1, IsoSlice& s2);

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
