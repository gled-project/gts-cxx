// Copyright (C) 1999 Stéphane Popinet (GTS: iso.c)
// Copyright (C) 2002 Gert Wollny (GTS: isotetra.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/Iso.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace gts {

PlaneFunc sample(std::function<double(Vec3)> f)
{
  return [f = std::move(f)](PlaneValues& a, const CartesianGrid& g, unsigned k) {
    for (unsigned i = 0; i < g.nx; ++i)
      for (unsigned j = 0; j < g.ny; ++j)
        a(i, j) = f(grid_point(g, i, j, k));
  };
}

namespace {

void check_grid(const CartesianGrid& g, const char* who)
{
  if (g.nx < 2 || g.ny < 2 || g.nz < 2)
    throw std::invalid_argument(std::string("gts::") + who + ": the grid needs at least two points in each direction");
}

EdgeId edge_between(Geometrium& g, VertexId a, VertexId b)
{
  const EdgeId e = g.find_edge(a, b);
  return e.valid() ? e : g.add_edge(a, b);
}

//==============================================================================
// Marching cubes, from iso.c
//==============================================================================

bool crosses(double v1, double v2) { return (v1 >= 0. && v2 < 0.) || (v1 < 0. && v2 >= 0.); }

} // namespace

IsoSlice::IsoSlice(unsigned nx, unsigned ny) : m_nx(nx), m_ny(ny)
{
  if (nx < 2 || ny < 2)
    throw std::invalid_argument("gts::IsoSlice: a slice needs at least two points in x and y");
  m_c[0].resize(std::size_t(nx) * ny);
  m_c[1].resize(std::size_t(nx - 1) * ny);
  m_c[2].resize(std::size_t(nx) * (ny - 1));
}

// point(i, j, a, t): the point at t between point (i, j) and its neighbour
// along array a: the next plane for 0, i + 1 for 1, j + 1 for 2. GTS made
// each crossing's vertex at once; here it is made when a face first uses
// it, so that no vertex is left without edges.
template <class Point>
void IsoSlice::fill_with(const Point& point, bool next, const PlaneValues& f1, const PlaneValues* f2, double iso)
{
  if (next)
    for (unsigned i = 0; i < m_nx; i++)
      for (unsigned j = 0; j < m_ny; j++)
      {
        const double v1 = f1(i, j) - iso, v2 = (*f2)(i, j) - iso;
        Crossing& c = at(0, i, j);
        c = {};
        if (crosses(v1, v2))
          c = { true, v2 >= 0. ? 1 : 0, point(i, j, 0, v1 / (v1 - v2)), {} };
      }
  for (unsigned i = 0; i < m_nx - 1; i++)
    for (unsigned j = 0; j < m_ny; j++)
    {
      const double v1 = f1(i, j) - iso, v2 = f1(i + 1, j) - iso;
      Crossing& c = at(1, i, j);
      c = {};
      if (crosses(v1, v2))
        c = { true, v2 >= 0. ? 1 : 0, point(i, j, 1, v1 / (v1 - v2)), {} };
    }
  for (unsigned i = 0; i < m_nx; i++)
    for (unsigned j = 0; j < m_ny - 1; j++)
    {
      const double v1 = f1(i, j) - iso, v2 = f1(i, j + 1) - iso;
      Crossing& c = at(2, i, j);
      c = {};
      if (crosses(v1, v2))
        c = { true, v2 >= 0. ? 1 : 0, point(i, j, 2, v1 / (v1 - v2)), {} };
    }
}

void IsoSlice::fill(const GridPlane& p1, const GridPlane* p2, const PlaneValues& f1, const PlaneValues* f2, double iso)
{
  if (f2 && ! p2)
    throw std::invalid_argument("gts::IsoSlice::fill: values of the next plane without its points");
  auto point = [&](unsigned i, unsigned j, unsigned a, double c2) {
    const double c1 = 1. - c2;
    const Vec3 x = p1(i, j);
    const Vec3 y = a == 0 ? (*p2)(i, j) : a == 1 ? p1(i + 1, j) : p1(i, j + 1);
    return Vec3{ c1*x.x + c2*y.x, c1*x.y + c2*y.y, c1*x.z + c2*y.z };
  };
  fill_with(point, f2 != nullptr, f1, f2, iso);
}

void IsoSlice::fill_cartesian(const CartesianGrid& g, unsigned k, const PlaneValues& f1, const PlaneValues* f2, double iso)
{
  auto point = [&](unsigned i, unsigned j, unsigned a, double t) {
    return a == 0 ? grid_point(g, i, j, k + t) : a == 1 ? grid_point(g, i + t, j, k) : grid_point(g, i, j + t, k);
  };
  fill_with(point, f2 != nullptr, f1, f2, iso);
}

namespace {

// The edges of a cube: which array, which slice, and the offsets in x and y
// (see GTS's doc/isocube.fig).
constexpr unsigned kCube[12][4] = {
  {0, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 1}, {0, 0, 1, 0},
  {1, 0, 0, 0}, {1, 0, 0, 1}, {1, 1, 0, 1}, {1, 1, 0, 0},
  {2, 0, 0, 0}, {2, 1, 0, 0}, {2, 1, 1, 0}, {2, 0, 1, 0}};

// For each edge and orientation, the edges a crossing on it may connect to,
// in order.
constexpr unsigned kNext[12][2][3] = {
  {{9, 1, 8}, {4, 3, 7}},   {{6, 2, 5}, {8, 0, 9}},   {{10, 3, 11}, {5, 1, 6}}, {{7, 0, 4}, {11, 2, 10}},
  {{3, 7, 0}, {8, 5, 11}},  {{11, 4, 8}, {1, 6, 2}},  {{2, 5, 1}, {9, 7, 10}},  {{10, 6, 9}, {0, 4, 3}},
  {{5, 11, 4}, {0, 9, 1}},  {{1, 8, 0}, {7, 10, 6}},  {{6, 9, 7}, {3, 11, 2}},  {{2, 10, 3}, {4, 8, 5}}};

// The faces of one cube, whose crossing on edge e is at(e). The crossings
// are walked round in loops, and each loop of more than two is fanned into
// faces from its first crossing.
template <class At>
void cube_faces(Surface& s, At&& at)
{
  Geometrium& g = s.geometrium();
  // GTS marked the visited vertices with their reserved pointer.
  bool visited[12] = {};
  for (unsigned l = 0; l < 12; l++)
  {
    std::remove_reference_t<decltype(at(0))>* va[12];
    unsigned nv = 0, e = l;
    auto* ov = &at(e);
    while (ov->exists && ! visited[e])
    {
      va[nv++] = ov;
      visited[e] = true;
      const unsigned* ne = kNext[e][ov->orientation];
      unsigned m = 0;
      do { e = ne[m++]; ov = &at(e); } while (m < 3 && ! ov->exists);
    }
    if (nv > 2)
    {
      for (unsigned m = 0; m < nv; ++m)
        if ( ! va[m]->v.valid()) va[m]->v = g.add_vertex(va[m]->p);
      EdgeId e1 = edge_between(g, va[0]->v, va[1]->v);
      for (unsigned m = 1; m < nv - 1; m++)
      {
        const EdgeId e2 = edge_between(g, va[m]->v, va[m + 1]->v);
        const EdgeId e3 = edge_between(g, va[m + 1]->v, va[0]->v);
        // Towards f < iso: GTS's order, (e1, e2, e3), reverted.
        s.add_face(g.add_face(e2, e1, e3));
        e1 = e3;
      }
    }
  }
}

} // namespace

void isosurface_slice(Surface& s, IsoSlice& s1, IsoSlice& s2)
{
  if (s1.m_nx != s2.m_nx || s1.m_ny != s2.m_ny)
    throw std::invalid_argument("gts::isosurface_slice: slices of different sizes");
  IsoSlice* slices[2] = { &s1, &s2 };
  for (unsigned j = 0; j + 1 < s1.m_nx; j++)
    for (unsigned k = 0; k + 1 < s1.m_ny; k++)
      cube_faces(s, [&](unsigned e) -> IsoSlice::Crossing& {
        return slices[kCube[e][1]]->at(kCube[e][0], j + kCube[e][2], k + kCube[e][3]);
      });
}

namespace {

//==============================================================================
// Tetrahedra, from isotetra.c
//==============================================================================

// A corner of a tetrahedron: grid indices, whether it is a cube centre
// (body-centred variant only), and f - iso there.
struct TetraVertex
{
  int    x, y, z;
  bool   mid;
  double d;
};

using Slice = PlaneValues;

void iso_sub(Slice& s, double iso)
{
  for (unsigned x = 0; x < s.nx(); ++x)
    for (unsigned y = 0; y < s.ny(); ++y)
      s(x, y) -= iso;
}

void slice_init(Slice& s, double v)
{
  for (unsigned x = 0; x < s.nx(); ++x)
    for (unsigned y = 0; y < s.ny(); ++y)
      s(x, y) = v;
}

// helper_t: the vertices made so far, by edge of the grid, for the planes
// at the bottom (vtop, named so in GTS) and the top (vbot) of the current
// layer of cubes, and for the edges between them (vmid).
// off is the index of the grid's first point, 1 in the enlarged grid of the
// bounded method.
struct Helper
{
  int nx, ny, off;
  std::vector<VertexId> vtop, vmid, vbot;

  Helper(int nx_, int ny_, int off_ = 0) :
    nx(nx_), ny(ny_), off(off_), vtop(4 * nx * ny), vmid(4 * nx * ny), vbot(4 * nx * ny) {}

  // helper_advance: the top plane becomes the bottom one.
  void advance()
  {
    std::swap(vbot, vtop);
    std::fill(vmid.begin(), vmid.end(), VertexId());
    std::fill(vbot.begin(), vbot.end(), VertexId());
  }
};

// get_vertex: the crossing on the line from v1 to v2, made once.
VertexId get_vertex(Geometrium& geo, int mz, const TetraVertex& v1, const TetraVertex& v2, Helper& h, const CartesianGrid& g)
{
  if (v1.d - v2.d == 0.)
    throw std::logic_error("gts::isosurface_tetra: no crossing between two corners");
  double dx = 0.0, dy = 0.0, dz = 0.0;
  const double d = v1.d / (v1.d - v2.d);
  int index = 0;
  if (v1.x != v2.x) { index |= 1; dx = d; }
  if (v1.y != v2.y) { index |= 2; dy = d; }
  if (v1.z != v2.z) { dz = d; }

  int x = v1.x;
  if (v1.x > v2.x) { x = v2.x; dx = 1.0 - dx; }
  int y = v1.y;
  if (v1.y > v2.y) { y = v2.y; dy = 1.0 - dy; }
  int z = v1.z;
  if (v1.z > v2.z) { z = v2.z; dz = 1.0 - dz; }

  const std::size_t idx2 = std::size_t(4 * (x + y * h.nx) + index);
  VertexId& v = v1.z == v2.z ? (mz == z ? h.vtop[idx2] : h.vbot[idx2]) : h.vmid[idx2];
  if ( ! v.valid())
    v = geo.add_vertex(grid_point(g, x - h.off + dx, y - h.off + dy, z - h.off + dz));
  return v;
}

// add_face: the face through the crossings on a1-a2, b1-b2 and c1-c2.
void add_tetra_face(Surface& s, const TetraVertex& a1, const TetraVertex& a2, const TetraVertex& b1,
                    const TetraVertex& b2, const TetraVertex& c1, const TetraVertex& c2, bool rev, Helper& h,
                    int z, const CartesianGrid& g)
{
  Geometrium& geo = s.geometrium();
  const VertexId v1 = get_vertex(geo, z, a1, a2, h, g);
  const VertexId v2 = get_vertex(geo, z, b1, b2, h, g);
  const VertexId v3 = get_vertex(geo, z, c1, c2, h, g);
  if (v1 == v2 || v2 == v3 || v1 == v3)
    throw std::logic_error("gts::isosurface_tetra: a face with two crossings at one vertex");
  EdgeId e1, e2, e3;
  if ( ! rev) { e1 = edge_between(geo, v1, v2); e2 = edge_between(geo, v2, v3); e3 = edge_between(geo, v1, v3); }
  else        { e1 = edge_between(geo, v1, v3); e2 = edge_between(geo, v2, v3); e3 = edge_between(geo, v1, v2); }
  s.add_face(geo.add_face(e1, e2, e3));
}

// analyze_tetrahedra: the faces in one tetrahedron, by which corners have
// d >= 0. The cases of the complement have the opposite orientation.
template <class AddFace>
void analyze_tetrahedra(const TetraVertex& a, const TetraVertex& b, const TetraVertex& c, const TetraVertex& d,
                        bool parity, AddFace&& add)
{
  bool rev = parity;
  int code = 0;
  if (a.d >= 0.) code |= 1;
  if (b.d >= 0.) code |= 2;
  if (c.d >= 0.) code |= 4;
  if (d.d >= 0.) code |= 8;

  switch (code)
  {
    case 15:
    case 0: return;
    case 14: rev = ! parity; [[fallthrough]];
    case 1:  add(a, b, a, d, a, c, rev); break;
    case 13: rev = ! parity; [[fallthrough]];
    case 2:  add(a, b, b, c, b, d, rev); break;
    case 12: rev = ! parity; [[fallthrough]];
    case 3:  add(a, d, a, c, b, c, rev); add(a, d, b, c, b, d, rev); break;
    case 11: rev = ! parity; [[fallthrough]];
    case 4:  add(a, c, c, d, b, c, rev); break;
    case 10: rev = ! parity; [[fallthrough]];
    case 5:  add(a, b, a, d, c, d, rev); add(a, b, c, d, b, c, rev); break;
    case 9:  rev = ! parity; [[fallthrough]];
    case 6:  add(a, b, a, c, c, d, rev); add(a, b, c, d, b, d, rev); break;
    case 7:  rev = ! parity; [[fallthrough]];
    case 8:  add(a, d, b, d, c, d, rev); break;
  }
}

// iso_slice_evaluate: each cube between two slices, split into five
// tetrahedra, in one of two ways by the parity of the cube.
void iso_slice_evaluate(const Slice& s1, const Slice& s2, const CartesianGrid& g, int z, Surface& s, Helper& h)
{
  for (int y = 0; y < int(g.ny) - 1; y++)
    for (int x = 0; x < int(g.nx) - 1; x++)
    {
      const bool parity = ((x ^ y) ^ z) & 1;
      const TetraVertex v0 = { x,     y,     z,     false, s1(x,     y    ) };
      const TetraVertex v1 = { x,     y + 1, z,     false, s1(x,     y + 1) };
      const TetraVertex v2 = { x + 1, y,     z,     false, s1(x + 1, y    ) };
      const TetraVertex v3 = { x + 1, y + 1, z,     false, s1(x + 1, y + 1) };
      const TetraVertex v4 = { x,     y,     z + 1, false, s2(x,     y    ) };
      const TetraVertex v5 = { x,     y + 1, z + 1, false, s2(x,     y + 1) };
      const TetraVertex v6 = { x + 1, y,     z + 1, false, s2(x + 1, y    ) };
      const TetraVertex v7 = { x + 1, y + 1, z + 1, false, s2(x + 1, y + 1) };
      auto add = [&](const TetraVertex& a1, const TetraVertex& a2, const TetraVertex& b1, const TetraVertex& b2,
                     const TetraVertex& c1, const TetraVertex& c2, bool rev) {
        add_tetra_face(s, a1, a2, b1, b2, c1, c2, rev, h, z, g);
      };
      if ( ! parity)
      {
        analyze_tetrahedra(v0, v1, v2, v4, parity, add);
        analyze_tetrahedra(v7, v1, v4, v2, parity, add);
        analyze_tetrahedra(v1, v7, v3, v2, parity, add);
        analyze_tetrahedra(v1, v7, v4, v5, parity, add);
        analyze_tetrahedra(v2, v6, v4, v7, parity, add);
      }
      else
      {
        analyze_tetrahedra(v4, v5, v6, v0, parity, add);
        analyze_tetrahedra(v3, v5, v0, v6, parity, add);
        analyze_tetrahedra(v5, v3, v7, v6, parity, add);
        analyze_tetrahedra(v5, v3, v0, v1, parity, add);
        analyze_tetrahedra(v6, v2, v0, v3, parity, add);
      }
    }
}

//==============================================================================
// Body-centred tetrahedra, from isotetra.c
//==============================================================================

// helper_bcl: GTS keyed the vertices by strings of the corners' indices,
// "x y z mid", in hash tables. Here the key is the corner, or the two
// corners in a fixed order; any fixed order names a pair the same way.
using CornerKey = std::array<int, 4>;
using BclKey    = std::pair<CornerKey, CornerKey>; // the second is {} for a single corner

struct HelperBcl
{
  std::map<BclKey, VertexId> vtop, vbot;
  void advance() { vbot = std::move(vtop); vtop.clear(); }
};

VertexId get_vertex_bcl(Geometrium& geo, int mz, const TetraVertex& v1, const TetraVertex& v2, HelperBcl& h,
                        const CartesianGrid& g)
{
  if (v1.d - v2.d == 0.)
    throw std::logic_error("gts::isosurface_tetra_bcl: no crossing between two corners");
  auto& table = v1.z > mz && v2.z > mz ? h.vtop : h.vbot;
  const double d = v1.d / (v1.d - v2.d);

  const CornerKey k1 = { v1.x, v1.y, v1.z, v1.mid }, k2 = { v2.x, v2.y, v2.z, v2.mid };
  const BclKey key = d == 0.0 ? BclKey(k1, {}) : d == 1.0 ? BclKey(k2, {}) :
                     k1 < k2 ? BclKey(k1, k2) : BclKey(k2, k1);
  auto [it, made] = table.try_emplace(key);
  if (made)
  {
    // Between the positions of the two corners; a cube centre is half a
    // step from its corner.
    const Vec3 p1 = grid_point(g, v1.x + v1.mid / 2.0, v1.y + v1.mid / 2.0, v1.z + v1.mid / 2.0);
    const Vec3 p2 = grid_point(g, v2.x + v2.mid / 2.0, v2.y + v2.mid / 2.0, v2.z + v2.mid / 2.0);
    it->second = geo.add_vertex({ p1.x * (1.0 - d) + d * p2.x, p1.y * (1.0 - d) + d * p2.y, p1.z * (1.0 - d) + d * p2.z });
  }
  return it->second;
}

void add_face_bcl(Surface& s, const TetraVertex& a1, const TetraVertex& a2, const TetraVertex& b1,
                  const TetraVertex& b2, const TetraVertex& c1, const TetraVertex& c2, bool rev, HelperBcl& h,
                  int z, const CartesianGrid& g)
{
  Geometrium& geo = s.geometrium();
  const VertexId v1 = get_vertex_bcl(geo, z, a1, a2, h, g);
  const VertexId v2 = get_vertex_bcl(geo, z, b1, b2, h, g);
  const VertexId v3 = get_vertex_bcl(geo, z, c1, c2, h, g);
  if (v1 == v2 || v2 == v3 || v1 == v3)
    return;
  EdgeId e1, e2, e3;
  if ( ! rev) { e1 = edge_between(geo, v1, v2); e2 = edge_between(geo, v2, v3); e3 = edge_between(geo, v1, v3); }
  else        { e1 = edge_between(geo, v1, v3); e2 = edge_between(geo, v2, v3); e3 = edge_between(geo, v1, v2); }
  s.add_face(geo.add_face(e1, e2, e3));
}

// iso_slice_evaluate_bcl for the cubes between s1 and s2: each cube's
// centre is joined to the centres of the next cubes in x, y and z, where
// they exist; s3, the slice above s2, is null for the last layer.
void iso_slice_evaluate_bcl(const Slice& s1, const Slice& s2, const Slice* s3, const CartesianGrid& g, int z,
                            Surface& s, HelperBcl& h)
{
  auto centre = [](const Slice& a, const Slice& b, int x, int y) {
    return (a(x, y) + b(x, y) + a(x, y + 1) + b(x, y + 1) +
            a(x + 1, y) + b(x + 1, y) + a(x + 1, y + 1) + b(x + 1, y + 1)) / 8.0;
  };
  const int nx = int(g.nx) - 1, ny = int(g.ny) - 1; // cubes in x and y
  for (int y = 0; y < ny; y++)
    for (int x = 0; x < nx; x++)
    {
      const TetraVertex v0 = { x, y, z, true, centre(s1, s2, x, y) };
      const TetraVertex v4 = { x + 1, y,     z,     false, s1(x + 1, y    ) };
      const TetraVertex v5 = { x,     y + 1, z,     false, s1(x,     y + 1) };
      const TetraVertex v6 = { x + 1, y + 1, z,     false, s1(x + 1, y + 1) };
      const TetraVertex v7 = { x + 1, y,     z + 1, false, s2(x + 1, y    ) };
      const TetraVertex v8 = { x,     y + 1, z + 1, false, s2(x,     y + 1) };
      const TetraVertex v9 = { x + 1, y + 1, z + 1, false, s2(x + 1, y + 1) };
      const TetraVertex w0 = { x,     y,     z + 1, false, s2(x,     y    ) };
      auto add = [&](const TetraVertex& a1, const TetraVertex& a2, const TetraVertex& b1, const TetraVertex& b2,
                     const TetraVertex& c1, const TetraVertex& c2, bool rev) {
        add_face_bcl(s, a1, a2, b1, b2, c1, c2, rev, h, z, g);
      };
      // analyze_tetrahedra_bcl starts from rev = 0 instead of the parity.
      if (x + 1 < nx)
      {
        const TetraVertex v1 = { x + 1, y, z, true, centre(s1, s2, x + 1, y) };
        analyze_tetrahedra(v0, v9, v6, v1, false, add);
        analyze_tetrahedra(v0, v6, v4, v1, false, add);
        analyze_tetrahedra(v0, v4, v7, v1, false, add);
        analyze_tetrahedra(v0, v7, v9, v1, false, add);
      }
      if (y + 1 < ny)
      {
        const TetraVertex v2 = { x, y + 1, z, true, centre(s1, s2, x, y + 1) };
        analyze_tetrahedra(v0, v5, v6, v2, false, add);
        analyze_tetrahedra(v0, v6, v9, v2, false, add);
        analyze_tetrahedra(v0, v9, v8, v2, false, add);
        analyze_tetrahedra(v0, v8, v5, v2, false, add);
      }
      if (s3)
      {
        const TetraVertex v3 = { x, y, z + 1, true, centre(s2, *s3, x, y) };
        analyze_tetrahedra(v0, v8, v9, v3, false, add);
        analyze_tetrahedra(v0, v9, v7, v3, false, add);
        analyze_tetrahedra(v0, v7, w0, v3, false, add);
        analyze_tetrahedra(v0, w0, v8, v3, false, add);
      }
    }
}

// copy_to_bounded: src minus iso inside a border of fill.
void copy_to_bounded(Slice& dest, const Slice& src, double iso, double fill)
{
  slice_init(dest, fill);
  for (unsigned x = 0; x < src.nx(); ++x)
    for (unsigned y = 0; y < src.ny(); ++y)
      dest(x + 1, y + 1) = src(x, y) - iso;
}

} // namespace

//==============================================================================
// The four methods
//==============================================================================

void isosurface_cartesian(Surface& s, CartesianGrid g, const PlaneFunc& f, double iso)
{
  check_grid(g, "isosurface_cartesian");
  IsoSlice slice_a(g.nx, g.ny), slice_b(g.nx, g.ny);
  IsoSlice* slice1 = &slice_a;
  IsoSlice* slice2 = &slice_b;
  PlaneValues fa(g.nx, g.ny), fb(g.nx, g.ny);
  PlaneValues* f1 = &fa;
  PlaneValues* f2 = &fb;

  f(*f1, g, 0);
  f(*f2, g, 1);
  slice1->fill_cartesian(g, 0, *f1, f2, iso);
  for (unsigned k = 2; k < g.nz; k++)
  {
    f(*f1, g, k);
    std::swap(f1, f2);
    slice2->fill_cartesian(g, k - 1, *f1, f2, iso);
    isosurface_slice(s, *slice1, *slice2);
    std::swap(slice1, slice2);
  }
  slice2->fill_cartesian(g, g.nz - 1, *f2, nullptr, iso);
  isosurface_slice(s, *slice1, *slice2);
}

//==============================================================================
// Marching cubes from seeds
//==============================================================================

std::vector<GridCell> cells_near(const CartesianGrid& g, std::span<const Vec3> points, unsigned band)
{
  check_grid(g, "cells_near");
  // The cells from c - band to c + band along one axis, clipped to the grid.
  auto range = [band](double x, double x0, double dx, unsigned n, long& lo, long& hi) {
    const double c = std::floor((x - x0) / dx);
    if ( ! std::isfinite(c)) return false;
    const double last = n - 2;
    if (c + band < 0. || c - band > last) return false;
    lo = long(std::max(c - band, 0.));
    hi = long(std::min(c + band, last));
    return true;
  };
  std::vector<GridCell> cells;
  for (const Vec3& p : points)
  {
    long i0, i1, j0, j1, k0, k1;
    if ( ! range(p.x, g.x, g.dx, g.nx, i0, i1) || ! range(p.y, g.y, g.dy, g.ny, j0, j1) ||
         ! range(p.z, g.z, g.dz, g.nz, k0, k1))
      continue;
    for (long i = i0; i <= i1; ++i)
      for (long j = j0; j <= j1; ++j)
        for (long k = k0; k <= k1; ++k)
          cells.push_back({ unsigned(i), unsigned(j), unsigned(k) });
  }
  std::ranges::sort(cells);
  cells.erase(std::ranges::unique(cells).begin(), cells.end());
  return cells;
}

std::vector<GridCell> cells_near(const CartesianGrid& g, const Surface& s, unsigned band)
{
  std::vector<Vec3> points;
  for (VertexId v : s.vertices()) points.push_back(s.geometrium().vertex(v).p);
  return cells_near(g, points, band);
}

std::vector<GridCell> cells_crossing(const CartesianGrid& g, const std::function<double(Vec3)>& f, double iso,
                                     Vec3 a, Vec3 b)
{
  double fa = f(a) - iso;
  if ( ! crosses(fa, f(b) - iso)) return {};
  const double step = std::min({ g.dx, g.dy, g.dz }) / 4;
  for (int n = 0; n < 100 && distance(a, b) > step; ++n)
  {
    const Vec3 m = 0.5 * (a + b);
    const double fm = f(m) - iso;
    if (crosses(fa, fm)) b = m;
    else                 { a = m; fa = fm; }
  }
  const Vec3 m = 0.5 * (a + b);
  return cells_near(g, std::span(&m, 1), 1);
}

// The cells are found first, by a flood fill through the cube faces with
// crossings, and then made in the order of isosurface_cartesian().
void isosurface_cartesian(Surface& s, CartesianGrid g, const std::function<double(Vec3)>& f, double iso,
                          std::span<const GridCell> seeds)
{
  check_grid(g, "isosurface_cartesian");
  const std::uint64_t nx = g.nx, ny = g.ny;
  auto point_id = [&](unsigned i, unsigned j, unsigned k) { return (k * ny + j) * nx + i; };

  // f - iso at the points, each evaluated once.
  std::unordered_map<std::uint64_t, double> values;
  auto value = [&](unsigned i, unsigned j, unsigned k) {
    auto [it, made] = values.try_emplace(point_id(i, j, k));
    if (made) it->second = f(grid_point(g, i, j, k)) - iso;
    return it->second;
  };

  std::vector<GridCell> cells;
  std::unordered_set<std::uint64_t> seen;
  std::vector<GridCell> queue;
  auto push = [&](GridCell c) {
    if (c.i < g.nx - 1 && c.j < g.ny - 1 && c.k < g.nz - 1 && seen.insert(point_id(c.i, c.j, c.k)).second)
      queue.push_back(c);
  };
  for (const GridCell& c : seeds) push(c);
  while ( ! queue.empty())
  {
    const GridCell c = queue.back();
    queue.pop_back();
    // Corner n: bit 0 is its offset in x, bit 1 in y, bit 2 in z.
    bool above[8];
    for (unsigned n = 0; n < 8; ++n)
      above[n] = value(c.i + (n & 1), c.j + (n >> 1 & 1), c.k + (n >> 2)) >= 0.;
    if (std::ranges::all_of(above, [&](bool a) { return a == above[0]; })) continue;
    cells.push_back(c);
    // Each face of the cube: the corners with bit a equal to side.
    for (unsigned a = 0; a < 3; ++a)
      for (unsigned side = 0; side < 2; ++side)
      {
        bool any_above = false, any_below = false;
        for (unsigned n = 0; n < 8; ++n)
          if ((n >> a & 1) == side) (above[n] ? any_above : any_below) = true;
        if ( ! (any_above && any_below)) continue;
        GridCell d = c;
        unsigned& x = a == 0 ? d.i : a == 1 ? d.j : d.k;
        if (side == 0 && x == 0) continue;
        x = side ? x + 1 : x - 1;
        push(d);
      }
  }
  std::ranges::sort(cells, {}, [](const GridCell& c) { return std::tuple(c.k, c.i, c.j); });

  // The crossings by grid edge: the edge's lower point and kCube's array,
  // 0 along z, 1 along x, 2 along y; computed as fill_cartesian() does.
  struct Crossing
  {
    bool     exists = false;
    int      orientation = 0;
    Vec3     p;
    VertexId v;
  };
  std::unordered_map<std::uint64_t, Crossing> crossings;
  auto crossing = [&](unsigned a, unsigned i, unsigned j, unsigned k) -> Crossing& {
    auto [it, made] = crossings.try_emplace(point_id(i, j, k) * 3 + a);
    if (made)
    {
      const double v1 = value(i, j, k);
      const double v2 = a == 0 ? value(i, j, k + 1) : a == 1 ? value(i + 1, j, k) : value(i, j + 1, k);
      if (crosses(v1, v2))
      {
        const double t = v1 / (v1 - v2);
        it->second = { true, v2 >= 0. ? 1 : 0,
                       a == 0 ? grid_point(g, i, j, k + t) : a == 1 ? grid_point(g, i + t, j, k) : grid_point(g, i, j + t, k),
                       {} };
      }
    }
    return it->second;
  };
  for (const GridCell& c : cells)
    cube_faces(s, [&](unsigned e) -> Crossing& {
      return crossing(kCube[e][0], c.i + kCube[e][2], c.j + kCube[e][3], c.k + kCube[e][1]);
    });
}

void isosurface_tetra(Surface& s, CartesianGrid g, const PlaneFunc& f, double iso)
{
  check_grid(g, "isosurface_tetra");
  Slice sa(g.nx, g.ny), sb(g.nx, g.ny);
  Slice* slice1 = &sa;
  Slice* slice2 = &sb;
  Helper helper(int(g.nx), int(g.ny));

  f(*slice1, g, 0);
  iso_sub(*slice1, iso);
  for (unsigned z = 1; z < g.nz; ++z)
  {
    f(*slice2, g, z);
    iso_sub(*slice2, iso);
    iso_slice_evaluate(*slice1, *slice2, g, int(z) - 1, s, helper);
    std::swap(slice1, slice2);
    helper.advance();
  }
}

void isosurface_tetra_bounded(Surface& s, CartesianGrid g, const PlaneFunc& f, double iso)
{
  check_grid(g, "isosurface_tetra_bounded");
  Slice sa(g.nx + 2, g.ny + 2), sb(g.nx + 2, g.ny + 2);
  Slice* slice1 = &sa;
  Slice* slice2 = &sb;
  slice_init(*slice1, -1.0);
  Slice transfer(g.nx, g.ny);

  // The enlarged grid: index 0 is the layer, a step before g's first point.
  CartesianGrid gi = g;
  gi.nx = g.nx + 2;
  gi.ny = g.ny + 2;
  Helper helper(int(gi.nx), int(gi.ny), 1);

  unsigned z = 0;
  for (; z < g.nz; ++z)
  {
    f(transfer, g, z);
    copy_to_bounded(*slice2, transfer, iso, -1.);
    iso_slice_evaluate(*slice1, *slice2, gi, int(z), s, helper);
    std::swap(slice1, slice2);
    helper.advance();
  }
  slice_init(*slice2, -1.0);
  iso_slice_evaluate(*slice1, *slice2, gi, int(z), s, helper);
}

void isosurface_tetra_bcl(Surface& s, CartesianGrid g, const PlaneFunc& f, double iso)
{
  check_grid(g, "isosurface_tetra_bcl");
  Slice sa(g.nx, g.ny), sb(g.nx, g.ny), sc(g.nx, g.ny);
  Slice* slice1 = &sa;
  Slice* slice2 = &sb;
  Slice* slice3 = &sc;
  HelperBcl helper;

  f(*slice1, g, 0);
  iso_sub(*slice1, iso);
  f(*slice2, g, 1);
  iso_sub(*slice2, iso);
  unsigned z = 2;
  for (; z < g.nz; ++z)
  {
    f(*slice3, g, z);
    iso_sub(*slice3, iso);
    iso_slice_evaluate_bcl(*slice1, *slice2, slice3, g, int(z) - 2, s, helper);
    Slice* hs = slice1;
    slice1 = slice2;
    slice2 = slice3;
    slice3 = hs;
    helper.advance();
  }
  iso_slice_evaluate_bcl(*slice1, *slice2, nullptr, g, int(z) - 2, s, helper);
}

} // namespace gts
