// Copyright (C) 1999 Stéphane Popinet (GTS: iso.c)
// Copyright (C) 2002 Gert Wollny (GTS: isotetra.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/Iso.h"

#include <array>
#include <map>
#include <stdexcept>
#include <utility>

namespace gts {

PlaneFunc sample(std::function<double(Vec3)> f)
{
  return [f = std::move(f)](PlaneValues& a, const CartesianGrid& g, unsigned) {
    double x = g.x;
    for (unsigned i = 0; i < g.nx; ++i, x += g.dx)
    {
      double y = g.y;
      for (unsigned j = 0; j < g.ny; ++j, y += g.dy)
        a(i, j) = f({ x, y, g.z });
    }
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

// A crossing of the surface on an edge of the grid. GTS made the vertex of
// every crossing at once; here it is made when a face first uses it, so
// that no vertex is left without edges.
struct Crossing
{
  bool     exists = false;
  int      orientation = 0; // LEFT 0, RIGHT 1: whether f >= iso at the far end
  Vec3     p;
  VertexId v;
};

// GtsIsoSlice: the crossings on the z edges from one plane to the next
// (array 0), and on the x and y edges in the plane (arrays 1 and 2).
struct IsoSlice
{
  unsigned nx, ny;
  std::array<std::vector<Crossing>, 3> c;

  IsoSlice(unsigned nx_, unsigned ny_) : nx(nx_), ny(ny_)
  {
    c[0].resize(std::size_t(nx) * ny);
    c[1].resize(std::size_t(nx - 1) * ny);
    c[2].resize(std::size_t(nx) * (ny - 1));
  }
  // Arrays 0 and 1 have ny points per row, array 2 has ny - 1.
  Crossing& at(unsigned a, unsigned i, unsigned j) { return c[a][std::size_t(i) * (a == 2 ? ny - 1 : ny) + j]; }
};

bool crosses(double v1, double v2) { return (v1 >= 0. && v2 < 0.) || (v1 < 0. && v2 >= 0.); }

// gts_iso_slice_fill_cartesian, with f2 null for the last plane.
void slice_fill(IsoSlice& s, const CartesianGrid& g, const PlaneValues& f1, const PlaneValues* f2, double iso)
{
  unsigned i, j;
  double x, y;
  if (f2)
    for (i = 0, x = g.x; i < g.nx; i++, x += g.dx)
      for (j = 0, y = g.y; j < g.ny; j++, y += g.dy)
      {
        const double v1 = f1(i, j) - iso, v2 = (*f2)(i, j) - iso;
        Crossing& c = s.at(0, i, j);
        c = {};
        if (crosses(v1, v2))
          c = { true, v2 >= 0. ? 1 : 0, { x, y, g.z + g.dz * v1 / (v1 - v2) }, {} };
      }
  for (i = 0, x = g.x; i < g.nx - 1; i++, x += g.dx)
    for (j = 0, y = g.y; j < g.ny; j++, y += g.dy)
    {
      const double v1 = f1(i, j) - iso, v2 = f1(i + 1, j) - iso;
      Crossing& c = s.at(1, i, j);
      c = {};
      if (crosses(v1, v2))
        c = { true, v2 >= 0. ? 1 : 0, { x + g.dx * v1 / (v1 - v2), y, g.z }, {} };
    }
  for (i = 0, x = g.x; i < g.nx; i++, x += g.dx)
    for (j = 0, y = g.y; j < g.ny - 1; j++, y += g.dy)
    {
      const double v1 = f1(i, j) - iso, v2 = f1(i, j + 1) - iso;
      Crossing& c = s.at(2, i, j);
      c = {};
      if (crosses(v1, v2))
        c = { true, v2 >= 0. ? 1 : 0, { x, y + g.dy * v1 / (v1 - v2), g.z }, {} };
    }
}

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

// gts_isosurface_slice: the faces between two slices. In each cube, the
// crossings are walked round in loops, and each loop of more than two is
// fanned into faces from its first crossing.
void isosurface_slice(Surface& s, IsoSlice& s1, IsoSlice& s2)
{
  Geometrium& g = s.geometrium();
  IsoSlice* slices[2] = { &s1, &s2 };
  for (unsigned j = 0; j + 1 < s1.nx; j++)
    for (unsigned k = 0; k + 1 < s1.ny; k++)
    {
      auto at = [&](unsigned e) -> Crossing& {
        return slices[kCube[e][1]]->at(kCube[e][0], j + kCube[e][2], k + kCube[e][3]);
      };
      // GTS marked the visited vertices with their reserved pointer.
      bool visited[12] = {};
      for (unsigned l = 0; l < 12; l++)
      {
        Crossing* va[12];
        unsigned nv = 0, e = l;
        Crossing* ov = &at(e);
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
            s.add_face(g.add_face(e1, e2, e3));
            e1 = e3;
          }
        }
      }
    }
}

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
struct Helper
{
  int nx, ny;
  std::vector<VertexId> vtop, vmid, vbot;

  Helper(int nx_, int ny_) : nx(nx_), ny(ny_), vtop(4 * nx * ny), vmid(4 * nx * ny), vbot(4 * nx * ny) {}

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
    v = geo.add_vertex({ g.dx * (x + dx) + g.x, g.dy * (y + dy) + g.y, g.dz * (z + dz) + g.z });
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
    const double x1 = g.dx * (v1.x + (v1.mid / 2.0)) + g.x;
    const double x2 = g.dx * (v2.x + (v2.mid / 2.0)) + g.x;
    const double y1 = g.dy * (v1.y + (v1.mid / 2.0)) + g.y;
    const double y2 = g.dy * (v2.y + (v2.mid / 2.0)) + g.y;
    const double z1 = g.dz * (v1.z + (v1.mid / 2.0)) + g.z;
    const double z2 = g.dz * (v2.z + (v2.mid / 2.0)) + g.z;
    it->second = geo.add_vertex({ x1 * (1.0 - d) + d * x2, y1 * (1.0 - d) + d * y2, z1 * (1.0 - d) + d * z2 });
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

void iso_slice_evaluate_bcl(const Slice& s1, const Slice& s2, const Slice& s3, const CartesianGrid& g, int z,
                            Surface& s, HelperBcl& h)
{
  for (int y = 0; y < int(g.ny) - 2; y++)
    for (int x = 0; x < int(g.nx) - 2; x++)
    {
      const TetraVertex v0 = { x, y, z, true,
        (s1(x, y) + s2(x, y) + s1(x, y + 1) + s2(x, y + 1) +
         s1(x + 1, y) + s2(x + 1, y) + s1(x + 1, y + 1) + s2(x + 1, y + 1)) / 8.0 };
      const TetraVertex v1 = { x + 1, y, z, true,
        (s1(x + 1, y) + s2(x + 1, y) + s1(x + 1, y + 1) + s2(x + 1, y + 1) +
         s1(x + 2, y) + s2(x + 2, y) + s1(x + 2, y + 1) + s2(x + 2, y + 1)) / 8.0 };
      const TetraVertex v2 = { x, y + 1, z, true,
        (s1(x, y + 1) + s2(x, y + 1) + s1(x, y + 2) + s2(x, y + 2) +
         s1(x + 1, y + 1) + s2(x + 1, y + 1) + s1(x + 1, y + 2) + s2(x + 1, y + 2)) / 8.0 };
      const TetraVertex v3 = { x, y, z + 1, true,
        (s2(x, y) + s3(x, y) + s2(x, y + 1) + s3(x, y + 1) +
         s2(x + 1, y) + s3(x + 1, y) + s2(x + 1, y + 1) + s3(x + 1, y + 1)) / 8.0 };
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
      analyze_tetrahedra(v0, v9, v6, v1, false, add);
      analyze_tetrahedra(v0, v6, v4, v1, false, add);
      analyze_tetrahedra(v0, v4, v7, v1, false, add);
      analyze_tetrahedra(v0, v7, v9, v1, false, add);
      analyze_tetrahedra(v0, v5, v6, v2, false, add);
      analyze_tetrahedra(v0, v6, v9, v2, false, add);
      analyze_tetrahedra(v0, v9, v8, v2, false, add);
      analyze_tetrahedra(v0, v8, v5, v2, false, add);
      analyze_tetrahedra(v0, v8, v9, v3, false, add);
      analyze_tetrahedra(v0, v9, v7, v3, false, add);
      analyze_tetrahedra(v0, v7, w0, v3, false, add);
      analyze_tetrahedra(v0, w0, v8, v3, false, add);
    }
}

// copy_to_bounded, with GTS's two indexing errors kept; see QUIRKS.md.
void copy_to_bounded(Slice& dest, const Slice& src, double iso, double fill)
{
  for (unsigned y = 0; y < dest.ny(); ++y) dest(0, y) = fill;
  unsigned y = dest.ny();
  for (unsigned x = 1; x + 1 < src.nx(); ++x)
  {
    dest(x, 0) = fill;
    for (y = 0; y < src.ny(); ++y) dest(x, y + 1) = src(x - 1, y) - iso;
    dest(x, y + 1) = fill;
  }
  if (y < dest.nx())
    for (unsigned k = 0; k < dest.ny(); ++k) dest(y, k) = fill;
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

  // The steps of g.z are GTS's, so that z rounds as in GTS.
  f(*f1, g, 0);
  g.z += g.dz;
  f(*f2, g, 1);
  g.z -= g.dz;
  slice_fill(*slice1, g, *f1, f2, iso);
  g.z += g.dz;
  for (unsigned i = 2; i < g.nz; i++)
  {
    g.z += g.dz;
    f(*f1, g, i);
    std::swap(f1, f2);
    g.z -= g.dz;
    slice_fill(*slice2, g, *f1, f2, iso);
    g.z += g.dz;
    isosurface_slice(s, *slice1, *slice2);
    std::swap(slice1, slice2);
  }
  slice_fill(*slice2, g, *f2, nullptr, iso);
  isosurface_slice(s, *slice1, *slice2);
}

void isosurface_tetra(Surface& s, CartesianGrid g, const PlaneFunc& f, double iso)
{
  check_grid(g, "isosurface_tetra");
  const CartesianGrid gi = g;
  Slice sa(g.nx, g.ny), sb(g.nx, g.ny);
  Slice* slice1 = &sa;
  Slice* slice2 = &sb;
  Helper helper(int(g.nx), int(g.ny));

  unsigned z = 0;
  f(*slice1, g, z);
  iso_sub(*slice1, iso);
  z++;
  g.z += g.dz;
  while (z < g.nz)
  {
    f(*slice2, g, z);
    iso_sub(*slice2, iso);
    g.z += g.dz;
    iso_slice_evaluate(*slice1, *slice2, gi, int(z) - 1, s, helper);
    std::swap(slice1, slice2);
    helper.advance();
    ++z;
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

  CartesianGrid gi = g;
  gi.x -= g.dx;
  gi.y -= g.dy;
  gi.z -= g.dz;
  gi.nx = g.nx + 2;
  gi.ny = g.ny + 2;
  gi.nz = g.nz;
  Helper helper(int(gi.nx), int(gi.ny));

  unsigned z = 0;
  while (z < g.nz)
  {
    f(transfer, g, z);
    g.z += g.dz;
    copy_to_bounded(*slice2, transfer, iso, -1.);
    iso_slice_evaluate(*slice1, *slice2, gi, int(z), s, helper);
    std::swap(slice1, slice2);
    helper.advance();
    ++z;
  }
  slice_init(*slice2, -1.0);
  iso_slice_evaluate(*slice1, *slice2, gi, int(z), s, helper);
}

void isosurface_tetra_bcl(Surface& s, CartesianGrid g, const PlaneFunc& f, double iso)
{
  check_grid(g, "isosurface_tetra_bcl");
  const CartesianGrid gi = g;
  Slice sa(g.nx, g.ny), sb(g.nx, g.ny), sc(g.nx, g.ny);
  Slice* slice1 = &sa;
  Slice* slice2 = &sb;
  Slice* slice3 = &sc;
  HelperBcl helper;

  unsigned z = 0;
  f(*slice1, g, z);
  iso_sub(*slice1, iso);
  z++;
  g.z += g.dz;
  f(*slice2, g, z);
  // GTS subtracts iso from slice1 again here, and never from slice2.
  iso_sub(*slice1, iso);
  z++;
  g.z += g.dz;
  while (z < g.nz)
  {
    f(*slice3, g, z);
    iso_sub(*slice3, iso);
    g.z += g.dz;
    iso_slice_evaluate_bcl(*slice1, *slice2, *slice3, gi, int(z) - 2, s, helper);
    Slice* hs = slice1;
    slice1 = slice2;
    slice2 = slice3;
    slice3 = hs;
    helper.advance();
    ++z;
  }
}

} // namespace gts
