// Copyright (C) 1999 Stéphane Popinet (GTS: surface.c, segment.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/Refine.h"

#include "gts/EHeap.h"
#include "gts/Stats.h"

// See Geometrium.cxx: Surface::edges() is a std::generator.
#if defined(__GNUC__) && ! defined(__clang__)
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif

#include <cmath>
#include <stdexcept>
#include <vector>

namespace gts {

//==============================================================================
// Splitting functions
//==============================================================================

VertexId midvertex(Geometrium& g, EdgeId e)
{
  const auto [v1, v2] = g.edge(e).v;
  const Vec3 p1 = g.vertex(v1).p, p2 = g.vertex(v2).p;
  return g.add_interpolated_vertex({ (p1.x + p2.x) / 2., (p1.y + p2.y) / 2., (p1.z + p2.z) / 2. }, v1, v2, 0.5);
}

VertexId unit_sphere_arc_midvertex(Geometrium& g, EdgeId e)
{
  const auto [v1, v2] = g.edge(e).v;
  const Vec3 p1 = g.vertex(v1).p, p2 = g.vertex(v2).p;
  Vec3 m = { 0.5 * (p1.x + p2.x), 0.5 * (p1.y + p2.y), 0.5 * (p1.z + p2.z) };
  const double n = norm(m);
  m = { m.x / n, m.y / n, m.z / n };
  return g.add_interpolated_vertex(m, v1, v2, 0.5);
}

//==============================================================================
// Tessellation and refinement
//==============================================================================

void tessellate(Surface& s, const Splitter& split)
{
  Geometrium& g = s.geometrium();
  const Splitter& sp = split ? split : Splitter(unit_sphere_arc_midvertex);
  const std::vector<FaceId> faces(s.faces().begin(), s.faces().end());
  SplitMemo memo;
  for (FaceId f : faces)
    g.split_face4(f, sp, memo);
  g.split_remaining(memo);
}

void refine(Surface& s, const RefineStop& stop, const EdgeCost& cost, const Splitter& split)
{
  Geometrium& g = s.geometrium();
  const Splitter& sp = split ? split : Splitter(midvertex);
  EHeap<EdgeId> heap;
  if (cost)
    heap = EHeap<EdgeId>([&](EdgeId e) { return cost(g, e); });
  else
    heap = EHeap<EdgeId>([&](EdgeId e) { const Segment x = g.segment(e); return - distance2(x.a, x.b); });

  heap.freeze();
  for (EdgeId e : s.edges()) heap.insert(e);
  heap.thaw();
  std::vector<EdgeId> fresh;
  while (auto top = heap.remove_top())
  {
    const auto [e, c] = *top;
    if (stop(c, heap.size() + face_number(s, e) + 2)) break;
    // The new edges are costed once they have their faces.
    fresh.clear();
    g.split_edge(e, sp(g, e), [&](EdgeId x) { fresh.push_back(x); });
    for (EdgeId x : fresh) heap.insert(x);
  }
}

//==============================================================================
// Sphere
//==============================================================================

namespace {

constexpr double kX = 0.850650808352039932181540497063011072240401406; // sqrt(sqrt(5)+1)/sqrt(2*sqrt(5))
constexpr double kY = 0.525731112119133606025669084847876607285497935; // sqrt(2)/sqrt(5+sqrt(5))
constexpr double kZ = 0.0;

void generate_icosahedron(Surface& s)
{
  Geometrium& g = s.geometrium();
  const VertexId v01 = g.add_vertex({+kZ, +kX, -kY});
  const VertexId v02 = g.add_vertex({+kX, +kY, +kZ});
  const VertexId v03 = g.add_vertex({+kY, +kZ, -kX});
  const VertexId v04 = g.add_vertex({+kY, +kZ, +kX});
  const VertexId v05 = g.add_vertex({+kX, -kY, +kZ});
  const VertexId v06 = g.add_vertex({+kZ, +kX, +kY});
  const VertexId v07 = g.add_vertex({-kY, +kZ, +kX});
  const VertexId v08 = g.add_vertex({+kZ, -kX, -kY});
  const VertexId v09 = g.add_vertex({-kX, +kY, +kZ});
  const VertexId v10 = g.add_vertex({-kY, +kZ, -kX});
  const VertexId v11 = g.add_vertex({-kX, -kY, +kZ});
  const VertexId v12 = g.add_vertex({+kZ, -kX, +kY});
  const EdgeId e01 = g.add_edge(v01, v02);
  const EdgeId e02 = g.add_edge(v03, v02);
  const EdgeId e03 = g.add_edge(v01, v03);
  const EdgeId e04 = g.add_edge(v04, v05);
  const EdgeId e05 = g.add_edge(v02, v05);
  const EdgeId e06 = g.add_edge(v04, v02);
  const EdgeId e07 = g.add_edge(v06, v07);
  const EdgeId e08 = g.add_edge(v04, v07);
  const EdgeId e09 = g.add_edge(v06, v04);
  const EdgeId e10 = g.add_edge(v08, v03);
  const EdgeId e11 = g.add_edge(v03, v05);
  const EdgeId e12 = g.add_edge(v08, v05);
  const EdgeId e13 = g.add_edge(v06, v09);
  const EdgeId e14 = g.add_edge(v07, v09);
  const EdgeId e15 = g.add_edge(v08, v10);
  const EdgeId e16 = g.add_edge(v03, v10);
  const EdgeId e17 = g.add_edge(v06, v01);
  const EdgeId e18 = g.add_edge(v01, v09);
  const EdgeId e19 = g.add_edge(v08, v11);
  const EdgeId e20 = g.add_edge(v10, v11);
  const EdgeId e21 = g.add_edge(v06, v02);
  const EdgeId e22 = g.add_edge(v12, v11);
  const EdgeId e23 = g.add_edge(v12, v08);
  const EdgeId e24 = g.add_edge(v12, v07);
  const EdgeId e25 = g.add_edge(v07, v11);
  const EdgeId e26 = g.add_edge(v12, v04);
  const EdgeId e27 = g.add_edge(v09, v11);
  const EdgeId e28 = g.add_edge(v10, v09);
  const EdgeId e29 = g.add_edge(v12, v05);
  const EdgeId e30 = g.add_edge(v01, v10);
  s.add_face(g.add_face(e01, e02, e03));
  s.add_face(g.add_face(e04, e05, e06));
  s.add_face(g.add_face(e07, e08, e09));
  s.add_face(g.add_face(e10, e11, e12));
  s.add_face(g.add_face(e13, e14, e07));
  s.add_face(g.add_face(e15, e16, e10));
  s.add_face(g.add_face(e17, e18, e13));
  s.add_face(g.add_face(e19, e20, e15));
  s.add_face(g.add_face(e21, e01, e17));
  s.add_face(g.add_face(e22, e19, e23));
  s.add_face(g.add_face(e09, e06, e21));
  s.add_face(g.add_face(e24, e25, e22));
  s.add_face(g.add_face(e26, e08, e24));
  s.add_face(g.add_face(e20, e27, e28));
  s.add_face(g.add_face(e29, e04, e26));
  s.add_face(g.add_face(e14, e27, e25));
  s.add_face(g.add_face(e23, e12, e29));
  s.add_face(g.add_face(e02, e05, e11));
  s.add_face(g.add_face(e30, e28, e18));
  s.add_face(g.add_face(e03, e16, e30));
}

} // namespace

Surface& generate_sphere(Surface& s, unsigned order)
{
  if (order == 0)
    throw std::invalid_argument("gts::generate_sphere: the order must be at least 1");
  generate_icosahedron(s);
  for (unsigned i = 1; i < order; ++i)
    tessellate(s);
  return s;
}

} // namespace gts
