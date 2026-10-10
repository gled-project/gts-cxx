// Copyright (C) 1999 Stéphane Popinet (GTS: surface.c, vopt.c, matrix.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/Coarsen.h"

#include "gts/EHeap.h"
#include "gts/Stats.h"

// See Geometrium.cxx: Surface::edges() is a std::generator.
#if defined(__GNUC__) && ! defined(__clang__)
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace gts {

//==============================================================================
// Constrained quadratic optimization, from GTS's matrix.c
//==============================================================================

namespace {

using V3 = std::array<double, 3>;
using M3 = std::array<V3, 3>;

double scalar(const V3& a, const V3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

V3 cross(const V3& a, const V3& b)
{
  return { a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] };
}

// The constraints A x = b found so far, n of them.
struct Constraints
{
  M3       A{};
  V3       b{};
  unsigned n = 0;
};

// gts_matrix_compatible_row: adds the row A1 x = b1 unless it is zero or
// nearly dependent on the rows already there (within 1 degree). The test
// for a zero row is GTS's; see QUIRKS.md.
void compatible_row(Constraints& c, V3 A1, double b1)
{
  constexpr double kCosAlpha2 = 0.999695413509;   // alpha = 1 degree
  constexpr double kSinAlpha2 = 3.04586490453e-4; // alpha = 1 degree

  double na1 = scalar(A1, A1);
  if (na1 == 0.0)
    return;

  na1 = std::sqrt(na1);
  A1[0] /= na1; A1[1] /= na1; A1[2] /= na1; b1 /= na1;

  if (c.n == 1)
  {
    const double a0a1 = scalar(c.A[0], A1);
    if (a0a1 * a0a1 >= kCosAlpha2)
      return;
  }
  else if (c.n == 2)
  {
    const V3 V = cross(c.A[0], c.A[1]);
    const double s = scalar(V, A1);
    if (s * s <= scalar(V, V) * kSinAlpha2)
      return;
  }
  c.A[c.n] = A1;
  c.b[c.n] = b1;
  ++c.n;
}

// gts_matrix_quadratic_optimization: adds the constraints that minimize
// x H x / 2 + c x in the space left free by the constraints so far.
void quadratic_optimization(Constraints& k, const M3& H, const V3& c)
{
  switch (k.n)
  {
    case 0:
      compatible_row(k, H[0], - c[0]);
      compatible_row(k, H[1], - c[1]);
      compatible_row(k, H[2], - c[2]);
      return;
    case 1:
    {
      V3 Q0 = { 0., 0., 0. };
      const V3& a = k.A[0];
      double max = a[0] * a[0];
      int d = 0;
      // A vector orthogonal to the constraint.
      if (a[1] * a[1] > max) { max = a[1] * a[1]; d = 1; }
      if (a[2] * a[2] > max) { max = a[2] * a[2]; d = 2; }
      switch (d)
      {
        case 0: Q0[0] = - a[2] / a[0]; Q0[2] = 1.0; break;
        case 1: Q0[1] = - a[2] / a[1]; Q0[2] = 1.0; break;
        case 2: Q0[2] = - a[0] / a[2]; Q0[0] = 1.0; break;
      }
      // A second one, orthogonal to the first and to the constraint.
      const V3 Q1 = cross(a, Q0);
      compatible_row(k, { scalar(Q0, H[0]), scalar(Q0, H[1]), scalar(Q0, H[2]) }, - scalar(Q0, c));
      compatible_row(k, { scalar(Q1, H[0]), scalar(Q1, H[1]), scalar(Q1, H[2]) }, - scalar(Q1, c));
      return;
    }
    case 2:
    {
      // A vector orthogonal to both constraints.
      const V3 Q = cross(k.A[0], k.A[1]);
      compatible_row(k, { scalar(Q, H[0]), scalar(Q, H[1]), scalar(Q, H[2]) }, - scalar(Q, c));
      return;
    }
    default:
      throw std::logic_error("gts: quadratic_optimization with three constraints already");
  }
}

// gts_matrix3_inverse applied to b: the solution of A x = b.
Vec3 solve(const Constraints& k)
{
  const M3& m = k.A;
  const double det = (m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) -
                      m[0][1] * (m[1][0] * m[2][2] - m[2][0] * m[1][2]) +
                      m[0][2] * (m[1][0] * m[2][1] - m[2][0] * m[1][1]));
  if (k.n != 3 || det == 0.0)
    throw std::logic_error("gts::volume_optimized_vertex: the constraints do not fix the position");
  M3 mi;
  mi[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / det;
  mi[0][1] = (m[2][1] * m[0][2] - m[0][1] * m[2][2]) / det;
  mi[0][2] = (m[0][1] * m[1][2] - m[1][1] * m[0][2]) / det;
  mi[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / det;
  mi[1][1] = (m[0][0] * m[2][2] - m[2][0] * m[0][2]) / det;
  mi[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) / det;
  mi[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) / det;
  mi[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) / det;
  mi[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / det;
  const V3& b = k.b;
  return { mi[0][0] * b[0] + mi[0][1] * b[1] + mi[0][2] * b[2],
           mi[1][0] * b[0] + mi[1][1] * b[1] + mi[1][2] * b[2],
           mi[2][0] * b[0] + mi[2][1] * b[1] + mi[2][2] * b[2] };
}

//==============================================================================
// Neighbourhoods in GTS's list order
//==============================================================================

// GTS built these lists by prepending what it had not seen yet, so they come
// in the reverse order of discovery. The sums below run over them in that
// order, so that they round as in GTS.
template <class Id>
void add_unique(std::vector<Id>& v, Id x)
{
  if (std::find(v.begin(), v.end(), x) == v.end()) v.push_back(x);
}

// gts_vertex_triangles (v1) followed by gts_vertex_triangles (v2).
std::vector<FaceId> vertex_triangles(const Geometrium& g, VertexId v1, VertexId v2)
{
  std::vector<FaceId> r;
  for (VertexId v : { v1, v2 })
    for (EdgeId e : g.edges_of(v))
      for (FaceId t : g.faces_of(e)) add_unique(r, t);
  std::reverse(r.begin(), r.end());
  return r;
}

// gts_vertex_neighbors (v1) followed by gts_vertex_neighbors (v2).
std::vector<VertexId> vertex_neighbors(const Geometrium& g, VertexId v1, VertexId v2)
{
  std::vector<VertexId> r;
  for (VertexId v : { v1, v2 })
    for (EdgeId e : g.edges_of(v)) add_unique(r, g.other_vertex(e, v));
  std::reverse(r.begin(), r.end());
  return r;
}

//==============================================================================
// Volume optimization, from GTS's vopt.c
//==============================================================================

// triangle_normal in vopt.c: the normal by gts_triangle_normal's choice of
// vertices, and nt = (p1 x p2) . p3.
struct NormalNt { double n1, n2, n3, nt; };

NormalNt triangle_normal(const Geometrium& g, FaceId t)
{
  const Edge& e1 = g.edge(g.face(t).e[0]);
  const Edge& e2 = g.edge(g.face(t).e[1]);
  const Vec3 p1 = g.vertex(e1.v[0]).p;
  Vec3 p2, p3;
  if      (e1.v[0] == e2.v[0]) { p2 = g.vertex(e2.v[1]).p; p3 = g.vertex(e1.v[1]).p; }
  else if (e1.v[1] == e2.v[1]) { p2 = g.vertex(e1.v[1]).p; p3 = g.vertex(e2.v[0]).p; }
  else if (e1.v[0] == e2.v[1]) { p2 = g.vertex(e2.v[0]).p; p3 = g.vertex(e1.v[1]).p; }
  else                         { p2 = g.vertex(e1.v[1]).p; p3 = g.vertex(e2.v[1]).p; }

  const double x1 = p2.x - p1.x, y1 = p2.y - p1.y, z1 = p2.z - p1.z;
  const double x2 = p3.x - p1.x, y2 = p3.y - p1.y, z2 = p3.z - p1.z;
  return { y1 * z2 - z1 * y2, z1 * x2 - x1 * z2, x1 * y2 - y1 * x2,
           ((p1.y * p2.z - p1.z * p2.y) * p3.x +
            (p1.z * p2.x - p1.x * p2.z) * p3.y +
            (p1.x * p2.y - p1.y * p2.x) * p3.z) };
}

// The two ends of a boundary edge in the orientation of its face.
std::array<Vec3, 2> boundary_ends(const Geometrium& g, EdgeId edge, FaceId f)
{
  const auto& E = g.face(f).e;
  const EdgeId edge2 = edge == E[0] ? E[1] : edge == E[1] ? E[2] : E[0];
  VertexId v1 = g.edge(edge).v[0], v2 = g.edge(edge).v[1];
  if (v2 != g.edge(edge2).v[0] && v2 != g.edge(edge2).v[1]) { v2 = v1; v1 = g.edge(edge).v[1]; }
  return { g.vertex(v1).p, g.vertex(v2).p };
}

struct Boundary
{
  V3 e1{}, e2{}, cb{};
  M3 Hb{};
};

void boundary_preservation(const Geometrium& g, EdgeId edge, FaceId f, Boundary& B)
{
  const auto [p1, p2] = boundary_ends(g, edge, f);
  const V3 e = { p2.x - p1.x, p2.y - p1.y, p2.z - p1.z };
  B.e1[0] += e[0]; B.e1[1] += e[1]; B.e1[2] += e[2];

  const V3 e3 = { p2.y * p1.z - p2.z * p1.y, p2.z * p1.x - p2.x * p1.z, p2.x * p1.y - p2.y * p1.x };
  B.e2[0] += e3[0]; B.e2[1] += e3[1]; B.e2[2] += e3[2];

  M3& H = B.Hb;
  H[0][0] += e[1] * e[1] + e[2] * e[2];
  H[0][1] -= e[0] * e[1];
  H[0][2] -= e[0] * e[2];
  H[1][0] = H[0][1];
  H[1][1] += e[0] * e[0] + e[2] * e[2];
  H[1][2] -= e[1] * e[2];
  H[2][0] = H[0][2];
  H[2][1] = H[1][2];
  H[2][2] += e[0] * e[0] + e[1] * e[1];

  B.cb[0] += e[1] * e3[2] - e[2] * e3[1];
  B.cb[1] += e[2] * e3[0] - e[0] * e3[2];
  B.cb[2] += e[0] * e3[1] - e[1] * e3[0];
}

double boundary_cost(const Geometrium& g, EdgeId edge, FaceId f, Vec3 p)
{
  const auto [p1, p2] = boundary_ends(g, edge, f);
  const double e0 = (p2.y - p1.y) * (p.z - p2.z) - (p2.z - p1.z) * (p.y - p2.y);
  const double e1 = (p2.z - p1.z) * (p.x - p2.x) - (p2.x - p1.x) * (p.z - p2.z);
  const double e2 = (p2.x - p1.x) * (p.y - p2.y) - (p2.y - p1.y) * (p.x - p2.x);
  return e0 * e0 + e1 * e1 + e2 * e2;
}

// Calls f(edge, face) for the boundary edges around the two ends of e: all
// those of the first end, e included, and those of the second but e.
template <class F>
void foreach_boundary(const Geometrium& g, EdgeId e, F&& f)
{
  for (EdgeId x : g.edges_of(g.edge(e).v[0]))
    if (const FaceId t = boundary_face(g, x); t.valid()) f(x, t);
  for (EdgeId x : g.edges_of(g.edge(e).v[1]))
    if (x != e)
      if (const FaceId t = boundary_face(g, x); t.valid()) f(x, t);
}

} // namespace

Vec3 volume_optimized_vertex(const Geometrium& g, EdgeId edge, const VolumeOptimizedParams& params)
{
  const VertexId v1 = g.edge(edge).v[0], v2 = g.edge(edge).v[1];
  Constraints k;

  // Boundary preservation.
  Boundary B;
  unsigned nb = 0;
  foreach_boundary(g, edge, [&](EdgeId x, FaceId t) { boundary_preservation(g, x, t, B); ++nb; });
  if (nb > 0)
  {
    const V3& e1 = B.e1;
    const V3& e2 = B.e2;
    const M3 H = { V3{ e1[2] * e1[2] + e1[1] * e1[1], - e1[0] * e1[1], - e1[0] * e1[2] },
                   V3{ - e1[0] * e1[1], e1[2] * e1[2] + e1[0] * e1[0], - e1[1] * e1[2] },
                   V3{ - e1[0] * e1[2], - e1[1] * e1[2], e1[1] * e1[1] + e1[0] * e1[0] } };
    const V3 c = { e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0] };
    quadratic_optimization(k, H, c);
  }

  // Volume preservation.
  double sn1 = 0., sn2 = 0., sn3 = 0.;
  double sn11 = 0., sn22 = 0., sn33 = 0.;
  double sn12 = 0., sn13 = 0., sn23 = 0.;
  double st = 0., stn1 = 0., stn2 = 0., stn3 = 0.;
  for (FaceId t : vertex_triangles(g, v1, v2))
  {
    const auto [n1, n2, n3, nt] = triangle_normal(g, t);
    sn1 += n1; sn2 += n2; sn3 += n3;
    sn11 += n1 * n1; sn22 += n2 * n2; sn33 += n3 * n3;
    sn12 += n1 * n2; sn13 += n1 * n3; sn23 += n2 * n3;
    st += nt; stn1 += nt * n1; stn2 += nt * n2; stn3 += nt * n3;
  }
  compatible_row(k, { sn1, sn2, sn3 }, st);

  // Volume and boundary optimization, weighted.
  if (k.n < 3)
  {
    M3 H = { V3{ sn11, sn12, sn13 }, V3{ sn12, sn22, sn23 }, V3{ sn13, sn23, sn33 } };
    V3 c = { - stn1, - stn2, - stn3 };
    const Vec3 a = g.vertex(v1).p, b = g.vertex(v2).p;
    const double le = 9. * params.boundary_weight * distance2(a, b);
    if (nb > 0)
      for (int i = 0; i < 3; ++i)
      {
        for (int j = 0; j < 3; ++j)
          H[i][j] = params.volume_weight * H[i][j] + le * B.Hb[i][j];
        c[i] = params.volume_weight * c[i] + le * B.cb[i];
      }
    quadratic_optimization(k, H, c);
  }

  // Triangle shape optimization.
  if (k.n < 3)
  {
    double nv = 0.0;
    V3 c = { 0., 0., 0. };
    for (VertexId w : vertex_neighbors(g, v1, v2))
      if (w != v1 && w != v2)
      {
        const Vec3 p = g.vertex(w).p;
        nv += 1.0;
        c[0] -= p.x; c[1] -= p.y; c[2] -= p.z;
      }
    const M3 H = { V3{ nv, 0., 0. }, V3{ 0., nv, 0. }, V3{ 0., 0., nv } };
    quadratic_optimization(k, H, c);
  }

  return solve(k);
}

double volume_optimized_cost(const Geometrium& g, EdgeId e, const VolumeOptimizedParams& params)
{
  const Vec3 v = volume_optimized_vertex(g, e, params);
  const VertexId v1 = g.edge(e).v[0], v2 = g.edge(e).v[1];
  const double length2 = distance2(g.vertex(v1).p, g.vertex(v2).p);

  // edge_volume_cost
  double volume = 0.0;
  for (FaceId t : vertex_triangles(g, v1, v2))
  {
    const auto [n1, n2, n3, nt] = triangle_normal(g, t);
    const double a = v.x * n1 + v.y * n2 + v.z * n3 - nt;
    volume += a * a;
  }
  volume /= 36.;

  // edge_boundary_cost
  double boundary = 0.;
  foreach_boundary(g, e, [&](EdgeId x, FaceId t) { boundary += boundary_cost(g, x, t, v); });
  boundary /= 4.;

  // edge_shape_cost
  double shape = 0.;
  for (VertexId w : vertex_neighbors(g, v1, v2))
    if (w != v1 && w != v2) shape += distance2(g.vertex(w).p, v);

  return params.volume_weight * volume +
         params.boundary_weight * length2 * boundary +
         params.shape_weight * length2 * length2 * shape;
}

Vec3 midpoint(const Geometrium& g, EdgeId e)
{
  const Vec3 p1 = g.vertex(g.edge(e).v[0]).p, p2 = g.vertex(g.edge(e).v[1]).p;
  return { (p1.x + p2.x) / 2., (p1.y + p2.y) / 2., (p1.z + p2.z) / 2. };
}

//==============================================================================
// Coarsening
//==============================================================================

CoarsenStop stop_below_edges(std::size_t min_edges)
{
  return [min_edges](double, std::size_t n) { return n < min_edges; };
}

CoarsenStop stop_above_cost(double max_cost)
{
  return [max_cost](double cost, std::size_t) { return cost > max_cost; };
}

void coarsen(Surface& s, const CoarsenStop& stop, const EdgeCost& cost, const CollapsePosition& position, double min_angle)
{
  Geometrium& g = s.geometrium();
  const CollapsePosition& pos = position ? position : CollapsePosition(midpoint);
  using Heap = EHeap<EdgeId>;
  Heap heap;
  if (cost)
    heap = Heap([&](EdgeId e) { return cost(g, e); });
  else
    heap = Heap([&](EdgeId e) { const Segment x = g.segment(e); return distance2(x.a, x.b); });

  // The heap handle of each edge, as GTS kept it in the edge's reserved
  // pointer. Collapses make no edges, so the size stays.
  constexpr Heap::Handle kOut = Heap::Handle(-1);
  std::vector<Heap::Handle> handle(g.edge_capacity(), kOut);
  auto insert = [&](EdgeId e) { handle[e.i] = heap.insert(e); };
  auto remove = [&](EdgeId e) {
    if (handle[e.i] != kOut) { heap.remove(handle[e.i]); handle[e.i] = kOut; }
  };

  double max_cos2 = std::cos(min_angle);
  max_cos2 *= max_cos2;

  heap.freeze();
  for (EdgeId e : s.edges()) insert(e);
  heap.thaw();

  std::vector<EdgeId> second;
  while (auto top = heap.remove_top())
  {
    const auto [e, c] = *top;
    handle[e.i] = kOut;
    if ( ! (c < DBL_MAX) || stop(c, heap.size() - face_number(s, e)))
      break;

    // edge_collapse: a refused edge goes to the end.
    if ( ! g.collapsible(e))
    {
      handle[e.i] = heap.insert(e, DBL_MAX);
      continue;
    }
    const Vec3 p = pos(g, e);
    if (g.collapse_folds(e, p, max_cos2))
    {
      handle[e.i] = heap.insert(e, DBL_MAX);
      continue;
    }
    const VertexId v = g.collapse_edge(e, p, remove);
    if ( ! v.valid())
      continue;

    // update_2nd_closest_neighbors: the edges within two edges of v, in
    // GTS's list order, are taken out and put back with new costs. GTS put
    // in an edge that was not in the heap, after a CRITICAL message.
    second.clear();
    for (EdgeId x : g.edges_of(v))
    {
      const VertexId w = g.other_vertex(x, v);
      for (EdgeId y : g.edges_of(w)) add_unique(second, y);
    }
    std::reverse(second.begin(), second.end());
    for (EdgeId y : second)
    {
      remove(y);
      insert(y);
    }
  }
}

} // namespace gts
