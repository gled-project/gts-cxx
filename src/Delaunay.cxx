// Copyright (C) 1999 Stéphane Popinet (GTS: cdt.c, refine.c, point.c, triangle.c, segment.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-11.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/Delaunay.h"
#include "gts/EHeap.h"
#include "gts/Predicates.h"
#include "gts/Refine.h"
#include "gts/Stats.h"
#include "gts/Topology.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <deque>
#include <stdexcept>

namespace gts {

//==============================================================================
// Planar tests and constructions
//==============================================================================

Intersect point_in_triangle(Vec3 p, const Triangle& t)
{
  const double d1 = orient2d(t.a, t.b, p);
  if (d1 < 0.0) return Intersect::Out;
  const double d2 = orient2d(t.b, t.c, p);
  if (d2 < 0.0) return Intersect::Out;
  const double d3 = orient2d(t.c, t.a, p);
  if (d3 < 0.0) return Intersect::Out;
  if (d1 == 0.0 || d2 == 0.0 || d3 == 0.0) return Intersect::On;
  return Intersect::In;
}

double in_triangle_circle(Vec3 p, const Triangle& t)
{
  return in_circle(p, t.a, t.b, t.c);
}

Intersect segments_intersect(const Segment& s1, const Segment& s2)
{
  const double d1 = orient2d(s1.a, s1.b, s2.a);
  const double d2 = orient2d(s1.a, s1.b, s2.b);
  if ((d1 > 0.0 && d2 > 0.0) || (d1 < 0.0 && d2 < 0.0)) return Intersect::Out;
  const double d3 = orient2d(s2.a, s2.b, s1.a);
  const double d4 = orient2d(s2.a, s2.b, s1.b);
  if ((d3 > 0.0 && d4 > 0.0) || (d3 < 0.0 && d4 < 0.0)) return Intersect::Out;
  if (d1 == 0.0 || d2 == 0.0 || d3 == 0.0 || d4 == 0.0) return Intersect::On;
  return Intersect::In;
}

std::optional<Vec3> circumcircle_center(const Triangle& t)
{
  const double xa = t.a.x, ya = t.a.y, xb = t.b.x, yb = t.b.y, xc = t.c.x, yc = t.c.y;
  const double xd = (xa + xb)/2., yd = (ya + yb)/2.;
  const double xe = (xa + xc)/2., ye = (ya + yc)/2.;
  const double xad = xd - xa, yad = yd - ya;
  const double xae = xe - xa, yae = ye - ya;
  const double det = xad*yae - xae*yad;
  if (det == 0.) return std::nullopt;
  return Vec3{ (yae*yad*(yd - ye) + xad*yae*xd - xae*yad*xe)/det,
               -(xae*xad*(xd - xe) + yad*xae*yd - yae*xad*ye)/det,
               0. };
}

double interpolate_height(const Geometrium& g, FaceId f, Vec3 p)
{
  const Edge& e1 = g.edge(g.face(f).e[0]);
  const Vec3 p1 = g.vertex(e1.v[0]).p, p2 = g.vertex(e1.v[1]).p, p3 = g.vertex(g.opposite_vertex(f)).p;
  const double x1 = p2.x - p1.x, y1 = p2.y - p1.y;
  const double x2 = p3.x - p1.x, y2 = p3.y - p1.y;
  const double det = x1*y2 - x2*y1;
  if (det == 0.) return (p1.z + p2.z + p3.z)/3.;
  const double x = p.x - p1.x, y = p.y - p1.y;
  const double a = (x*y2 - y*x2)/det;
  const double b = (y*x1 - x*y1)/det;
  return (1. - a - b)*p1.z + a*p2.z + b*p3.z;
}

FaceId enclosing_triangle(Geometrium& g, std::span<const Vec3> points, double scale)
{
  // GTS's value, to 12 digits.
  constexpr double kSqrt3 = 1.73205080757;
  if (points.empty()) return {};
  double xmax = points[0].x, xmin = xmax, ymax = points[0].y, ymin = ymax;
  for (std::size_t i = 1; i < points.size(); ++i)
  {
    const Vec3& p = points[i];
    if (p.x > xmax) xmax = p.x;
    else if (p.x < xmin) xmin = p.x;
    if (p.y > ymax) ymax = p.y;
    else if (p.y < ymin) ymin = p.y;
  }
  const double xo = (xmax + xmin)/2., yo = (ymax + ymin)/2.;
  double r = scale*std::sqrt((xmax - xo)*(xmax - xo) + (ymax - yo)*(ymax - yo));
  if (r == 0.0) r = scale;
  const VertexId v1 = g.add_vertex({ xo + r*kSqrt3, yo - r, 0.0 });
  const VertexId v2 = g.add_vertex({ xo, yo + 2.*r, 0.0 });
  const VertexId v3 = g.add_vertex({ xo - r*kSqrt3, yo - r, 0.0 });
  const EdgeId e1 = g.add_edge(v1, v2);
  const EdgeId e2 = g.add_edge(v2, v3);
  const EdgeId e3 = g.add_edge(v3, v1);
  return g.add_face(e1, e2, e3);
}

//==============================================================================
// Point location
//==============================================================================

namespace {

Vec3 pos(const Geometrium& g, VertexId v) { return g.vertex(v).p; }

// neighbor in cdt.c: a face of s on e other than f.
FaceId neighbor(const Surface& s, FaceId f, EdgeId e)
{
  for (FaceId t : s.geometrium().faces_of(e))
    if (t != f && s.contains(t)) return t;
  return FaceId();
}

// closest_face in cdt.c: of the first n^1/3 faces of s, the counterclockwise
// one whose first vertex is closest to p. GTS took them in the order of its
// hash table; gts-cxx takes them in index order.
FaceId closest_face(const Surface& s, Vec3 p)
{
  const Geometrium& g = s.geometrium();
  double dmin = DBL_MAX;
  FaceId closest;
  int stop = int(std::exp(std::log(double(s.face_number()))/3.));
  for (FaceId f : s.faces())
  {
    if (g.orientation(f) > 0.)
    {
      const Vec3 p1 = pos(g, g.edge(g.face(f).e[0]).v[0]);
      const double d = (p.x - p1.x)*(p.x - p1.x) + (p.y - p1.y)*(p.y - p1.y);
      if (d < dmin)
      {
        dmin = d;
        closest = f;
      }
    }
    if (--stop <= 0) break;
  }
  return closest;
}

// The edge of t that the segment from o, inside t, to p crosses; invalid
// when p is in t too (on_summit false) or when the segment leaves t through
// a vertex (on_summit true).
EdgeId triangle_next_edge(const Geometrium& g, FaceId t, Vec3 o, Vec3 p, bool& on_summit)
{
  const auto ve = g.vertices_edges(t);
  const Vec3 p1 = pos(g, ve.v[0]), p2 = pos(g, ve.v[1]), p3 = pos(g, ve.v[2]);
  on_summit = false;
  double orient = orient2d(o, p1, p);
  if (orient > 0.0)
  {
    orient = orient2d(o, p2, p);
    if (orient > 0.0)
    {
      if (orient2d(p2, p3, p) >= 0.0) return {};
      return ve.e[1];
    }
    if (orient < 0.0)
    {
      if (orient2d(p1, p2, p) >= 0.0) return {};
      return ve.e[0];
    }
    if (orient2d(p1, p2, p) < 0.0) on_summit = true;
    return {};
  }
  if (orient < 0.0)
  {
    orient = orient2d(o, p3, p);
    if (orient > 0.0)
    {
      if (orient2d(p3, p1, p) >= 0.0) return {};
      return ve.e[2];
    }
    if (orient < 0.0)
    {
      if (orient2d(p2, p3, p) >= 0.0) return {};
      return ve.e[1];
    }
    if (orient2d(p3, p1, p) < 0.0) on_summit = true;
    return {};
  }
  if (orient2d(p2, p3, p) < 0.0) return ve.e[1];
  if (orient2d(p1, p2, p) < 0.0) on_summit = true;
  return {};
}

// triangle_barycenter: the x and y of o.
void barycenter(const Geometrium& g, FaceId t, Vec3& o)
{
  const Vec3 p = pos(g, g.opposite_vertex(t));
  const Edge& e1 = g.edge(g.face(t).e[0]);
  const Vec3 a = pos(g, e1.v[0]), b = pos(g, e1.v[1]);
  o.x = (p.x + a.x + b.x)/3.;
  o.y = (p.y + a.y + b.y)/3.;
}

// point_locate in cdt.c: the walk from o in f towards p. A walk that meets a
// vertex restarts from the centre of a neighbouring face.
FaceId walk(const Surface& s, Vec3 o, Vec3 p, FaceId f)
{
  const Geometrium& g = s.geometrium();
  for (;;)
  {
    bool on_summit;
    EdgeId prev = triangle_next_edge(g, f, o, p, on_summit);
    if ( ! prev.valid())
    {
      if ( ! on_summit) return f; // p is inside f
      // The segment leaves f exactly through a vertex: restart from a neighbour.
      const auto& E = g.face(f).e;
      FaceId f1;
      if ((f1 = neighbor(s, f, E[0])).valid() || (f1 = neighbor(s, f, E[1])).valid() ||
          (f1 = neighbor(s, f, E[2])).valid())
      {
        barycenter(g, f1, o);
        f = f1;
        continue;
      }
      return {};
    }

    f = neighbor(s, f, prev);
    VertexId v1, v2, v3;
    EdgeId e2, e3;
    if (f.valid())
    {
      const auto ve = g.vertices_edges(f, prev);
      v1 = ve.v[0]; v2 = ve.v[1]; v3 = ve.v[2];
      prev = ve.e[0]; e2 = ve.e[1]; e3 = ve.e[2];
    }
    bool restart = false;
    while (f.valid())
    {
      const double orient = orient2d(o, pos(g, v3), p);
      if (orient < 0.0)
      {
        if (orient2d(pos(g, v2), pos(g, v3), p) >= 0.0) return f;
        f = neighbor(s, f, e2);
        prev = e2;
        v1 = v3;
      }
      else if (orient > 0.0)
      {
        if (orient2d(pos(g, v3), pos(g, v1), p) >= 0.0) return f;
        f = neighbor(s, f, e3);
        prev = e3;
        v2 = v3;
      }
      else
      {
        if (orient2d(pos(g, v2), pos(g, v3), p) >= 0.0) return f;
        // The segment goes exactly through v3: restart from a neighbour.
        FaceId f1;
        if ((f1 = neighbor(s, f, e2)).valid() || (f1 = neighbor(s, f, e3)).valid())
        {
          barycenter(g, f1, o);
          f = f1;
          restart = true;
          break;
        }
        return {};
      }
      // e2, e3 and v3 of the new face
      if (f.valid())
      {
        const auto& E = g.face(f).e;
        if      (prev == E[0]) { e2 = E[1]; e3 = E[2]; }
        else if (prev == E[1]) { e2 = E[2]; e3 = E[0]; }
        else                   { e2 = E[0]; e3 = E[1]; }
        const Edge& X = g.edge(e2);
        v3 = (X.v[0] == v1 || X.v[0] == v2) ? X.v[1] : X.v[0];
      }
    }
    if ( ! restart) return {};
  }
}

} // namespace

FaceId locate(const Surface& s, Vec3 p, FaceId guess)
{
  const Geometrium& g = s.geometrium();
  if ( ! guess.valid())
    guess = closest_face(s, p);
  else if ( ! s.contains(guess) || ! (g.orientation(guess) > 0.))
    return {};
  if ( ! guess.valid()) return {};
  Vec3 o;
  barycenter(g, guess, o);
  return walk(s, o, p, guess);
}

//==============================================================================
// Insertion and removal
//==============================================================================

namespace {

// A triangulation being changed. Every face that enters or leaves s goes
// through add() and remove(), which tell the refinement's heap, as GTS's
// heap surface class did.
struct Cdt
{
  Surface&    s;
  Geometrium& g;
  std::function<void(FaceId)> on_add, on_remove;

  explicit Cdt(Surface& s_) : s(s_), g(s_.geometrium()) {}

  void add(FaceId f)
  {
    s.add_face(f);
    if (on_add) on_add(f);
  }
  void remove(FaceId f)
  {
    if (on_remove) on_remove(f);
    s.remove_face(f);
  }
  // gts_object_destroy on an edge: its faces leave s through remove().
  void destroy(EdgeId e)
  {
    std::vector<FaceId> faces;
    for (FaceId t : g.faces_of(e)) if (s.contains(t)) faces.push_back(t);
    for (FaceId t : faces) if (on_remove) on_remove(t);
    g.destroy(e);
  }
  // An edge from a to b, the existing one if there is.
  EdgeId edge(VertexId a, VertexId b)
  {
    const EdgeId e = g.find_edge(a, b);
    return e.valid() ? e : g.add_edge(a, b);
  }
  FaceId face(EdgeId a, EdgeId b, EdgeId c, FaceId attributes_from)
  {
    const FaceId f = g.add_face(a, b, c);
    g.split_attributes(f, attributes_from);
    return f;
  }
};

// swap_if_in_circle in cdt.c, after misc/swap.fig: swaps e1, the edge of f1
// from v1 to v2 opposite v3, when the vertex across it is inside the
// circle of f1.
void swap_if_in_circle(Cdt& c, FaceId f1, VertexId v1, VertexId v2, VertexId v3,
                       EdgeId e1, EdgeId e2, EdgeId e3)
{
  Geometrium& g = c.g;
  if (is_constraint(g, e1)) return;
  const FaceId f2 = neighbor(c.s, f1, e1);
  if ( ! f2.valid()) return; // e1 is on the boundary of s

  const auto& E = g.face(f2).e;
  EdgeId e4, e5;
  if      (E[0] == e1) { e4 = E[1]; e5 = E[2]; }
  else if (E[1] == e1) { e4 = E[2]; e5 = E[0]; }
  else                 { e4 = E[0]; e5 = E[1]; }
  const Edge& X4 = g.edge(e4);
  const Edge& X1 = g.edge(e1);
  const VertexId v4 = (X4.v[0] == X1.v[0] || X4.v[0] == X1.v[1]) ? X4.v[1] : X4.v[0];

  if (in_circle(pos(g, v4), pos(g, v1), pos(g, v2), pos(g, v3)) > 0.0)
  {
    const EdgeId en = c.edge(v3, v4);
    const FaceId f3 = c.face(en, e5, e2, f1);
    const FaceId f4 = c.face(en, e3, e4, f2);
    c.remove(f1);
    c.remove(f2);
    c.add(f3);
    c.add(f4);
    swap_if_in_circle(c, f3, v4, v2, v3, e5, e2, en);
    swap_if_in_circle(c, f4, v1, v4, v3, e4, en, e3);
  }
}

VertexId add_vertex_to_face(Cdt& c, VertexId v, FaceId f)
{
  Geometrium& g = c.g;
  const auto ve = g.vertices_edges(f);
  const VertexId v1 = ve.v[0], v2 = ve.v[1], v3 = ve.v[2];
  const EdgeId e1 = ve.e[0], e2 = ve.e[1], e3 = ve.e[2];
  if (v == v1 || v == v2 || v == v3) return {}; // v is in s already
  const Vec3 p = pos(g, v);
  for (VertexId w : { v1, v2, v3 })
    if (p.x == pos(g, w).x && p.y == pos(g, w).y) return w;

  const EdgeId e4 = c.edge(v, v1);
  const EdgeId e5 = c.edge(v, v2);
  const EdgeId e6 = c.edge(v, v3);
  // cf. figure misc/swap.fig
  const FaceId nf0 = c.face(e4, e1, e5, f);
  const FaceId nf1 = c.face(e5, e2, e6, f);
  const FaceId nf2 = c.face(e6, e3, e4, f);
  c.remove(f);
  c.add(nf0);
  c.add(nf1);
  c.add(nf2);
  swap_if_in_circle(c, nf0, v1, v2, v, e1, e5, e4);
  swap_if_in_circle(c, nf1, v2, v3, v, e2, e6, e5);
  swap_if_in_circle(c, nf2, v3, v1, v, e3, e4, e6);
  return {};
}

VertexId add_vertex(Cdt& c, VertexId v, FaceId guess)
{
  const FaceId f = locate(c.s, pos(c.g, v), guess);
  if ( ! f.valid()) return v;
  return add_vertex_to_face(c, v, f);
}

// polygon_in_circle in cdt.c: whether a vertex of poly, other than p1, p2
// and p3, lies inside the circle through them. A vertex of the previous
// edge is not tested again.
bool polygon_in_circle(const Geometrium& g, const std::vector<EdgeId>& poly, VertexId p1, VertexId p2, VertexId p3)
{
  VertexId a, b;
  for (EdgeId e : poly)
  {
    const Edge& S = g.edge(e);
    for (VertexId v : S.v)
      if (v != a && v != b && v != p1 && v != p2 && v != p3 &&
          in_circle(pos(g, v), pos(g, p1), pos(g, p2), pos(g, p3)) > 0.)
        return true;
    a = S.v[0];
    b = S.v[1];
  }
  return false;
}

// triangulate_polygon in cdt.c: fills the polygon of the edges in poly,
// a chain, with faces whose circumcircles hold no vertex of the polygon.
// Each face takes the attributes of ref.
void triangulate_polygon(Cdt& c, std::vector<EdgeId> poly, FaceId ref)
{
  Geometrium& g = c.g;
  if (poly.size() < 2) return;

  const EdgeId s = poly[0];
  const Edge& S = g.edge(s);
  const Edge& S1 = g.edge(poly[1]);
  VertexId v1, v2;
  if (S.v[0] == S1.v[0] || S.v[0] == S1.v[1]) { v1 = S.v[1]; v2 = S.v[0]; }
  else                                         { v1 = S.v[0]; v2 = S.v[1]; }

  std::size_t i = 1;
  VertexId v3 = v2;
  bool found = false;
  while (i < poly.size() && ! found)
  {
    const Edge& X = g.edge(poly[i]);
    v3 = X.v[0] == v3 ? X.v[1] : X.v[0];
    if (v3 != v1 && orient2d(pos(g, v1), pos(g, v2), pos(g, v3)) >= 0. &&
        ! polygon_in_circle(g, poly, v1, v2, v3))
      found = true;
    else
      ++i;
  }
  if ( ! found) return;

  const EdgeId e1 = c.edge(v2, v3);
  const EdgeId e2 = c.edge(v3, v1);
  c.add(c.face(s, e1, e2, ref));

  // The polygon from poly[1] to poly[i], closed by e1, and the one from e2
  // to the end.
  std::vector<EdgeId> poly1(poly.begin() + 1, poly.begin() + i + 1);
  if (e1 != poly[i]) poly1.push_back(e1);
  std::vector<EdgeId> poly2;
  if (i + 1 < poly.size())
  {
    if (e2 != poly[i + 1]) poly2.push_back(e2);
    poly2.insert(poly2.end(), poly.begin() + i + 1, poly.end());
  }
  triangulate_polygon(c, std::move(poly1), ref);
  triangulate_polygon(c, std::move(poly2), ref);
}

void remove_vertex(Cdt& c, VertexId v)
{
  const std::vector<FaceId> triangles = vertex_faces(c.g, v);
  FaceId ref;
  for (FaceId t : triangles)
    if (c.s.contains(t)) { ref = t; break; }
  if ( ! ref.valid()) return;
  triangulate_polygon(c, fan_oriented(c.s, v), ref);
  for (FaceId t : triangles)
    if (c.s.contains(t)) c.remove(t);
}

} // namespace

VertexId delaunay_add_vertex_to_face(Surface& s, VertexId v, FaceId f)
{
  Cdt c(s);
  return add_vertex_to_face(c, v, f);
}

VertexId delaunay_add_vertex(Surface& s, VertexId v, FaceId guess)
{
  Cdt c(s);
  return add_vertex(c, v, guess);
}

void delaunay_remove_vertex(Surface& s, VertexId v)
{
  Cdt c(s);
  remove_vertex(c, v);
}

//==============================================================================
// Constraints
//==============================================================================

namespace {

// remove_triangles in cdt.c
void remove_triangles(Cdt& c, EdgeId e)
{
  std::vector<FaceId> faces;
  for (FaceId t : c.g.faces_of(e)) if (c.s.contains(t)) faces.push_back(t);
  for (FaceId t : faces) c.remove(t);
}

// remove_intersected_edge in cdt.c, without recursion: removes the faces
// from f on that the segment from a to b crosses, from its edge e on,
// collecting the edges left on each side. GTS prepended them to its lists;
// here they are appended. Returns the constraints crossed.
std::vector<EdgeId> remove_intersected_edge(Cdt& c, VertexId a, VertexId b, EdgeId e, FaceId f,
                                            std::vector<EdgeId>& left, std::vector<EdgeId>& right)
{
  Geometrium& g = c.g;
  std::vector<EdgeId> constraints;
  for (;;)
  {
    if (is_constraint(g, e)) constraints.push_back(e);

    const auto ve = g.vertices_edges(f, e);
    const VertexId v1 = ve.v[0], v2 = ve.v[1], v3 = ve.v[2];
    const EdgeId e1 = ve.e[1], e2 = ve.e[2];
    const double o1 = orient2d(pos(g, v2), pos(g, v3), pos(g, b));
    const double o2 = orient2d(pos(g, v3), pos(g, v1), pos(g, b));

    auto drop = [&] {
      remove_triangles(c, e);
      if ( ! is_constraint(g, e) && ! g.edge(e).face.valid()) c.destroy(e);
    };
    // NEXT_CUT: continue across edge, leaving edge1 on its side.
    auto next_cut = [&](EdgeId edge, EdgeId edge1, std::vector<EdgeId>& side) {
      const FaceId next = neighbor(c.s, f, edge);
      drop();
      if ( ! next.valid())
        throw std::logic_error("gts::delaunay_add_constraint: the constraint leaves the triangulation");
      side.push_back(edge1);
      e = edge;
      f = next;
    };

    if (o1 == 0.)
    {
      // b is v3.
      drop();
      left.push_back(e2);
      right.push_back(e1);
      return constraints;
    }
    else if (o1 > 0.)
      next_cut(e2, e1, right);
    else if (o2 >= 0.)
      next_cut(e1, e2, left);
    else if (orient2d(pos(g, a), pos(g, b), pos(g, v3)) > 0.)
      next_cut(e1, e2, left);
    else
      next_cut(e2, e1, right);
  }
}

// remove_intersected_vertex in cdt.c: from the face round a that the
// segment from a to b leaves through, removes the faces it crosses. ref is
// the first of them, kept for its attributes.
std::vector<EdgeId> remove_intersected_vertex(Cdt& c, VertexId a, VertexId b,
                                              std::vector<EdgeId>& left, std::vector<EdgeId>& right, FaceId& ref)
{
  Geometrium& g = c.g;
  const VertexId v = a;
  for (FaceId t : vertex_faces(g, v))
  {
    if ( ! c.s.contains(t)) continue;
    auto vs = g.vertices_of(t);
    VertexId v2 = vs[1], v3 = vs[2];
    if (v == vs[1])      { v2 = vs[2]; v3 = vs[0]; }
    else if (v == vs[2]) { v3 = vs[1]; v2 = vs[0]; }

    if (orient2d(pos(g, v), pos(g, v2), pos(g, b)) >= 0. && orient2d(pos(g, v3), pos(g, v), pos(g, b)) >= 0.)
    {
      const double o3 = orient2d(pos(g, v2), pos(g, v3), pos(g, b));
      EdgeId e = g.opposite_edge(t, v);
      const FaceId next = neighbor(c.s, t, e);
      ref = t;
      const auto ve = g.vertices_edges(t, e);
      e = ve.e[0];
      const EdgeId e2 = ve.e[1], e1 = ve.e[2];

      if (o3 >= 0.) return {}; // b is inside t, or on its edge

      {
        KeepOrphans keep(g, { .faces = true, .edges = true });
        c.remove(t);
      }
      left.push_back(e2);
      right.push_back(e1);
      if ( ! next.valid())
        throw std::logic_error("gts::delaunay_add_constraint: the constraint leaves the triangulation");
      return remove_intersected_edge(c, a, b, e, next, left, right);
    }
  }
  throw std::logic_error("gts::delaunay_add_constraint: no face of the first vertex holds the constraint");
}

} // namespace

std::vector<EdgeId> delaunay_add_constraint(Surface& s, EdgeId ce)
{
  Cdt c(s);
  Geometrium& g = c.g;
  set_constraint(g, ce);
  const VertexId v1 = g.edge(ce).v[0], v2 = g.edge(ce).v[1];
  std::vector<EdgeId> left, right, constraints;
  FaceId ref;
  {
    KeepOrphans keep(g, { .edges = true });
    constraints = remove_intersected_vertex(c, v1, v2, left, right, ref);
  }
  std::vector<EdgeId> poly{ ce };
  poly.insert(poly.end(), right.begin(), right.end());
  triangulate_polygon(c, std::move(poly), ref);
  poly = { ce };
  poly.insert(poly.end(), left.rbegin(), left.rend());
  triangulate_polygon(c, std::move(poly), ref);
  if (ref.valid() && ! g.dead(ref) && g.face(ref).surfaces == 0)
  {
    KeepOrphans keep(g, { .edges = true });
    g.destroy(ref);
  }
  return constraints;
}

FaceId delaunay_check(const Surface& s)
{
  const Geometrium& g = s.geometrium();
  for (FaceId t : s.faces())
  {
    const auto v = g.vertices_of(t);
    std::vector<VertexId> list;
    for (VertexId w : v)
      for (VertexId n : vertex_neighbors(g, w, &s))
        if (std::find(list.begin(), list.end(), n) == list.end()) list.push_back(n);
    for (VertexId n : list)
      if (n != v[0] && n != v[1] && n != v[2] &&
          in_circle(pos(g, n), pos(g, v[0]), pos(g, v[1]), pos(g, v[2])) > 0.)
        return t;
  }
  return {};
}

void delaunay_remove_hull(Surface& s)
{
  Geometrium& g = s.geometrium();
  std::vector<EdgeId> b = boundary(s);
  // GTS's list, used as a stack: the front of the list is the back here.
  std::vector<EdgeId> stack(b.rbegin(), b.rend());
  KeepOrphans keep(g, { .edges = true });
  while ( ! stack.empty())
  {
    const EdgeId e = stack.back();
    stack.pop_back();
    if (is_constraint(g, e)) continue;
    if (const FaceId t = boundary_face(g, e, &s); t.valid())
    {
      for (EdgeId x : g.face(t).e)
        if (x != e && ! is_constraint(g, x) && ! boundary_face(g, x, &s).valid())
          stack.push_back(x);
      s.remove_face(t);
    }
    if ( ! g.edge(e).face.valid()) g.destroy(e);
  }
}

//==============================================================================
// Conforming and refinement, from GTS's refine.c
//==============================================================================

bool vertex_encroaches_edge(const Geometrium& g, VertexId v, EdgeId e)
{
  const Vec3 p = pos(g, v), p1 = pos(g, g.edge(e).v[0]), p2 = pos(g, g.edge(e).v[1]);
  return (p1.x - p.x)*(p2.x - p.x) + (p1.y - p.y)*(p2.y - p.y) < 0.0;
}

VertexId edge_is_encroached(const Surface& s, EdgeId e, const Encroaches& encroaches)
{
  const Geometrium& g = s.geometrium();
  for (FaceId f : g.faces_of(e))
    if (s.contains(f))
    {
      const VertexId v = g.opposite_vertex(f, e);
      if (encroaches ? encroaches(s, v, e) : vertex_encroaches_edge(g, v, e)) return v;
    }
  return {};
}

namespace {

// The queue of encroached constraints, each in it at most once while it
// lives (ALREADY_ENCROACHED).
struct Encroached
{
  std::deque<EdgeId> fifo;
  std::vector<bool>  queued;

  void push(EdgeId e)
  {
    if (queued.size() <= e.i) queued.resize(e.i + 1);
    queued[e.i] = true;
    fifo.push_back(e);
  }
  bool contains(EdgeId e) const { return e.i < queued.size() && queued[e.i]; }
};

bool encroached_by(const Surface& s, VertexId v, EdgeId e, const Encroaches& encroaches)
{
  return encroaches ? encroaches(s, v, e) : vertex_encroaches_edge(s.geometrium(), v, e);
}

// vertex_encroaches: queues the constraints opposite v that v encroaches on.
void vertex_encroaches(Cdt& c, VertexId v, Encroached& q, const Encroaches& encroaches)
{
  Geometrium& g = c.g;
  for (FaceId f : vertex_faces(g, v))
    if (c.s.contains(f))
    {
      const EdgeId e = g.opposite_edge(f, v);
      if ( ! q.contains(e) && is_constraint(g, e) && encroached_by(c.s, v, e, encroaches))
        q.push(e);
    }
}

// split_edge in refine.c: where to split e to avoid infinite cycles, after
// Shewchuk's thesis. Next to another constraint at one end, on concentric
// shells of power-of-two radii round that end; otherwise in the middle.
VertexId split_point(Cdt& c, EdgeId e)
{
  Geometrium& g = c.g;
  constexpr double kSqrt2 = 1.41421356237309504880168872420969807856967187;
  EdgeId cn;
  for (FaceId t : g.faces_of(e))
  {
    if ( ! c.s.contains(t)) continue;
    const auto& E = g.face(t).e;
    EdgeId e1, e2;
    if      (E[0] == e) { e1 = E[1]; e2 = E[2]; }
    else if (E[1] == e) { e1 = E[0]; e2 = E[2]; }
    else                { e1 = E[0]; e2 = E[1]; }
    if (is_constraint(g, e1) && ! is_constraint(g, e2))      cn = e1;
    else if (is_constraint(g, e2) && ! is_constraint(g, e1)) cn = e2;
    if (cn.valid()) break;
  }
  if ( ! cn.valid()) return midvertex(g, e);

  // Use power of two concentric shells.
  const VertexId v1 = g.edge(e).v[0], v2 = g.edge(e).v[1];
  const Vec3 p1 = pos(g, v1), p2 = pos(g, v2);
  const double l = std::sqrt((p2.x - p1.x)*(p2.x - p1.x) + (p2.y - p1.y)*(p2.y - p1.y));
  double nearestpower = 1.;
  while (l > kSqrt2*nearestpower) nearestpower *= 2.;
  while (l < kSqrt2*nearestpower/2.) nearestpower /= 2.;
  double split = nearestpower/l/2.;
  if (g.edge(cn).v[0] == v2 || g.edge(cn).v[1] == v2) split = 1. - split;
  return g.add_interpolated_vertex({ (1. - split)*p1.x + split*p2.x,
                                     (1. - split)*p1.y + split*p2.y,
                                     (1. - split)*p1.z + split*p2.z }, v1, v2, split);
}

// split_encroached in refine.c: splits the queued constraints, each into
// two constraints at a new vertex, until the queue is empty or steiner_max
// vertices were added. Returns what is left of steiner_max.
int split_encroached(Cdt& c, Encroached& q, int steiner_max, const Encroaches& encroaches)
{
  Geometrium& g = c.g;
  while (steiner_max-- != 0)
  {
    if (q.fifo.empty()) break;
    const EdgeId s = q.fifo.front();
    q.fifo.pop_front();
    const VertexId v = split_point(c, s);
    const FaceId boundary = boundary_face(g, s, &c.s);
    FaceId f = boundary;

    // Two copies of s, from its ends to v.
    const EdgeId e1 = g.add_edge(g.edge(s).v[0], v);
    const EdgeId e2 = g.add_edge(v, g.edge(s).v[1]);
    for (EdgeId x : { e1, e2 })
    {
      g.flags(x) = g.edge(s).flags & ~kDead;
      g.split_attributes(x, s);
    }
    set_constraint(g, s, false);

    if ( ! f.valid())
      for (FaceId t : g.faces_of(s))
        if (c.s.contains(t)) { f = t; break; } // gts_edge_has_parent_surface
    if ( ! f.valid() || add_vertex_to_face(c, v, f).valid())
      throw std::logic_error("gts::delaunay_conform: cannot add the vertex that splits a constraint");

    if (boundary.valid() && ! g.dead(s)) c.destroy(s);

    vertex_encroaches(c, v, q, encroaches);
    if (edge_is_encroached(c.s, e1, encroaches).valid()) q.push(e1);
    if (edge_is_encroached(c.s, e2, encroaches).valid()) q.push(e2);
  }
  return steiner_max;
}

} // namespace

unsigned delaunay_conform(Surface& s, int steiner_max, const Encroaches& encroaches)
{
  Cdt c(s);
  Encroached q;
  std::vector<EdgeId> edges;
  for (EdgeId e : s.edges()) edges.push_back(e);
  for (EdgeId e : edges)
    if (is_constraint(c.g, e) && edge_is_encroached(s, e, encroaches).valid())
      q.push(e);
  split_encroached(c, q, steiner_max, encroaches);
  return unsigned(q.fifo.size());
}

unsigned delaunay_refine(Surface& s, int steiner_max, const FaceCost& cost, const Encroaches& encroaches)
{
  Cdt c(s);
  Geometrium& g = c.g;
  using Heap = EHeap<FaceId>;
  Heap heap;
  constexpr Heap::Handle kNoHandle = ~Heap::Handle(0);
  std::vector<Heap::Handle> handle; // EHEAP_PAIR, by face index
  auto handle_of = [&](FaceId f) -> Heap::Handle& {
    if (handle.size() <= f.i) handle.resize(f.i + 1, kNoHandle);
    return handle[f.i];
  };
  auto insert = [&](FaceId f) {
    const double key = cost(g, f);
    if (key != 0.) handle_of(f) = heap.insert(f, key);
  };
  for (FaceId f : s.faces()) insert(f);
  c.on_add = insert;
  c.on_remove = [&](FaceId f) {
    Heap::Handle& h = handle_of(f);
    if (h != kNoHandle) { heap.remove(h); h = kNoHandle; }
  };

  Encroached q;
  while (steiner_max-- != 0)
  {
    const auto top = heap.remove_top();
    if ( ! top) break;
    const FaceId f = top->first;
    handle_of(f) = kNoHandle;
    const auto center = circumcircle_center(g.triangle(f));
    if ( ! center)
      throw std::runtime_error("gts::delaunay_refine: a face has no circumcircle");
    const VertexId v = g.add_vertex(*center);
    if (add_vertex(c, v, f).valid())
      throw std::runtime_error("gts::delaunay_refine: a circumcentre is outside the triangulation or on a vertex");

    vertex_encroaches(c, v, q, encroaches);
    if ( ! q.fifo.empty())
    {
      remove_vertex(c, v);
      steiner_max = split_encroached(c, q, steiner_max, encroaches);
    }
  }
  return unsigned(heap.size());
}

} // namespace gts
