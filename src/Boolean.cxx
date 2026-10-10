// Copyright (C) 1999 Stéphane Popinet (GTS: boolean.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/Boolean.h"

#include "gts/BBTree.h"
#include "gts/Predicates.h"
#include "gts/Stats.h"

// See Geometrium.cxx: Surface::edges() is a std::generator.
#if defined(__GNUC__) && ! defined(__clang__)
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif

#include <algorithm>
#include <vector>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>

namespace gts {

namespace {

// GTS kept the state of the algorithm in the reserved pointers and the flags
// of its objects. Here it is in tables by index. Inter grows them when it
// makes an element, never when it reads them, so that a reference into a
// table stays valid within an expression such as next[a] = next[b].
template <class T>
struct Side
{
  std::vector<T> v;
  T&   operator[](Index i) { return v[i]; }
  void grow(std::size_t n) { if (n > v.size()) v.resize(std::max(n, 2 * v.size())); }
};

// A point for the orientation tests: a vertex, or the point o above a face
// that triangulate_face makes. Simulation of Simplicity ranks vertices by
// index, and o after every vertex.
struct Pt
{
  Vec3     p;
  unsigned rank;
};

constexpr unsigned kRankO = std::numeric_limits<unsigned>::max();

// SurfaceInter::sos_ties(), counted in this thread while t_origin gives the
// origin of each vertex: 1 for the first surface, 2 for the second, 4 for
// any other point.
thread_local std::size_t t_sos_ties = 0;
thread_local const std::vector<std::uint8_t>* t_origin = nullptr;

int orient_sos(const Pt& a, const Pt& b, const Pt& c, const Pt& d)
{
  if (t_origin && orient3d(a.p, b.p, c.p, d.p) == 0.)
  {
    auto origin = [](unsigned r) -> unsigned { return r < t_origin->size() ? (*t_origin)[r] : 4; };
    const unsigned m = origin(a.rank) | origin(b.rank) | origin(c.rank) | origin(d.rank);
    if (m == 3) ++t_sos_ties;
  }
  const Vec3* p[4] = { &a.p, &b.p, &c.p, &d.p };
  const unsigned r[4] = { a.rank, b.rank, c.rank, d.rank };
  return orient3d_sos(p, r);
}

double orient(const Pt& a, const Pt& b, const Pt& c, const Pt& d) { return orient3d(a.p, b.p, c.p, d.p); }

Pt vpt(const Geometrium& g, VertexId v) { return { g.vertex(v).p, v.i }; }

// segment_triangle_intersection: where the edge s crosses the face t, if it
// does, with the parameter of the point from s's second vertex.
struct Crossing { Vec3 p; VertexId from, to; double t; };

std::optional<Crossing> segment_triangle_intersection(const Geometrium& g, EdgeId s, FaceId t)
{
  const auto [abc, edges] = g.vertices_edges(t);
  const Pt A = vpt(g, abc[0]), B = vpt(g, abc[1]), C = vpt(g, abc[2]);
  VertexId vD = g.edge(s).v[0], vE = g.edge(s).v[1];
  Pt D = vpt(g, vD), E = vpt(g, vE);

  int ABCE = orient_sos(A, B, C, E);
  int ABCD = orient_sos(A, B, C, D);
  if (ABCE < 0 || ABCD > 0)
  {
    std::swap(D, E);
    std::swap(vD, vE);
    std::swap(ABCE, ABCD);
  }
  if (ABCE < 0 || ABCD > 0) return std::nullopt;
  if (orient_sos(A, D, C, E) < 0) return std::nullopt;
  if (orient_sos(A, B, D, E) < 0) return std::nullopt;
  if (orient_sos(B, C, D, E) < 0) return std::nullopt;
  const double a = orient(A, B, C, E), b = orient(A, B, C, D);
  if (a != b)
  {
    const double c = a / (a - b);
    return Crossing{ { E.p.x + c * (D.p.x - E.p.x), E.p.y + c * (D.p.y - E.p.y), E.p.z + c * (D.p.z - E.p.z) }, vE, vD, c };
  }
  // D and E are in the plane of ABC.
  return Crossing{ { (E.p.x + D.p.x) / 2., (E.p.y + D.p.y) / 2., (E.p.z + D.p.z) / 2. }, vE, vD, 0.5 };
}

// The state of one gts_surface_inter_new.
struct Inter
{
  Inter(Geometrium& geo, Surface* a, Surface* b) : g(geo), cut{ a, b } { grow(); }

  // The tables to the sizes of the element arrays.
  void grow()
  {
    const std::size_t nv = g.vertex_capacity(), ne = g.edge_capacity(), nf = g.face_capacity();
    vertex_face.grow(nv); vflags.grow(nv);
    points.grow(ne); sub_edges.grow(ne); is_inter.grow(ne); next.grow(ne); eflags.grow(ne);
    face_inter.grow(nf);
  }
  EdgeId new_edge(VertexId a, VertexId b) { const EdgeId e = g.add_edge(a, b); grow(); return e; }
  FaceId new_face(EdgeId a, EdgeId b, EdgeId c) { const FaceId f = g.add_face(a, b, c); grow(); return f; }

  Geometrium& g;
  Surface*    cut[2]; // si->s1 and si->s2 before triangulation: the cut faces

  Side<std::vector<VertexId>> points;     // of an edge: the crossings along it, from v1
  Side<std::vector<EdgeId>>   sub_edges;  // of an edge: the pieces from create_edges
  Side<FaceId>                vertex_face; // of a crossing: the face whose crossing made it
  Side<std::vector<EdgeId>>   face_inter; // of a face: the curve edges across it, newest first
  Side<std::uint8_t>          is_inter;   // of an edge: an edge of the curve (EdgeInter)
  std::vector<EdgeId>         edges;      // the curve, newest first

  // The chains of segments while triangulating a face.
  Side<EdgeId>       next;
  Side<std::uint8_t> eflags;  // INTERIOR, RELEVANT
  Side<std::uint8_t> vflags;  // INTERIOR
  static constexpr std::uint8_t INTERIOR = 1, RELEVANT = 2;

  Pt pt(VertexId v) const { return vpt(g, v); }
  VertexId v1(EdgeId e) const { return g.edge(e).v[0]; }
  VertexId v2(EdgeId e) const { return g.edge(e).v[1]; }
  bool is_set(EdgeId e, std::uint8_t f) { return eflags[e.i] & f; }
  void set(EdgeId e, std::uint8_t f)    { eflags[e.i] |= f; }
  void unset(EdgeId e, std::uint8_t f)  { eflags[e.i] &= std::uint8_t(~f); }

  //--- The curve ------------------------------------------------------------

  VertexId intersects(EdgeId e, FaceId t);
  void     add_edge_inter(EdgeId e, FaceId t, VertexId v);
  int      triangle_point_orientation(FaceId t1, FaceId t2, int o1, VertexId p) const;
  int      intersection_orientation(FaceId t1, EdgeId e, FaceId t2) const;
  void     intersect_edges(FaceId t1, FaceId t2);

  //--- Cutting the faces ----------------------------------------------------

  void create_edges(EdgeId s);
  void add_boundary(EdgeId s, EdgeId nxt, std::vector<EdgeId>& boundary);
  void triangulate_face(FaceId t, Surface& surface);

  EdgeId prev_flag(EdgeId s, std::uint8_t flag);
  EdgeId next_flag(EdgeId s, std::uint8_t flag);
  EdgeId next_interior(VertexId v);
  EdgeId prev_interior(VertexId v);
  EdgeId reverse(EdgeId start, bool interior, bool& isloop);
  std::vector<EdgeId> interior_loops(const std::vector<EdgeId>& interior);
  std::vector<EdgeId> boundary_loops(const std::vector<EdgeId>& boundary);

  struct Context { FaceId t; Surface* surface; Pt o; };
  bool   segment_intersects(VertexId p1, VertexId p2, const std::vector<EdgeId>& interior, const Pt& o) const;
  bool   is_inside_wedge(EdgeId s1, EdgeId s2, VertexId p, const Pt& o) const;
  EdgeId connection(VertexId p, const std::vector<EdgeId>& interior, const std::vector<EdgeId>& bloops, const Pt& o);
  double loop_orientation(EdgeId start, VertexId p, const Pt& o);
  void   connect_interior_loop(EdgeId start, std::vector<EdgeId>& interior, std::vector<EdgeId>& bloops, const Pt& o);
  void   triangulate_boundary_interior(const std::vector<EdgeId>& boundary, std::vector<EdgeId> interior,
                                       const Context& c);

  struct Ear { VertexId v1, v2, v3; EdgeId s1, s2, s3; };
  bool   point_in_wedge(VertexId p1, VertexId p2, VertexId p3, VertexId p, bool closed, const Pt& o) const;
  bool   segment_intersects1(VertexId p1, VertexId p2, VertexId p3, VertexId p4, const Pt& o) const;
  EdgeId triangle_intersects_segments(VertexId p1, VertexId p2, VertexId p3, bool closed, EdgeId start, const Pt& o);
  bool   new_ear(EdgeId s, Ear& e, EdgeId start, unsigned sloppy, const Pt& o);
  void   triangulate_loop(EdgeId start, const Context& c);

  // ORIENTED_VERTICES: the ends of s, the one it shares with nxt second.
  std::pair<VertexId, VertexId> oriented(EdgeId s, EdgeId nxt) const
  {
    if (v1(s) == v1(nxt) || v1(s) == v2(nxt)) return { v2(s), v1(s) };
    return { v1(s), v2(s) };
  }
  // ORIENTATION and ORIENTATION_SOS: with o as the third point.
  double ORIENT(VertexId a, VertexId b, VertexId c, const Pt& o) const { return orient(pt(a), pt(b), o, pt(c)); }
  int    ORIENT_SOS(VertexId a, VertexId b, VertexId c, const Pt& o) const { return orient_sos(pt(a), pt(b), o, pt(c)); }
};

//==============================================================================
// The curve
//==============================================================================

int triangle_triangle_orientation(const Pt* p[6], const VertexId id[6])
{
  int o4 = 0, o5 = 0, o6 = 0;
  if (id[3] != id[0] && id[3] != id[1] && id[3] != id[2]) o4 = orient_sos(*p[0], *p[1], *p[2], *p[3]);
  if (id[4] != id[0] && id[4] != id[1] && id[4] != id[2]) o5 = orient_sos(*p[0], *p[1], *p[2], *p[4]);
  if (o4 * o5 < 0) return 0;
  if (id[5] != id[0] && id[5] != id[1] && id[5] != id[2]) o6 = orient_sos(*p[0], *p[1], *p[2], *p[5]);
  if (o4 * o6 < 0 || o5 * o6 < 0) return 0;
  if (o4) return o4;
  if (o5) return o5;
  if ( ! o6) throw std::logic_error("gts::SurfaceInter: two faces with the same vertices");
  return o6;
}

int Inter::triangle_point_orientation(FaceId t1, FaceId t2, int o1, VertexId vp) const
{
  const VertexId id[6] = { v1(g.face(t1).e[0]), v2(g.face(t1).e[0]), g.opposite_vertex(t1),
                           v1(g.face(t2).e[0]), v2(g.face(t2).e[0]), g.opposite_vertex(t2) };
  const Pt P[6] = { pt(id[0]), pt(id[1]), pt(id[2]), pt(id[3]), pt(id[4]), pt(id[5]) };
  const Pt* a[6] = { &P[0], &P[1], &P[2], &P[3], &P[4], &P[5] };
  int o = triangle_triangle_orientation(a, id);
  if (o != 0) return o;
  const Pt* b[6] = { &P[3], &P[4], &P[5], &P[0], &P[1], &P[2] };
  const VertexId idb[6] = { id[3], id[4], id[5], id[0], id[1], id[2] };
  o = triangle_triangle_orientation(b, idb);
  if (o != 0)
  {
    const int o2 = orient_sos(P[3], P[4], P[5], pt(vp));
    return - o * o1 * o2;
  }
  return 0;
}

// Inserts the crossing v of e with t into e's list, in order along e.
void Inter::add_edge_inter(EdgeId e, FaceId t, VertexId v)
{
  VertexId ev1 = v1(e), ev2 = v2(e);
  auto& list = points[e.i];
  vertex_face[v.i] = t;
  if (list.empty())
  {
    list.push_back(v);
    return;
  }
  const Pt p1 = pt(v1(g.face(t).e[0])), p2 = pt(v2(g.face(t).e[0])), p3 = pt(g.opposite_vertex(t));
  const int oref = orient_sos(p1, p2, p3, pt(ev1));
  int o1 = oref;
  std::size_t i = 0;
  for (; i < list.size(); ++i)
  {
    int o2 = triangle_point_orientation(t, vertex_face[list[i].i], oref, ev1);
    if (o2 == 0) o2 = orient_sos(p1, p2, p3, pt(list[i])); // safe sign evaluation failed
    if (o1 * o2 < 0) break;
    ev1 = list[i];
    o1 = o2;
  }
  if (i < list.size())
    list.insert(list.begin() + std::ptrdiff_t(i), v);
  else
  {
    if ( ! (o1 * orient_sos(p1, p2, p3, pt(ev2)) < 0))
      throw std::logic_error("gts::SurfaceInter: cannot place a crossing along its edge");
    list.push_back(v);
  }
}

VertexId Inter::intersects(EdgeId e, FaceId t)
{
  for (VertexId v : points[e.i])
    if (vertex_face[v.i] == t) return v;
  const auto c = segment_triangle_intersection(g, e, t);
  if ( ! c) return VertexId();
  const VertexId v = g.add_interpolated_vertex(c->p, c->from, c->to, c->t);
  grow();
  add_edge_inter(e, t, v);
  return v;
}

// See GTS's misc/orientation.fig.
int Inter::intersection_orientation(FaceId t1, EdgeId e, FaceId t2) const
{
  const auto [v, es] = g.vertices_edges(t1, e);
  const auto w = g.vertices_of(t2);
  return orient_sos(pt(w[0]), pt(w[1]), pt(w[2]), pt(v[1]));
}

// The segment of the curve across t1, a face of the first surface, and t2,
// a face of the second, if they cross.
void Inter::intersect_edges(FaceId t1, FaceId t2)
{
  VertexId vi1, vi2, v;
  EdgeId e1, e2, e;
  auto update = [&](int o) {
    if (o > 0) { vi2 = v; e2 = e; }
    else       { vi2 = vi1; e2 = e1; vi1 = v; e1 = e; }
  };
  const auto E1 = g.face(t1).e, E2 = g.face(t2).e;

  vi1 = intersects(E2[0], t1);
  e1 = E2[0];
  v = intersects(E2[1], t1);
  e = E2[1];
  if ( ! vi1.valid()) { vi1 = v; e1 = e; }
  else if (v.valid()) update(intersection_orientation(t2, E2[1], t1));
  if ( ! vi2.valid())
  {
    v = intersects(E2[2], t1);
    e = E2[2];
    if ( ! vi1.valid()) { vi1 = v; e1 = e; }
    else if (v.valid()) update(intersection_orientation(t2, E2[2], t1));
  }
  for (int k = 0; k < 3; ++k)
  {
    if (vi2.valid()) break;
    v = intersects(E1[k], t2);
    e = E1[k];
    if ( ! vi1.valid()) { vi1 = v; e1 = e; }
    else if (v.valid()) update(- intersection_orientation(t1, E1[k], t2));
  }

  if (vi1.valid() != vi2.valid())
    throw std::runtime_error("gts::SurfaceInter: two faces cross at a single point");
  if (vi1.valid())
  {
    const EdgeId ne = new_edge(vi1, vi2);
    is_inter[ne.i] = 1;
    cut[0]->add_face(t1);
    cut[1]->add_face(t2);
    edges.insert(edges.begin(), ne);
    face_inter[t1.i].insert(face_inter[t1.i].begin(), ne);
    face_inter[t2.i].insert(face_inter[t2.i].begin(), ne);
  }
}

//==============================================================================
// Chains of segments
//==============================================================================

EdgeId Inter::prev_flag(EdgeId s, std::uint8_t flag)
{
  for (EdgeId x : g.edges_of(v1(s)))
    if (x != s && is_set(x, flag)) return x;
  return EdgeId();
}

EdgeId Inter::next_flag(EdgeId s, std::uint8_t flag)
{
  for (EdgeId x : g.edges_of(v2(s)))
    if (x != s && is_set(x, flag)) return x;
  return EdgeId();
}

EdgeId Inter::next_interior(VertexId v)
{
  for (EdgeId x : g.edges_of(v))
    if (v1(x) == v && is_set(x, INTERIOR)) return x;
  return EdgeId();
}

EdgeId Inter::prev_interior(VertexId v)
{
  for (EdgeId x : g.edges_of(v))
    if (v2(x) == v && is_set(x, INTERIOR)) return x;
  return EdgeId();
}

// A reversed copy of the chain from start, linked to it at its ends.
// Returns the copy of the second segment.
EdgeId Inter::reverse(EdgeId start, bool interior, bool& isloop)
{
  EdgeId s = start, prev, rprev, rstart, rstart1;
  do
  {
    if ( ! is_inter[s.i])
      throw std::logic_error("gts::SurfaceInter: reversing a segment that is not on the curve");
    const EdgeId rs = new_edge(v2(s), v1(s));
    is_inter[rs.i] = 1;
    if ( ! rstart.valid())       rstart = rs;
    else if ( ! rstart1.valid()) rstart1 = rs;
    if (interior) set(rs, INTERIOR);
    next[rs.i] = rprev;
    rprev = rs;
    prev = s;
    s = next[s.i];
  } while (s.valid() && s != start);
  if (s == start)
  {
    next[rstart.i] = rprev;
    isloop = true;
  }
  else
  {
    next[rstart.i] = start;
    next[prev.i] = rprev;
    isloop = false;
  }
  return rstart1;
}

std::vector<EdgeId> Inter::interior_loops(const std::vector<EdgeId>& interior)
{
  std::vector<EdgeId> loops; // newest first
  for (EdgeId s : interior)
  {
    if ( ! is_set(s, RELEVANT)) continue;
    EdgeId start = s, end;
    do
    {
      const EdgeId nxt = next_flag(s, INTERIOR);
      unset(s, RELEVANT);
      end = s;
      s = next[s.i] = nxt;
    } while (s.valid() && s != start);

    if (s == start)
      loops.insert(loops.begin(), start);
    else
    {
      s = prev_flag(start, INTERIOR);
      while (s.valid())
      {
        unset(s, RELEVANT);
        next[s.i] = start;
        start = s;
        s = prev_flag(s, INTERIOR);
      }
      const EdgeId nxt = next_flag(end, RELEVANT);
      const EdgeId prv = prev_flag(start, RELEVANT);
      if (prv.valid()) vflags[v1(start).i] |= INTERIOR;
      if (nxt.valid()) vflags[v2(end).i] |= INTERIOR;
      if ( ! nxt.valid() && ! prv.valid())
        loops.insert(loops.begin(), start);
      else
      {
        bool isloop;
        reverse(start, true, isloop);
      }
    }
  }
  return loops;
}

std::vector<EdgeId> Inter::boundary_loops(const std::vector<EdgeId>& boundary)
{
  const EdgeId first = boundary.front();
  for (std::size_t i = 0; i < boundary.size(); ++i)
  {
    const EdgeId s = boundary[i];
    const EdgeId nxt = i + 1 < boundary.size() ? boundary[i + 1] : first;
    const VertexId v = v1(s) == v1(nxt) || v1(s) == v2(nxt) ? v1(s) : v2(s);
    if (vflags[v.i] & INTERIOR)
    {
      const EdgeId intprev = prev_interior(v);
      next[intprev.i] = nxt;
      next[s.i] = next_interior(v);
      vflags[v.i] &= std::uint8_t(~INTERIOR);
    }
    else
      next[s.i] = nxt;
  }

  std::vector<EdgeId> loops; // newest first
  for (EdgeId start : boundary)
  {
    if ( ! is_set(start, RELEVANT)) continue;
    EdgeId s = start;
    do
    {
      unset(s, RELEVANT);
      unset(s, INTERIOR);
      s = next[s.i];
    } while (s != start);
    loops.insert(loops.begin(), start);
  }
  return loops;
}

//==============================================================================
// Connecting the inner loops and triangulating
//==============================================================================

bool Inter::segment_intersects(VertexId p1, VertexId p2, const std::vector<EdgeId>& interior, const Pt& o) const
{
  for (EdgeId s : interior)
  {
    const VertexId p3 = v1(s), p4 = v2(s);
    if (p3 != p1 && p3 != p2 && p4 != p1 && p4 != p2)
    {
      int o1 = ORIENT_SOS(p3, p4, p1, o), o2 = ORIENT_SOS(p3, p4, p2, o);
      if (o1 * o2 < 0)
      {
        o1 = ORIENT_SOS(p1, p2, p3, o);
        o2 = ORIENT_SOS(p1, p2, p4, o);
        if (o1 * o2 < 0) return true;
      }
    }
  }
  return false;
}

bool Inter::is_inside_wedge(EdgeId s1, EdgeId s2, VertexId p, const Pt& o) const
{
  const auto [w1, w2] = oriented(s1, s2);
  const VertexId w3 = v1(s2) != w2 ? v1(s2) : v2(s2);
  if (ORIENT(w1, w2, w3, o) >= 0.)
  {
    if (ORIENT(w1, w2, p, o) <= 0. || ORIENT(w2, w3, p, o) <= 0.)
      return false;
  }
  else if (ORIENT(w1, w2, p, o) <= 0. && ORIENT(w2, w3, p, o) <= 0.)
    return false;
  return true;
}

EdgeId Inter::connection(VertexId p, const std::vector<EdgeId>& interior, const std::vector<EdgeId>& bloops, const Pt& o)
{
  for (EdgeId start : bloops)
  {
    EdgeId s = start;
    do
    {
      const EdgeId nxt = next[s.i];
      const VertexId w2 = v1(s) == v1(nxt) || v1(s) == v2(nxt) ? v1(s) : v2(s);
      if (is_inside_wedge(s, nxt, p, o) && ! segment_intersects(p, w2, interior, o))
        return s;
      s = nxt;
    } while (s != start);
  }
  return EdgeId();
}

double Inter::loop_orientation(EdgeId start, VertexId p, const Pt& o)
{
  EdgeId s = start;
  double r = 0.;
  do
  {
    const EdgeId nxt = next[s.i];
    const auto [w1, w2] = oriented(s, nxt);
    r += ORIENT(p, w1, w2, o);
    s = nxt;
  } while (s != start);
  return r;
}

void Inter::connect_interior_loop(EdgeId start, std::vector<EdgeId>& interior, std::vector<EdgeId>& bloops, const Pt& o)
{
  EdgeId s = start, c;
  do
  {
    if ( ! (c = connection(v2(s), interior, bloops, o)).valid())
      s = next[s.i];
  } while (s != start && ! c.valid());
  if ( ! c.valid())
    throw std::runtime_error("gts::SurfaceInter: an inner loop of the curve cannot be connected to the boundary of its face");
  const EdgeId nxt = next[c.i];
  const VertexId v = v1(c) == v1(nxt) || v1(c) == v2(nxt) ? v1(c) : v2(c);
  const VertexId cv = v2(s);
  bool isloop;
  EdgeId rs = reverse(s, false, isloop);
  if (isloop)
  {
    if (loop_orientation(rs, v, o) < 0.) std::swap(s, rs);
    bloops.insert(bloops.begin(), rs);
  }
  const EdgeId s1 = new_edge(v, cv);
  const EdgeId rs1 = new_edge(cv, v);
  next[c.i] = s1;
  next[rs1.i] = nxt;
  interior.insert(interior.begin(), s1);
  next[s1.i] = next[s.i];
  next[s.i] = rs1;
}

bool Inter::point_in_wedge(VertexId p1, VertexId p2, VertexId p3, VertexId p, bool closed, const Pt& o) const
{
  if (p == p2 || p == p3) return false;
  double o1 = ORIENT(p1, p2, p, o);
  if ((closed && o1 < 0.) || ( ! closed && o1 <= 0.)) return false;
  o1 = ORIENT(p3, p1, p, o);
  if ((closed && o1 < 0.) || ( ! closed && o1 <= 0.)) return false;
  return true;
}

bool Inter::segment_intersects1(VertexId p1, VertexId p2, VertexId p3, VertexId p4, const Pt& o) const
{
  int o1 = ORIENT_SOS(p3, p4, p1, o), o2 = ORIENT_SOS(p3, p4, p2, o);
  if (o1 * o2 > 0) return false;
  o1 = ORIENT_SOS(p1, p2, p3, o);
  o2 = ORIENT_SOS(p1, p2, p4, o);
  if (o1 * o2 > 0) return false;
  return true;
}

EdgeId Inter::triangle_intersects_segments(VertexId p1, VertexId p2, VertexId p3, bool closed, EdgeId start, const Pt& o)
{
  EdgeId s = start;
  do
  {
    const VertexId p4 = v1(s), p5 = v2(s);
    if      (p4 == p1) { if (point_in_wedge(p1, p2, p3, p5, closed, o)) return s; }
    else if (p4 == p2) { if (point_in_wedge(p2, p3, p1, p5, closed, o)) return s; }
    else if (p4 == p3) { if (point_in_wedge(p3, p1, p2, p5, closed, o)) return s; }
    else if (p5 == p1) { if (point_in_wedge(p1, p2, p3, p4, closed, o)) return s; }
    else if (p5 == p2) { if (point_in_wedge(p2, p3, p1, p4, closed, o)) return s; }
    else if (p5 == p3) { if (point_in_wedge(p3, p1, p2, p4, closed, o)) return s; }
    else if (segment_intersects1(p1, p2, p4, p5, o) ||
             segment_intersects1(p2, p3, p4, p5, o) ||
             segment_intersects1(p3, p1, p4, p5, o))
      return s;
    s = next[s.i];
  } while (s != start);
  return EdgeId();
}

bool Inter::new_ear(EdgeId s, Ear& e, EdgeId start, unsigned sloppy, const Pt& o)
{
  e.s1 = s;
  e.s2 = next[s.i];
  if ( ! e.s2.valid() || e.s2 == e.s1) return false;

  std::tie(e.v1, e.v2) = oriented(e.s1, e.s2);
  e.v3 = v1(e.s2) != e.v2 ? v1(e.s2) : v2(e.s2);
  if (e.v3 == e.v1) return false;
  e.s3 = next[e.s2.i];
  if ( ! e.s3.valid())
    throw std::logic_error("gts::SurfaceInter: a loop of a cut face is broken");
  const bool connects = (v1(e.s3) == e.v1 && v2(e.s3) == e.v3) || (v1(e.s3) == e.v3 && v2(e.s3) == e.v1);
  if (connects)
  {
    if (next[e.s3.i] != e.s1) return false;
  }
  else if (g.find_edge(e.v1, e.v3).valid())
    return false;
  else
    e.s3 = EdgeId();

  const double r = ORIENT(e.v1, e.v2, e.v3, o);
  switch (sloppy)
  {
    case 0:
      if (r <= 0. || triangle_intersects_segments(e.v1, e.v2, e.v3, true, start, o).valid()) return false;
      break;
    case 1:
      if (r < 0. || (r > 0. && triangle_intersects_segments(e.v1, e.v2, e.v3, false, start, o).valid())) return false;
      break;
    case 2:
      if ((r > 0. && triangle_intersects_segments(e.v1, e.v2, e.v3, false, start, o).valid()) ||
          (r < 0. && triangle_intersects_segments(e.v2, e.v1, e.v3, false, start, o).valid()))
        return false;
      break;
    case 3:
      if (r < 0.) return false;
      break;
  }
  if ( ! (r > -1e-6))
    throw std::runtime_error("gts::SurfaceInter: an ear of a cut face is turned over");
  return true;
}

void Inter::triangulate_loop(EdgeId start, const Context& c)
{
  EdgeId prev = start, s = next[start.i];
  unsigned sloppy = 0;
  while (next[s.i] != s)
  {
    const EdgeId nxt = next[s.i];
    Ear e;
    if ( ! new_ear(s, e, start, sloppy, c.o))
    {
      if (s == start) sloppy++;
      prev = s;
      s = nxt;
    }
    else
    {
      if ( ! e.s3.valid()) e.s3 = new_edge(e.v1, e.v3);
      const FaceId f = new_face(e.s1, e.s2, e.s3);
      g.split_attributes(f, c.t);
      c.surface->add_face(f);
      unset(e.s1, RELEVANT); unset(e.s1, INTERIOR);
      unset(e.s2, RELEVANT); unset(e.s2, INTERIOR);
      next[prev.i] = e.s3;
      next[e.s3.i] = next[e.s2.i];
      next[e.s1.i] = next[e.s2.i] = EdgeId();
      start = prev;
      s = next[prev.i];
      sloppy = 0;
    }
  }
  unset(s, RELEVANT);
  unset(s, INTERIOR);
  next[s.i] = EdgeId();
}

void Inter::triangulate_boundary_interior(const std::vector<EdgeId>& boundary, std::vector<EdgeId> interior,
                                          const Context& c)
{
  for (EdgeId s : boundary) set(s, RELEVANT);
  for (EdgeId s : interior) { set(s, RELEVANT); set(s, INTERIOR); }

  const std::vector<EdgeId> iloops = interior_loops(interior);
  std::vector<EdgeId> bloops = boundary_loops(boundary);
  for (EdgeId l : iloops)
    connect_interior_loop(l, interior, bloops, c.o);
  for (EdgeId l : bloops)
    triangulate_loop(l, c);

  // merge_duplicate: the reversed copies of the curve take the place of
  // their originals in the faces.
  for (EdgeId e : interior)
  {
    const EdgeId dup = g.duplicate_of(e);
    if ( ! dup.valid())
      throw std::logic_error("gts::SurfaceInter: a segment of the curve has no reversed copy");
    g.replace_edge(dup, e);
    g.kill_edge(dup);
  }
}

// The pieces of each edge between its crossings, from the edge's first
// vertex. The pieces get the edge's attributes by the split hooks.
void Inter::create_edges(EdgeId s)
{
  // A copy: making edges grows the tables.
  const std::vector<VertexId> list = points[s.i];
  if (list.empty()) return;
  std::vector<EdgeId> pieces;
  VertexId a = list.front();
  pieces.push_back(new_edge(v1(s), a));
  g.split_attributes(pieces.back(), s);
  for (std::size_t i = 0; i < list.size(); ++i)
  {
    const VertexId b = i + 1 < list.size() ? list[i + 1] : v2(s);
    vertex_face[list[i].i] = FaceId();
    pieces.push_back(new_edge(a, b));
    g.split_attributes(pieces.back(), s);
    a = b;
  }
  sub_edges[s.i] = std::move(pieces);
}

// The pieces of s, in the direction that leads to nxt, prepended to the
// boundary.
void Inter::add_boundary(EdgeId s, EdgeId nxt, std::vector<EdgeId>& boundary)
{
  const auto& pieces = sub_edges[s.i];
  if (pieces.empty())
    boundary.insert(boundary.begin(), s);
  else if (v2(s) == v2(nxt) || v2(s) == v1(nxt))
    for (auto i = pieces.rbegin(); i != pieces.rend(); ++i) boundary.insert(boundary.begin(), *i);
  else
    for (EdgeId x : pieces) boundary.insert(boundary.begin(), x);
}

// The pieces of the cut face t, triangulated into surface.
void Inter::triangulate_face(FaceId t, Surface& surface)
{
  std::vector<EdgeId> interior = std::move(face_inter[t.i]);
  face_inter[t.i].clear();
  const Vec3 n = g.normal(t);
  if (n.x == 0. && n.y == 0. && n.z == 0.)
    throw std::runtime_error("gts::SurfaceInter: a cut face has no area");
  const Vec3 p = g.vertex(v1(g.face(t).e[0])).p;
  const Context c = { t, &surface, { p + n, kRankO } };

  const auto E = g.face(t).e;
  std::vector<EdgeId> boundary;
  add_boundary(E[2], E[0], boundary);
  add_boundary(E[1], E[2], boundary);
  add_boundary(E[0], E[1], boundary);
  triangulate_boundary_interior(boundary, std::move(interior), c);
}

//==============================================================================
// Checks and the walk over the faces
//==============================================================================

// triangle_orientation: 1 if t runs along e, -1 if against it.
int triangle_orientation(const Geometrium& g, FaceId t, EdgeId e)
{
  const auto& E = g.face(t).e;
  const EdgeId s = E[0] == e ? E[1] : E[1] == e ? E[2] : E[0];
  const VertexId v2 = g.edge(e).v[1];
  return g.edge(s).v[0] == v2 || g.edge(s).v[1] == v2 ? 1 : -1;
}

// check_orientation: e has exactly two faces of s, oriented oppositely.
bool check_orientation(const Geometrium& g, EdgeId e, const Surface& s)
{
  FaceId t1, t2;
  int o1 = 0;
  for (FaceId t : g.faces_of(e))
  {
    if ( ! s.contains(t)) continue;
    if ( ! t1.valid()) { t1 = t; o1 = triangle_orientation(g, t1, e); }
    else if ( ! t2.valid())
    {
      t2 = t;
      if ( ! (o1 * triangle_orientation(g, t2, e) < 0)) return false;
    }
    else
      throw std::logic_error("gts::SurfaceInter::check: a segment of the curve has more than two faces");
  }
  return t1.valid() && t2.valid();
}

// next_compatible_face: across e from f, a face of s1, or else one of s2
// oriented as f.
FaceId next_compatible_face(const Geometrium& g, EdgeId e, FaceId f, const Surface& s1, const Surface& s2)
{
  FaceId f2, f3;
  for (FaceId f1 : g.faces_of(e))
  {
    if (f1 == f) continue;
    if (s1.contains(f1)) return f1;
    if (s2.contains(f1))
    {
      if      ( ! f2.valid()) f2 = f1;
      else if ( ! f3.valid()) f3 = f1;
      else throw std::runtime_error("gts::SurfaceInter::boolean: a surface is not a manifold");
    }
  }
  if ( ! f3.valid())
  {
    // gts_edge_is_boundary (e, s2)
    if (face_number(s2, e) == 1) return FaceId();
    return f2;
  }
  if ( ! s1.contains(f))
    throw std::runtime_error("gts::SurfaceInter::boolean: three faces on an edge of an uncut face");
  return g.compatible(f, f2, e) ? f2 : f3;
}

} // namespace

//==============================================================================
// SurfaceInter
//==============================================================================

SurfaceInter::SurfaceInter(Surface& a, Surface& b) :
  m_geo(a.geometrium()), m_surface{ &a, &b }, m_pieces{ nullptr, nullptr }
{
  if (&b.geometrium() != &m_geo)
    throw std::invalid_argument("gts::SurfaceInter: the surfaces are in different Geometria");
  Geometrium& g = m_geo;
  std::vector<std::uint8_t> origin(g.vertex_capacity());
  for (VertexId v : a.vertices()) origin[v.i] |= 1;
  for (VertexId v : b.vertices()) origin[v.i] |= 2;
  for (auto& o : origin) if (o == 3) o = 4; // a vertex of both surfaces is neither
  t_sos_ties = 0;
  t_origin = &origin;
  struct Reset { ~Reset() { t_origin = nullptr; } } reset;
  Inter in(g, &g.new_surface(), &g.new_surface());
  try
  {
    // surface_inter_new
    const BBTree ta(a), tb(b);
    BBTree::traverse_overlapping(ta, tb, [&](FaceId t1, FaceId t2) { in.intersect_edges(t1, t2); });

    // create_edges for the edges of the cut faces, of each surface.
    for (int k = 0; k < 2; ++k)
    {
      std::vector<EdgeId> es;
      for (EdgeId e : in.cut[k]->edges()) es.push_back(e);
      for (EdgeId e : es) in.create_edges(e);
    }
    for (int k = 0; k < 2; ++k)
    {
      Surface& p = g.new_surface();
      m_pieces[k] = &p;
      const std::vector<FaceId> faces(in.cut[k]->faces().begin(), in.cut[k]->faces().end());
      for (FaceId t : faces) in.triangulate_face(t, p);
      g.remove_surface(*in.cut[k]);
      in.cut[k] = nullptr;
    }
  }
  catch (...)
  {
    for (Surface* s : in.cut)      if (s) g.remove_surface(*s);
    for (Surface* s : m_pieces)    if (s) g.remove_surface(*s);
    throw;
  }
  m_edges = std::move(in.edges);
  m_sos_ties = t_sos_ties;
}

SurfaceInter::~SurfaceInter()
{
  for (Surface* s : m_pieces) m_geo.remove_surface(*s);
}

SurfaceInter::Check SurfaceInter::check() const
{
  const Geometrium& g = m_geo;
  Check r;
  r.closed = ! m_edges.empty();
  std::vector<std::uint8_t> on_curve(g.edge_capacity());
  for (EdgeId e : m_edges) on_curve[e.i] = 1;

  // check_edge
  for (EdgeId s : m_edges)
  {
    unsigned nn = 0;
    const VertexId a = g.edge(s).v[0], b = g.edge(s).v[1];
    for (EdgeId s1 : g.edges_of(a))
    {
      if ( ! r.ok) break;
      if (s1 != s && on_curve[s1.i]) { if (g.edge(s1).v[1] != a) r.ok = false; nn++; }
    }
    for (EdgeId s1 : g.edges_of(b))
    {
      if ( ! r.ok) break;
      if (s1 != s && on_curve[s1.i]) { if (g.edge(s1).v[0] != b) r.ok = false; nn++; }
    }
    if (nn != 2) r.closed = false;
    if ( ! check_orientation(g, s, *m_pieces[0])) r.ok = false;
    if ( ! check_orientation(g, s, *m_pieces[1])) r.ok = false;
  }
  // As in GTS, a failure leaves the curve counted as not closed; see
  // QUIRKS.md.

  // check_surface_edge: each edge of the pieces has one or two of them,
  // and an edge with one also has a face of the uncut surface.
  if (r.closed)
    for (int k = 0; k < 2; ++k)
      for (EdgeId e : m_pieces[k]->edges())
      {
        const unsigned nf = face_number(*m_pieces[k], e);
        if (nf < 1 || nf > 2) r.ok = false;
        if (nf == 1 && face_number(*m_surface[k], e) == 0) r.ok = false;
      }
  return r;
}

void SurfaceInter::boolean(Surface& result, BooleanOp op) const
{
  const Geometrium& g = m_geo;
  int k, orient;
  switch (op)
  {
    case BooleanOp::A_out_B: k = 0; orient = 1;  break;
    case BooleanOp::A_in_B:  k = 0; orient = -1; break;
    case BooleanOp::B_out_A: k = 1; orient = -1; break;
    case BooleanOp::B_in_A:  k = 1; orient = 1;  break;
    default: throw std::invalid_argument("gts::SurfaceInter::boolean: unknown operation");
  }
  const Surface& s = *m_pieces[k];
  const Surface& orig = *m_surface[k];

  std::vector<std::uint8_t> on_curve(g.edge_capacity()), visited(g.face_capacity());
  for (EdgeId e : m_edges) on_curve[e.i] = 1;

  // walk_faces, breadth first from each segment of the curve. The queue is
  // GtsFifo's: a face and the edge it was reached across.
  std::vector<std::pair<FaceId, EdgeId>> queue;
  auto walk = [&](EdgeId e, FaceId f) {
    queue.assign(1, { f, e });
    for (std::size_t head = 0; head < queue.size(); ++head)
    {
      std::tie(f, e) = queue[head];
      if (visited[f.i]) continue;
      result.add_face(f);
      visited[f.i] = 1;
      for (EdgeId x : g.face(f).e)
      {
        if (x == e || on_curve[x.i]) continue;
        const FaceId f1 = next_compatible_face(g, x, f, s, orig);
        if (f1.valid()) queue.emplace_back(f1, x);
      }
    }
  };
  for (EdgeId e : m_edges)
    for (FaceId t : g.faces_of(e))
      if (s.contains(t) && orient * triangle_orientation(g, t, e) > 0)
      {
        walk(e, t);
        break;
      }
}

//==============================================================================
// Self-intersection
//==============================================================================

bool foreach_intersecting_face(const Surface& s, const std::function<void(FaceId, FaceId)>& f)
{
  const Geometrium& g = s.geometrium();
  auto touch = [&](EdgeId a, EdgeId b) {
    const Edge& A = g.edge(a);
    const Edge& B = g.edge(b);
    return A.v[0] == B.v[0] || A.v[0] == B.v[1] || A.v[1] == B.v[0] || A.v[1] == B.v[1];
  };
  bool any = false;
  const BBTree tree(s);
  BBTree::traverse_overlapping(tree, tree, [&](FaceId t1, FaceId t2) {
    if (t1 == t2) return;
    const auto E1 = g.face(t1).e;
    for (EdgeId x : g.face(t2).e)
      if ( ! touch(x, E1[0]) && ! touch(x, E1[1]) && ! touch(x, E1[2]) &&
          segment_triangle_intersection(g, x, t1))
      {
        any = true;
        f(t1, t2);
        return;
      }
  });
  return any;
}

bool is_self_intersecting(const Surface& s)
{
  return foreach_intersecting_face(s, [](FaceId, FaceId) {});
}

} // namespace gts
