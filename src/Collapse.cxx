// The collapse primitives of Geometrium: whether an edge may collapse, whether
// the collapse would fold the surface, and the collapse itself.
//
// Copyright (C) 1999 Stéphane Popinet (GTS: surface.c, edge.c, vertex.c, triangle.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/Geometrium.h"

#include <algorithm>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace gts {

namespace {

// gts_triangle_use_edges: a face on e1 that has e2 and e3 as edges.
bool face_uses_edges(const Geometrium& g, EdgeId e1, EdgeId e2, EdgeId e3)
{
  for (FaceId t : g.faces_of(e1))
  {
    const auto& E = g.face(t).e;
    const bool has2 = E[0] == e2 || E[1] == e2 || E[2] == e2;
    const bool has3 = E[0] == e3 || E[1] == e3 || E[2] == e3;
    if (has2 && has3) return true;
  }
  return false;
}

bool face_has_edge(const Face& f, EdgeId e) { return f.e[0] == e || f.e[1] == e || f.e[2] == e; }

// gts_edge_is_boundary with no surface: the face, if e has exactly one.
FaceId boundary_face(const Geometrium& g, EdgeId e)
{
  const FaceId f = g.edge(e).face;
  return f.valid() && g.next_around(f, e) == f ? f : FaceId();
}

bool vertex_is_boundary(const Geometrium& g, VertexId v)
{
  for (EdgeId e : g.edges_of(v))
    if (boundary_face(g, e).valid()) return true;
  return false;
}

// Whether the faces around v fall into sets that meet only at v, as
// vertex_contacts() counts them, without changing anything. The faces of
// one set are reached by turning round v across its edges, each with two
// faces; an edge with more takes a search.
bool is_contact(const Geometrium& g, VertexId v)
{
  auto at_v = [&](EdgeId e) { return g.edge(e).v[0] == v || g.edge(e).v[1] == v; };
  // Twice the number of faces around v, as each has two edges on v.
  std::size_t twice = 0;
  EdgeId start;
  bool manifold = true;
  for (EdgeId e : g.edges_of(v))
  {
    std::size_t n = 0;
    for (FaceId t : g.faces_of(e)) { (void) t; ++n; }
    if (n > 2) manifold = false;
    if (n > 0 && ! start.valid()) start = e;
    twice += n;
  }
  if ( ! start.valid()) return false;

  if (manifold)
  {
    // The other face on e, and the other edge of t on v.
    auto across = [&](EdgeId e, FaceId t) {
      for (FaceId x : g.faces_of(e)) if (x != t) return x;
      return FaceId();
    };
    auto other_edge = [&](FaceId t, EdgeId e) {
      for (EdgeId x : g.face(t).e) if (x != e && at_v(x)) return x;
      return EdgeId();
    };
    // The faces from t on, crossing the edges after e, until back at t0 or
    // at a boundary; none for a degenerate face.
    const FaceId t0 = *g.faces_of(start).begin();
    auto turn = [&](FaceId t, EdgeId e, bool& closed) -> std::optional<std::size_t> {
      std::size_t n = 0;
      for (; t.valid() && n <= twice; t = across(e, t))
      {
        ++n;
        if ( ! (e = other_edge(t, e)).valid()) return std::nullopt;
        if (across(e, t) == t0) { closed = true; break; }
      }
      return n;
    };
    bool closed = false;
    if (const auto n1 = turn(t0, start, closed))
    {
      std::optional<std::size_t> n2 = 0;
      if ( ! closed) n2 = turn(across(start, t0), start, closed);
      if (n2) return 2 * (*n1 + *n2) < twice;
    }
  }

  std::vector<FaceId> faces;
  for (EdgeId e : g.edges_of(v))
    for (FaceId t : g.faces_of(e))
      if (std::ranges::find(faces, t) == faces.end()) faces.push_back(t);
  std::vector<FaceId> reached{ faces[0] };
  for (std::size_t i = 0; i < reached.size(); ++i)
    for (EdgeId e : g.face(reached[i]).e)
      if (at_v(e))
        for (FaceId t : g.faces_of(e))
          if (std::ranges::find(reached, t) == reached.end()) reached.push_back(t);
  return reached.size() < faces.size();
}

// triangle_vertices_edges in edge.c: the edges of t other than e, the one
// on e's first vertex first, and the vertex opposite e.
struct OtherEdges { VertexId v; EdgeId e1, e2; };

OtherEdges other_edges(const Geometrium& g, FaceId t, EdgeId e)
{
  auto [e1, e2, e3] = g.face(t).e;
  const VertexId v1 = g.edge(e).v[0];
  if      (e1 == e) e1 = e3;
  else if (e2 == e) e2 = e3;
  if (g.edge(e2).v[0] == v1 || g.edge(e2).v[1] == v1) std::swap(e1, e2);
  return { g.edge(e1).v[0] == v1 ? g.edge(e1).v[1] : g.edge(e1).v[0], e1, e2 };
}

// gts_edge_belongs_to_tetrahedron
bool belongs_to_tetrahedron(const Geometrium& g, EdgeId e)
{
  std::vector<FaceId> faces;
  for (FaceId t : g.faces_of(e)) faces.push_back(t);
  for (std::size_t i = 0; i < faces.size(); ++i)
  {
    const OtherEdges a = other_edges(g, faces[i], e);
    for (std::size_t j = i + 1; j < faces.size(); ++j)
    {
      const OtherEdges b = other_edges(g, faces[j], e);
      const EdgeId s5 = g.find_edge(a.v, b.v);
      if (s5.valid() && face_uses_edges(g, a.e1, b.e1, s5) && face_uses_edges(g, a.e2, b.e2, s5))
        return true;
    }
  }
  return false;
}

// points_are_folded in triangle.c.
bool points_are_folded(Vec3 A, Vec3 B, Vec3 C, Vec3 D, double max)
{
  constexpr double kAreaRatioMax2 = 1e8;
  const Vec3 AB = B - A, AC = C - A, AD = D - A;
  const Vec3 n1 = cross(AB, AC), n2 = cross(AD, AB);
  const double nn1 = dot(n1, n1), nn2 = dot(n2, n2);
  if (nn1 >= kAreaRatioMax2 * nn2 || nn2 >= kAreaRatioMax2 * nn1)
    return true;
  const double n1n2 = dot(n1, n2);
  if (n1n2 > 0.)
    return false;
  return n1n2 * n1n2 / (nn1 * nn2) > max;
}

} // namespace

//==============================================================================
// Validity
//==============================================================================

bool Geometrium::collapsible(EdgeId e) const
{
  const VertexId v1 = m_edges[e.i].v[0], v2 = m_edges[e.i].v[1];

  // Every vertex joined to both v1 and v2 must be the third vertex of a face
  // on e; otherwise the collapse would pinch the surface.
  for (EdgeId e1 : edges_of(v1))
  {
    if (e1 == e) continue;
    EdgeId e2;
    for (EdgeId x : edges_of(other_vertex(e1, v1)))
      if (m_edges[x.i].v[0] == v2 || m_edges[x.i].v[1] == v2) { e2 = x; break; }
    if (e2.valid() && ! face_uses_edges(*this, e, e1, e2))
      return false;
  }

  // At a contact vertex the collapse may join two points of the other set
  // of faces, and an expansion cannot tell which set each face was in.
  if (is_contact(*this, v1) || is_contact(*this, v2))
    return false;

  if (const FaceId t = boundary_face(*this, e); t.valid())
  {
    // A lone triangle.
    const auto& E = m_faces[t.i].e;
    if (boundary_face(*this, E[0]).valid() && boundary_face(*this, E[1]).valid() && boundary_face(*this, E[2]).valid())
      return false;
  }
  else
  {
    // An inner edge joining two boundaries, or an edge of a tetrahedron.
    if (vertex_is_boundary(*this, v1) && vertex_is_boundary(*this, v2))
      return false;
    if (belongs_to_tetrahedron(*this, e))
      return false;
  }
  return true;
}

//==============================================================================
// Folds
//==============================================================================

// GTS moved the ends of every edge of v1 and v2 to the new vertex for the
// test and moved them back after. Here every vertex is read through S(),
// which maps v1 and v2 to v1, and positions through P(), which gives v1 the
// new position. The test is a yes or no over pairs of faces, so the order
// in which GTS listed the faces does not matter.
bool Geometrium::collapse_folds(EdgeId e, Vec3 p, double max) const
{
  const VertexId v1 = m_edges[e.i].v[0], v2 = m_edges[e.i].v[1];
  auto S = [&](VertexId v) { return v == v2 ? v1 : v; };
  auto P = [&](VertexId v) { return v == v1 ? p : m_vertices[v.i].p; };

  // The vertex of t other than A and B, after the move.
  auto third = [&](FaceId t, VertexId A, VertexId B) {
    for (VertexId v : vertices_of(t))
      if (S(v) != A && S(v) != B) return S(v);
    throw std::logic_error("gts::Geometrium::collapse_folds: a face does not use both vertices of its edge");
  };
  // gts_triangles_are_folded
  auto folded = [&](const std::vector<FaceId>& faces, VertexId A, VertexId B) {
    for (std::size_t i = 0; i < faces.size(); ++i)
    {
      const VertexId C = third(faces[i], A, B);
      for (std::size_t j = i + 1; j < faces.size(); ++j)
        if (points_are_folded(P(A), P(B), P(C), P(third(faces[j], A, B)), max))
          return true;
    }
    return false;
  };

  // edge_triangles in surface.c: the faces around e1 after the collapse. A
  // face on e disappears, and its other edge e2 merges with e1, bringing
  // its own faces.
  std::vector<FaceId> faces;
  auto around = [&](EdgeId e1) {
    faces.clear();
    for (FaceId t : faces_of(e1))
    {
      if ( ! face_has_edge(m_faces[t.i], e))
      {
        faces.push_back(t);
        continue;
      }
      const auto& E = m_faces[t.i].e;
      const EdgeId e2 = E[0] != e && E[0] != e1 ? E[0] : E[1] != e && E[1] != e1 ? E[1] : E[2];
      for (FaceId u : faces_of(e2))
        if ( ! face_has_edge(m_faces[u.i], e)) faces.push_back(u);
    }
  };

  for (VertexId v : { v1, v2 })
    for (EdgeId e1 : edges_of(v))
    {
      if (e1 == e) continue;
      around(e1);
      if (folded(faces, S(m_edges[e1.i].v[0]), S(m_edges[e1.i].v[1])))
        return true;
    }

  // The edges opposite the new vertex in the faces around it.
  for (VertexId v : { v1, v2 })
    for (EdgeId x : edges_of(v))
      for (FaceId t : faces_of(x))
      {
        if (face_has_edge(m_faces[t.i], e)) continue;
        EdgeId opposite;
        for (EdgeId y : m_faces[t.i].e)
          if (S(m_edges[y.i].v[0]) != v1 && S(m_edges[y.i].v[1]) != v1) { opposite = y; break; }
        faces.clear();
        for (FaceId u : faces_of(opposite)) faces.push_back(u);
        if (folded(faces, m_edges[opposite.i].v[0], m_edges[opposite.i].v[1]))
          return true;
      }
  return false;
}

//==============================================================================
// Collapse
//==============================================================================

void Geometrium::destroy_face(FaceId f)
{
  for (auto& s : m_surfaces)
    if (s->m_faces.erase(f))
    {
      --m_faces[f.i].surfaces;
      call(s->m_face_callouts, FaceEvent{ f, false });
    }
  kill_face(f);
}

void Geometrium::destroy_edge(EdgeId e)
{
  const auto vs = m_edges[e.i].v;
  kill_edge(e);
  for (VertexId v : vs)
    if ( ! dead(v) && ! m_vertices[v.i].edge.valid())
      m_vertices[v.i].flags |= kDead;
}

void Geometrium::replace_vertex(VertexId v, VertexId with)
{
  std::vector<EdgeId> ring;
  for (EdgeId e : edges_of(v)) ring.push_back(e);
  for (EdgeId e : ring)
  {
    Edge& E = m_edges[e.i];
    if (E.v[0] == with || E.v[1] == with)
      throw std::logic_error("gts::Geometrium::replace_vertex: an edge would join a vertex to itself");
    const int k = vertex_slot(E, v);
    E.v[k] = with;
    link_edge_slot(e, k);
  }
  m_vertices[v.i].edge = EdgeId();
  m_vertices[v.i].flags |= kDead;
}

void Geometrium::replace_edge(EdgeId e, EdgeId with)
{
  std::vector<FaceId> ring;
  for (FaceId t : faces_of(e)) ring.push_back(t);
  for (FaceId t : ring)
  {
    auto E = m_faces[t.i].e;
    for (EdgeId& x : E) if (x == e) x = with;
    set_face_edges(t, E[0], E[1], E[2]);
  }
}

// As edge_collapse in GTS's surface.c, after its validity and fold tests,
// with gts_allow_floating_edges set as gts_surface_coarsen sets it.
VertexId Geometrium::collapse_edge(EdgeId e, Vec3 p, const std::function<void(EdgeId)>& on_removed_edge)
{
  const VertexId v1 = m_edges[e.i].v[0], v2 = m_edges[e.i].v[1];
  // The duplicate would join mid to itself.
  if (duplicate_of(e).valid())
    throw std::logic_error("gts::Geometrium::collapse_edge: the edge has a duplicate");

  // Destroying e destroys its faces; their other edges stay, even without
  // faces, to be merged below.
  std::vector<FaceId> faces;
  for (FaceId t : faces_of(e)) faces.push_back(t);
  for (FaceId t : faces) destroy_face(t);
  destroy_edge(e);

  VertexId mid = add_merged_vertex(p, v1, v2);
  replace_vertex(v1, mid);
  replace_vertex(v2, mid);

  // Each edge of mid with a duplicate keeps it, in the order of mid's ring,
  // as GTS walked mid's list while removing from it.
  std::vector<EdgeId> ring;
  for (EdgeId x : edges_of(mid)) ring.push_back(x);
  std::vector<std::array<EdgeId, 2>> merged;
  std::vector<EdgeId> removed;
  for (std::size_t i = 0; i < ring.size(); ++i)
  {
    const EdgeId e1 = ring[i];
    if (dead(e1)) continue;
    for (EdgeId dup; (dup = duplicate_of(e1)).valid(); )
    {
      replace_edge(dup, e1);
      for (auto& c : m_edge_columns) c.merge(e1.i, e1.i, dup.i);
      if (on_removed_edge) on_removed_edge(dup);
      destroy_edge(dup);
      merged.push_back({ e1, dup });
    }
    if ( ! m_edges[e1.i].face.valid())
    {
      // Only with duplicate faces in the surface.
      if (on_removed_edge) on_removed_edge(e1);
      destroy_edge(e1);
      removed.push_back(e1);
    }
  }
  // destroy_edge() marks mid dead when it is left without edges, unless it
  // never had any, as for a lone edge.
  if ( ! dead(mid) && ! m_vertices[mid.i].edge.valid()) m_vertices[mid.i].flags |= kDead;
  if (dead(mid)) mid = VertexId();

  if ( ! m_collapse_callouts->callouts.empty())
    call(m_collapse_callouts, CollapseEvent{ e, { v1, v2 }, mid, faces, merged, removed });
  return mid;
}

} // namespace gts
