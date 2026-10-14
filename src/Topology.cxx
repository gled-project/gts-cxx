// Copyright (C) 1999 Stéphane Popinet (GTS: vertex.c, edge.c, face.c, triangle.c, segment.c, surface.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-11.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/Topology.h"
#include "gts/Stats.h"

#include <algorithm>

namespace gts {

namespace {

template <class T>
void add_unique(std::vector<T>& v, T x)
{
  if (std::find(v.begin(), v.end(), x) == v.end()) v.push_back(x);
}

bool in(const Surface* s, FaceId f) { return ! s || s->contains(f); }

// Whether e has a face of s.
bool has_face(const Geometrium& g, EdgeId e, const Surface& s)
{
  for (FaceId f : g.faces_of(e))
    if (s.contains(f)) return true;
  return false;
}

// The edges of f in its order, starting from the edge opposite v when v is
// the first vertex, rotated as gts_vertex_fan_oriented does: e2 is the edge
// after v, e3 the one before it.
void rotate_to(const Geometrium::VerticesEdges& ve, VertexId v, EdgeId& e2, EdgeId& e3)
{
  const EdgeId e1 = ve.e[0];
  e2 = ve.e[1];
  e3 = ve.e[2];
  if (v == ve.v[1])
  {
    e2 = e3;
    e3 = e1;
  }
  else if (v == ve.v[2])
  {
    e3 = e2;
    e2 = e1;
  }
}

} // namespace

std::vector<FaceId> vertex_faces(const Geometrium& g, VertexId v, const Surface* s)
{
  std::vector<FaceId> list;
  for (EdgeId e : g.edges_of(v))
    for (FaceId t : g.faces_of(e))
      if (in(s, t)) add_unique(list, t);
  std::reverse(list.begin(), list.end());
  return list;
}

std::vector<VertexId> vertex_neighbors(const Geometrium& g, VertexId v, const Surface* s)
{
  std::vector<VertexId> list;
  for (EdgeId e : g.edges_of(v))
    if ( ! s || has_face(g, e, *s)) add_unique(list, g.other_vertex(e, v));
  std::reverse(list.begin(), list.end());
  return list;
}

// Where GTS's recursion in edge_fan_list stopped on a non-manifold edge or a
// wrongly oriented face, the list ends there.
std::vector<EdgeId> fan_oriented(const Surface& s, VertexId v)
{
  const Geometrium& g = s.geometrium();
  FaceId f;
  unsigned d = 2;
  for (EdgeId e : g.edges_of(v))
  {
    FaceId f1;
    unsigned degree = 0;
    for (FaceId t : g.faces_of(e))
      if (s.contains(t)) { f1 = t; degree++; }
    if ( ! f1.valid()) continue;
    if (degree > 2) return {}; // non-manifold edge
    if (degree == 1)
    {
      EdgeId e2, e3;
      rotate_to(g.vertices_edges(f1), v, e2, e3);
      if (e3 != e)
      {
        d = 1;
        f = f1;
      }
    }
    else if (degree <= d)
      f = f1;
  }
  if ( ! f.valid()) return {};

  EdgeId e2, e3;
  rotate_to(g.vertices_edges(f), v, e2, e3);
  std::vector<EdgeId> list{ e2 };

  // edge_fan_list
  const FaceId first = f;
  EdgeId e = e3;
  for (;;)
  {
    FaceId neighbor;
    bool manifold = true;
    for (FaceId f1 : g.faces_of(e))
      if (f1 != f && s.contains(f1))
      {
        if (neighbor.valid()) { manifold = false; break; }
        neighbor = f1;
      }
    if ( ! manifold || ! neighbor.valid() || neighbor == first) break;

    const auto& E = g.face(neighbor).e;
    EdgeId next, enext;
    if      (E[0] == e) { next = E[1]; enext = E[2]; }
    else if (E[1] == e) { next = E[2]; enext = E[0]; }
    else                { next = E[0]; enext = E[1]; }
    // checking for correct orientation
    if (g.edge(enext).v[0] != v && g.edge(enext).v[1] != v) break;

    list.push_back(next);
    f = neighbor;
    e = enext;
  }
  return list;
}

std::vector<FaceId> face_neighbors(const Geometrium& g, FaceId f, const Surface* s)
{
  std::vector<FaceId> list;
  for (EdgeId e : g.face(f).e)
    for (FaceId t : g.faces_of(e))
      if (t != f && in(s, t)) list.push_back(t);
  std::reverse(list.begin(), list.end());
  return list;
}

bool is_compatible(const Surface& s, FaceId f)
{
  const Geometrium& g = s.geometrium();
  for (EdgeId e : g.face(f).e)
    for (FaceId t : g.faces_of(e))
      if (t != f && s.contains(t) && ! g.compatible(f, t, e)) return false;
  return true;
}

std::optional<std::array<FaceId, 2>> manifold_faces(const Surface& s, EdgeId e)
{
  const Geometrium& g = s.geometrium();
  std::array<FaceId, 2> r;
  int n = 0;
  for (FaceId t : g.faces_of(e))
    if (s.contains(t))
    {
      if (n == 2) return std::nullopt;
      r[n++] = t;
    }
  if (n < 2) return std::nullopt;
  return r;
}

bool belongs_to_tetrahedron(const Geometrium& g, EdgeId e)
{
  // triangle_vertices_edges in edge.c: the vertex of t opposite e and the
  // two other edges, the one on e's first vertex first.
  struct Rest { VertexId v; EdgeId e1, e2; };
  auto rest = [&](FaceId t) {
    const VertexId v1 = g.edge(e).v[0];
    Rest r;
    for (EdgeId x : g.face(t).e)
    {
      if (x == e) continue;
      const Edge& X = g.edge(x);
      if (X.v[0] == v1 || X.v[1] == v1) r.e1 = x; else r.e2 = x;
    }
    r.v = g.opposite_vertex(t, e);
    return r;
  };
  // gts_triangle_use_edges: whether a face has the three edges.
  auto face_with = [&](EdgeId a, EdgeId b, EdgeId c) {
    for (FaceId t : g.faces_of(a))
    {
      const auto& E = g.face(t).e;
      auto has = [&](EdgeId x) { return E[0] == x || E[1] == x || E[2] == x; };
      if (has(b) && has(c)) return true;
    }
    return false;
  };

  std::vector<FaceId> faces;
  for (FaceId t : g.faces_of(e)) faces.push_back(t);
  for (std::size_t i = 0; i < faces.size(); ++i)
  {
    const Rest a = rest(faces[i]);
    for (std::size_t j = i + 1; j < faces.size(); ++j)
    {
      const Rest b = rest(faces[j]);
      const EdgeId s5 = g.find_edge(a.v, b.v);
      if (s5.valid() && face_with(a.e1, b.e1, s5) && face_with(a.e2, b.e2, s5))
        return true;
    }
  }
  return false;
}

unsigned edge_contacts(const Geometrium& g, EdgeId e)
{
  const VertexId v1 = g.edge(e).v[0], v2 = g.edge(e).v[1];
  // The faces around both vertices, marked until reached.
  std::vector<FaceId> marked = vertex_faces(g, v1);
  for (FaceId t : vertex_faces(g, v2)) add_unique(marked, t);
  auto take = [&](FaceId t) {
    auto it = std::find(marked.begin(), marked.end(), t);
    if (it == marked.end()) return false;
    marked.erase(it);
    return true;
  };
  auto uses = [&](EdgeId x, VertexId v) { return g.edge(x).v[0] == v || g.edge(x).v[1] == v; };
  // next_edge in edge.c: the edge of t, other than e1 and e, on a vertex of e.
  auto next_edge = [&](FaceId t, EdgeId e1) {
    for (EdgeId x : g.face(t).e)
      if (x != e1 && x != e && (uses(x, v1) || uses(x, v2))) return x;
    return EdgeId();
  };
  // triangle_next, depth first as GTS's recursion.
  std::function<void(EdgeId)> next = [&](EdgeId e1) {
    if ( ! e1.valid()) return;
    for (FaceId t : g.faces_of(e1))
      if (take(t)) next(next_edge(t, e1));
  };

  unsigned n = 0;
  for (FaceId t : g.faces_of(e))
    if (take(t))
    {
      const EdgeId e1 = next_edge(t, EdgeId());
      next(e1);
      next(next_edge(t, e1));
      ++n;
    }
  return n;
}

std::vector<EdgeId> edges_from_vertices(const Geometrium& g, std::span<const VertexId> vertices, const Surface* s)
{
  std::vector<EdgeId> list;
  std::vector<bool> seen(g.edge_capacity());
  for (VertexId v : vertices)
    for (EdgeId e : g.edges_of(v))
      if ( ! seen[e.i] && ( ! s || has_face(g, e, *s)))
      {
        seen[e.i] = true;
        list.push_back(e);
      }
  std::reverse(list.begin(), list.end());
  return list;
}

std::vector<FaceId> faces_from_edges(const Geometrium& g, std::span<const EdgeId> edges, const Surface* s)
{
  std::vector<FaceId> list;
  std::vector<bool> seen(g.face_capacity());
  for (EdgeId e : edges)
    for (FaceId t : g.faces_of(e))
      if ( ! seen[t.i] && in(s, t))
      {
        seen[t.i] = true;
        list.push_back(t);
      }
  std::reverse(list.begin(), list.end());
  return list;
}

std::vector<VertexId> vertices_from_edges(const Geometrium& g, std::span<const EdgeId> edges)
{
  std::vector<VertexId> list;
  std::vector<bool> seen(g.vertex_capacity());
  for (EdgeId e : edges)
    for (VertexId v : g.edge(e).v)
      if ( ! seen[v.i])
      {
        seen[v.i] = true;
        list.push_back(v);
      }
  std::reverse(list.begin(), list.end());
  return list;
}

std::size_t remove_faces_if(Surface& s, const std::function<bool(FaceId)>& pred)
{
  std::vector<FaceId> remove;
  for (FaceId f : s.faces())
    if (pred(f)) remove.push_back(f);
  for (FaceId f : remove) s.remove_face(f);
  return remove.size();
}

//==============================================================================
// Traversal and components
//==============================================================================

SurfaceTraverse::SurfaceTraverse(const Surface& s, FaceId start)
  : m_s(s), m_level(s.geometrium().face_capacity(), 0)
{
  m_level[start.i] = 1;
  m_queue.push_back(start);
}

std::optional<std::pair<FaceId, unsigned>> SurfaceTraverse::next()
{
  if (m_queue.empty()) return std::nullopt;
  const FaceId u = m_queue.front();
  m_queue.pop_front();
  const Geometrium& g = m_s.geometrium();
  // gts_face_foreach_neighbor
  for (EdgeId e : g.face(u).e)
    for (FaceId t : g.faces_of(e))
      if (t != u && m_s.contains(t) && ! m_level[t.i])
      {
        m_level[t.i] = m_level[u.i] + 1;
        m_queue.push_back(t);
      }
  return std::pair{ u, m_level[u.i] - 1 };
}

namespace {

// The other face of s on e, when e has exactly two faces of s.
FaceId other_of_two(const Surface& s, EdgeId e, FaceId f)
{
  const auto m = manifold_faces(s, e);
  if ( ! m) return FaceId();
  return (*m)[0] == f ? (*m)[1] : (*m)[0];
}

// traverse_boundary and traverse_remaining: a new component from f, unless
// f is in one already.
void grow_component(Surface& orig, FaceId f, std::vector<bool>& taken, std::vector<Surface*>& components)
{
  if (taken[f.i]) return;
  Geometrium& g = orig.geometrium();
  Surface& s = g.new_surface();
  components.push_back(&s);

  // traverse_manifold, without recursion
  std::vector<FaceId> stack{ f };
  while ( ! stack.empty())
  {
    const FaceId t = stack.back();
    stack.pop_back();
    if (taken[t.i]) continue;
    taken[t.i] = true;
    s.add_face(t);
    const auto E = g.face(t).e;
    for (int k = 2; k >= 0; --k)
      if (const FaceId n = other_of_two(orig, E[k], t); n.valid() && ! taken[n.i])
        stack.push_back(n);
  }

  // non_manifold_edges: the faces of s on an edge with more than two of
  // them leave s, to make components of their own.
  std::vector<FaceId> non_manifold;
  for (EdgeId e : s.edges())
    if (face_number(s, e) > 2)
      for (FaceId t : g.faces_of(e))
        if (s.contains(t)) add_unique(non_manifold, t);
  for (FaceId t : non_manifold)
  {
    taken[t.i] = false;
    s.remove_face(t);
  }
}

} // namespace

std::vector<Surface*> split_components(Surface& s)
{
  const Geometrium& g = s.geometrium();
  std::vector<bool> taken(g.face_capacity());
  std::vector<Surface*> components;
  // Boundary components first.
  for (EdgeId e : s.edges())
    if (const FaceId f = boundary_face(g, e, &s); f.valid())
      grow_component(s, f, taken, components);
  for (FaceId f : s.faces())
    grow_component(s, f, taken, components);
  // A face that left a component after the loop passed it, which GTS lost.
  for (FaceId f : s.faces())
    if ( ! taken[f.i])
    {
      taken[f.i] = true;
      components.push_back(&s.geometrium().new_surface());
      components.back()->add_face(f);
    }
  std::reverse(components.begin(), components.end());
  return components;
}

} // namespace gts

//==============================================================================
// Contact vertices
//==============================================================================

namespace gts {

unsigned Geometrium::vertex_contacts(VertexId v, bool sever)
{
  // The faces around v, marked until reached.
  const std::vector<FaceId> faces = vertex_faces(*this, v);
  std::vector<bool> marked(m_faces.size());
  for (FaceId t : faces) marked[t.i] = true;

  // replace_vertex in vertex.c: the edge of t other than e1 on v, moved to
  // with.
  auto replace = [&](FaceId t, EdgeId e1, VertexId with) {
    EdgeId e;
    for (EdgeId x : m_faces[t.i].e)
      if (x != e1 && (m_edges[x.i].v[0] == v || m_edges[x.i].v[1] == v)) { e = x; break; }
    if ( ! e.valid()) return e;
    if (with != v)
    {
      const int k = vertex_slot(m_edges[e.i], v);
      unlink_edge_slot(e, k);
      m_edges[e.i].v[k] = with;
      link_edge_slot(e, k);
    }
    return e;
  };
  // triangle_next, depth first as GTS's recursion.
  std::function<void(EdgeId, VertexId)> next = [&](EdgeId e, VertexId with) {
    if ( ! e.valid()) return;
    std::vector<FaceId> ring;
    for (FaceId t : faces_of(e)) ring.push_back(t);
    for (FaceId t : ring)
      if (marked[t.i])
      {
        marked[t.i] = false;
        next(replace(t, e, with), with);
      }
  };

  unsigned n = 0;
  for (FaceId t : faces)
    if (marked[t.i])
    {
      VertexId with = v;
      if (n && sever)
      {
        with = add_vertex(m_vertices[v.i].p);
        split_attributes(with, v);
      }
      marked[t.i] = false;
      const EdgeId e = replace(t, EdgeId(), with);
      next(e, with);
      next(replace(t, e, with), with);
      ++n;
    }
  return n;
}

} // namespace gts
