// Copyright (C) 1999 Stéphane Popinet (GTS: triangle.c, edge.c, segment.c, vertex.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/Geometrium.h"

#include "gts/Predicates.h"

// GCC 15 warns about mismatched new and delete inside libstdc++'s
// std::generator promise when it inlines the coroutine frame allocation.
// The pairing is the library's own, not this file's.
#if defined(__GNUC__) && ! defined(__clang__)
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif

#include <algorithm>
#include <iterator>
#include <numbers>
#include <stdexcept>

namespace gts {

//==============================================================================
// Surface
//==============================================================================

void Surface::add_face(FaceId f)
{
  if (m_geo->dead(f)) throw std::invalid_argument("gts::Surface::add_face: the face is dead");
  if (m_faces.insert(f).second)
  {
    Face& F = m_geo->m_faces[f.i];
    if (F.surfaces == UINT16_MAX)
      throw std::length_error("gts::Surface::add_face: face is in too many surfaces");
    ++F.surfaces;
    Geometrium::call(m_face_callouts, FaceEvent{ f, true });
  }
}

void Surface::remove_face(FaceId f)
{
  if (m_faces.erase(f))
  {
    --m_geo->m_faces[f.i].surfaces;
    Geometrium::call(m_face_callouts, FaceEvent{ f, false });
    const Face& F = m_geo->m_faces[f.i];
    if (F.surfaces == 0 && ! (F.flags & kDead))
      m_geo->release_face(f);
  }
}

Subscription Surface::on_face(FaceCallout f)
{
  return Geometrium::subscribe(m_face_callouts, std::move(f));
}

std::generator<EdgeId> Surface::edges() const
{
  std::vector<bool> seen(m_geo->edge_capacity());
  for (FaceId f : m_faces)
    for (EdgeId e : m_geo->face(f).e)
      if ( ! seen[e.i]) { seen[e.i] = true; co_yield e; }
}

std::generator<VertexId> Surface::vertices() const
{
  std::vector<bool> seen(m_geo->vertex_capacity());
  for (FaceId f : m_faces)
  {
    const Edge& e1 = m_geo->edge(m_geo->face(f).e[0]);
    for (VertexId v : { e1.v[0], e1.v[1], m_geo->opposite_vertex(f) })
      if ( ! seen[v.i]) { seen[v.i] = true; co_yield v; }
  }
}

void Surface::merge(const Surface& with)
{
  if (&with.geometrium() != m_geo)
    throw std::invalid_argument("gts::Surface::merge: the surfaces are in different Geometria");
  for (FaceId f : with.m_faces) add_face(f);
}

std::size_t Surface::edge_number() const
{
  std::size_t n = 0;
  for (EdgeId e : edges()) { (void) e; ++n; }
  return n;
}

std::size_t Surface::vertex_number() const
{
  std::size_t n = 0;
  for (VertexId v : vertices()) { (void) v; ++n; }
  return n;
}

//==============================================================================
// Elements and rings
//==============================================================================

VertexId Geometrium::add_vertex(Vec3 p)
{
  Vertex v;
  v.p = p;
  m_vertices.push_back(v);
  for (auto& c : m_vertex_columns) c.resize(m_vertices.size());
  return VertexId(Index(m_vertices.size() - 1));
}

VertexId Geometrium::interpolate_vertex(VertexId a, VertexId b, double t)
{
  const Vec3 pa = m_vertices[a.i].p, pb = m_vertices[b.i].p;
  return add_interpolated_vertex(pa + t * (pb - pa), a, b, t);
}

VertexId Geometrium::add_merged_vertex(Vec3 p, VertexId a, VertexId b)
{
  const VertexId v = add_vertex(p);
  for (auto& c : m_vertex_columns) c.merge(v.i, a.i, b.i);
  return v;
}

VertexId Geometrium::add_interpolated_vertex(Vec3 p, VertexId a, VertexId b, double t)
{
  const VertexId v = add_vertex(p);
  for (auto& c : m_vertex_columns) c.interpolate(v.i, a.i, b.i, t);
  return v;
}

// Inserts e after the last edge of the ring, which makes it the first.
void Geometrium::link_edge_slot(EdgeId e, int k)
{
  Vertex& V = m_vertices[m_edges[e.i].v[k].i];
  if ( ! V.edge.valid()) { V.edge = e; m_edges[e.i].next[k] = e; return; }
  Edge& last = m_edges[V.edge.i];
  const int s = vertex_slot(last, m_edges[e.i].v[k]);
  m_edges[e.i].next[k] = last.next[s];
  last.next[s] = e;
}

void Geometrium::link_edge(EdgeId e)
{
  for (int k = 0; k < 2; ++k) link_edge_slot(e, k);
}

EdgeId Geometrium::add_edge(VertexId a, VertexId b)
{
  if (a == b)
    throw std::invalid_argument("gts::Geometrium::add_edge: both ends are the same vertex");
  Edge x;
  x.v = { a, b };
  m_edges.push_back(x);
  for (auto& c : m_edge_columns) c.resize(m_edges.size());
  const EdgeId e(Index(m_edges.size() - 1));
  link_edge(e);
  return e;
}

// Inserts f after the last face of the ring, which makes it the first.
void Geometrium::link_face_slot(FaceId f, int k)
{
  Edge& E = m_edges[m_faces[f.i].e[k].i];
  if ( ! E.face.valid()) { E.face = f; m_faces[f.i].next[k] = f; return; }
  Face& last = m_faces[E.face.i];
  const int s = edge_slot(last, m_faces[f.i].e[k]);
  m_faces[f.i].next[k] = last.next[s];
  last.next[s] = f;
}

void Geometrium::link_face(FaceId f)
{
  for (int k = 0; k < 3; ++k) link_face_slot(f, k);
}

// Removes f from the ring, keeping the order of the others.
void Geometrium::unlink_face_slot(FaceId f, int k)
{
  const EdgeId e = m_faces[f.i].e[k];
  Edge& E = m_edges[e.i];
  const FaceId next = m_faces[f.i].next[k];
  if (next == f) { E.face = FaceId(); return; }
  FaceId p = E.face; // the predecessor of f in the ring around e
  while (m_faces[p.i].next[edge_slot(m_faces[p.i], e)] != f)
    p = m_faces[p.i].next[edge_slot(m_faces[p.i], e)];
  m_faces[p.i].next[edge_slot(m_faces[p.i], e)] = next;
  if (E.face == f) E.face = p; // the predecessor becomes the last
}

bool Geometrium::forms_triangle(EdgeId a, EdgeId b, EdgeId c) const
{
  VertexId v[6] = { m_edges[a.i].v[0], m_edges[a.i].v[1], m_edges[b.i].v[0],
                    m_edges[b.i].v[1], m_edges[c.i].v[0], m_edges[c.i].v[1] };
  std::sort(std::begin(v), std::end(v));
  return a != b && b != c && a != c &&
         v[0] == v[1] && v[2] == v[3] && v[4] == v[5] && v[1] != v[2] && v[3] != v[4];
}

FaceId Geometrium::add_face(EdgeId a, EdgeId b, EdgeId c)
{
  if ( ! forms_triangle(a, b, c))
    throw std::invalid_argument("gts::Geometrium::add_face: the edges do not form a triangle");
  Face x;
  x.e = { a, b, c };
  m_faces.push_back(x);
  for (auto& col : m_face_columns) col.resize(m_faces.size());
  const FaceId f(Index(m_faces.size() - 1));
  link_face(f);
  return f;
}

void Geometrium::set_face_edges(FaceId f, EdgeId a, EdgeId b, EdgeId c)
{
  if ( ! forms_triangle(a, b, c))
    throw std::invalid_argument("gts::Geometrium::set_face_edges: the edges do not form a triangle");
  relink_face(f, { a, b, c });
}

void Geometrium::relink_face(FaceId f, const std::array<EdgeId, 3>& ne)
{
  auto find = [](const std::array<EdgeId, 3>& es, EdgeId e) {
    for (int k = 0; k < 3; ++k) if (es[k] == e) return k;
    return -1;
  };
  for (int k = 0; k < 3; ++k)
    if (find(ne, m_faces[f.i].e[k]) < 0) unlink_face_slot(f, k);

  Face& F = m_faces[f.i];
  const auto old_e = F.e;
  const auto old_next = F.next;
  F.e = ne;
  for (int k = 0; k < 3; ++k)
  {
    const int j = find(old_e, ne[k]);
    if (j >= 0) F.next[k] = old_next[j];
  }
  for (int k = 0; k < 3; ++k)
    if (find(old_e, ne[k]) < 0) link_face_slot(f, k);
}

void Geometrium::revert_face(FaceId f)
{
  Face& F = m_faces[f.i];
  std::swap(F.e[0], F.e[1]);
  std::swap(F.next[0], F.next[1]);
}

void Geometrium::kill_face(FaceId f)
{
  if (m_faces[f.i].surfaces != 0)
    throw std::logic_error("gts::Geometrium::kill_face: the face is in a surface");
  for (int k = 0; k < 3; ++k) unlink_face_slot(f, k);
  m_faces[f.i].flags |= kDead;
}

void Geometrium::unlink_edge_slot(EdgeId e, int k)
{
  const VertexId v = m_edges[e.i].v[k];
  Vertex& V = m_vertices[v.i];
  const EdgeId next = m_edges[e.i].next[k];
  if (next == e) { V.edge = EdgeId(); return; }
  EdgeId p = V.edge;
  while (m_edges[p.i].next[vertex_slot(m_edges[p.i], v)] != e)
    p = m_edges[p.i].next[vertex_slot(m_edges[p.i], v)];
  m_edges[p.i].next[vertex_slot(m_edges[p.i], v)] = next;
  if (V.edge == e) V.edge = p;
}

void Geometrium::kill_edge(EdgeId e)
{
  if (m_edges[e.i].face.valid())
    throw std::logic_error("gts::Geometrium::kill_edge: the edge has faces");
  for (int k = 0; k < 2; ++k) unlink_edge_slot(e, k);
  m_edges[e.i].flags |= kDead;
}

void Geometrium::park(FaceId f)
{
  if (dead(f)) throw std::invalid_argument("gts::Geometrium::park: the face is dead");
  destroy_face(f);
  m_faces[f.i].flags |= kParked;
  m_faces[f.i].e = { EdgeId(), EdgeId(), EdgeId() };
}

void Geometrium::unpark(FaceId f, EdgeId a, EdgeId b, EdgeId c)
{
  if ( ! parked(f)) throw std::invalid_argument("gts::Geometrium::unpark: the face is not parked");
  if ( ! forms_triangle(a, b, c)) throw std::invalid_argument("gts::Geometrium::unpark: the edges do not form a triangle");
  Face& F = m_faces[f.i];
  F.flags &= ~(kDead | kParked);
  F.e = { a, b, c };
  F.surfaces = 0;
  link_face(f);
}

void Geometrium::park(VertexId v)
{
  if (dead(v) || m_vertices[v.i].edge.valid())
    throw std::invalid_argument("gts::Geometrium::park: the vertex is dead or has edges");
  m_vertices[v.i].flags |= kDead | kParked;
}

void Geometrium::unpark(VertexId v)
{
  if ( ! parked(v)) throw std::invalid_argument("gts::Geometrium::unpark: the vertex is not parked");
  m_vertices[v.i].flags &= ~(kDead | kParked);
}

void Geometrium::reattach(EdgeId e, VertexId from, VertexId to)
{
  Edge& E = m_edges[e.i];
  if ((E.v[0] != from && E.v[1] != from) || E.v[0] == to || E.v[1] == to)
    throw std::invalid_argument("gts::Geometrium::reattach: from is not an end of the edge, or to is");
  const int k = vertex_slot(E, from);
  unlink_edge_slot(e, k);
  E.v[k] = to;
  link_edge_slot(e, k);
}

void Geometrium::kill_vertex(VertexId v)
{
  if (m_vertices[v.i].edge.valid())
    throw std::logic_error("gts::Geometrium::kill_vertex: the vertex has edges");
  m_vertices[v.i].flags |= kDead;
}

// As gts_object_destroy on a face: face_destroy, then triangle_destroy for
// edges left without faces, other than keep_e, and segment_destroy for
// vertices left without edges, other than keep_v.
void Geometrium::destroy_cascade(FaceId f, EdgeId keep_e, VertexId keep_v)
{
  const auto edges = m_faces[f.i].e;
  destroy_face(f);
  if (m_orphans.edges) return;
  for (EdgeId e : edges)
  {
    if (e == keep_e || dead(e) || m_edges[e.i].face.valid()) continue;
    const auto vs = m_edges[e.i].v;
    kill_edge(e);
    if (m_orphans.vertices) continue;
    for (VertexId v : vs)
      if (v != keep_v && ! dead(v) && ! m_vertices[v.i].edge.valid())
        m_vertices[v.i].flags |= kDead;
  }
}

// As gts_object_destroy on an edge: edge_destroy, then segment_destroy.
void Geometrium::destroy_cascade(EdgeId e, VertexId keep_v)
{
  std::vector<FaceId> faces;
  for (FaceId f : faces_of(e)) faces.push_back(f);
  for (FaceId f : faces) destroy_cascade(f, e, keep_v);
  const auto vs = m_edges[e.i].v;
  kill_edge(e);
  if (m_orphans.vertices) return;
  for (VertexId v : vs)
    if (v != keep_v && ! dead(v) && ! m_vertices[v.i].edge.valid())
      m_vertices[v.i].flags |= kDead;
}

void Geometrium::destroy(FaceId f)
{
  if ( ! dead(f)) destroy_cascade(f, EdgeId(), VertexId());
}

void Geometrium::destroy(EdgeId e)
{
  if ( ! dead(e)) destroy_cascade(e, VertexId());
}

void Geometrium::destroy(VertexId v)
{
  if (dead(v)) return;
  std::vector<EdgeId> ring;
  for (EdgeId e : edges_of(v)) ring.push_back(e);
  for (EdgeId e : ring)
    if ( ! dead(e)) destroy_cascade(e, v);
  m_vertices[v.i].flags |= kDead;
}

// As GTS's destruction cascade: face_destroy, then triangle_destroy for edges
// left without triangles, then segment_destroy for vertices left without
// segments, each unless the matching gts_allow_floating_* was set.
void Geometrium::release_face(FaceId f)
{
  if (m_orphans.faces) return;
  const auto edges = m_faces[f.i].e;
  kill_face(f);
  if (m_orphans.edges) return;
  for (EdgeId e : edges)
  {
    if (dead(e) || m_edges[e.i].face.valid()) continue;
    const auto vs = m_edges[e.i].v;
    kill_edge(e);
    if (m_orphans.vertices) continue;
    for (VertexId v : vs)
      if ( ! dead(v) && ! m_vertices[v.i].edge.valid())
        m_vertices[v.i].flags |= kDead;
  }
}

//==============================================================================
// Adjacency queries
//==============================================================================

std::array<VertexId, 3> Geometrium::vertices_of(FaceId f) const
{
  const Edge& e1 = m_edges[m_faces[f.i].e[0].i];
  const Edge& e2 = m_edges[m_faces[f.i].e[1].i];
  if (e1.v[1] == e2.v[0]) return { e1.v[0], e1.v[1], e2.v[1] };
  if (e1.v[1] == e2.v[1]) return { e1.v[0], e1.v[1], e2.v[0] };
  if (e1.v[0] == e2.v[0]) return { e1.v[1], e1.v[0], e2.v[1] };
  return { e1.v[1], e1.v[0], e2.v[0] };
}

VertexId Geometrium::opposite_vertex(FaceId f) const
{
  const Edge& e1 = m_edges[m_faces[f.i].e[0].i];
  const Edge& e2 = m_edges[m_faces[f.i].e[1].i];
  return (e1.v[0] == e2.v[0] || e1.v[1] == e2.v[0]) ? e2.v[1] : e2.v[0];
}

VertexId Geometrium::opposite_vertex(FaceId f, EdgeId e) const
{
  const Edge& E = m_edges[e.i];
  for (VertexId v : vertices_of(f))
    if (v != E.v[0] && v != E.v[1]) return v;
  return VertexId();
}

Geometrium::VerticesEdges Geometrium::vertices_edges(FaceId f, EdgeId e) const
{
  const auto& E = m_faces[f.i].e;
  VerticesEdges r;
  if      ( ! e.valid() || e == E[0]) r.e = { E[0], E[1], E[2] };
  else if (e == E[1])                 r.e = { E[1], E[2], E[0] };
  else if (e == E[2])                 r.e = { E[2], E[0], E[1] };
  else throw std::invalid_argument("gts::Geometrium::vertices_edges: the edge is not an edge of the face");
  const Edge& a = m_edges[r.e[0].i];
  const Edge& b = m_edges[r.e[1].i];
  if      (a.v[1] == b.v[0]) r.v = { a.v[0], a.v[1], b.v[1] };
  else if (a.v[1] == b.v[1]) r.v = { a.v[0], a.v[1], b.v[0] };
  else if (a.v[0] == b.v[0]) r.v = { a.v[1], a.v[0], b.v[1] };
  else                       r.v = { a.v[1], a.v[0], b.v[0] };
  return r;
}

EdgeId Geometrium::opposite_edge(FaceId f, VertexId v) const
{
  const auto& E = m_faces[f.i].e;
  auto touches = [&](EdgeId e) { return m_edges[e.i].v[0] == v || m_edges[e.i].v[1] == v; };
  if ( ! touches(E[0])) return touches(E[1]) ? E[0] : EdgeId();
  if ( ! touches(E[1])) return E[1];
  return E[2];
}

EdgeId Geometrium::common_edge(FaceId f1, FaceId f2) const
{
  const auto& A = m_faces[f1.i].e;
  const auto& B = m_faces[f2.i].e;
  for (EdgeId e : A)
    if (e == B[0] || e == B[1] || e == B[2]) return e;
  return EdgeId();
}

EdgeId Geometrium::find_edge(VertexId a, VertexId b) const
{
  for (EdgeId e : edges_of(a))
    if (other_vertex(e, a) == b) return e;
  return EdgeId();
}

EdgeId Geometrium::duplicate_of(EdgeId e) const
{
  const VertexId v1 = m_edges[e.i].v[0], v2 = m_edges[e.i].v[1];
  for (EdgeId s : edges_of(v1))
    if (s != e && other_vertex(s, v1) == v2) return s;
  return EdgeId();
}

FaceId Geometrium::duplicate_of(FaceId f) const
{
  const EdgeId e2 = m_faces[f.i].e[1], e3 = m_faces[f.i].e[2];
  auto has = [&](FaceId t, EdgeId e) { const auto& E = m_faces[t.i].e; return E[0] == e || E[1] == e || E[2] == e; };
  for (FaceId t : faces_of(m_faces[f.i].e[0]))
    if (t != f && has(t, e2) && has(t, e3)) return t;
  return FaceId();
}

bool Geometrium::compatible(FaceId f1, FaceId f2, EdgeId e) const
{
  auto after = [&](FaceId f) {
    const auto& E = m_faces[f.i].e;
    return E[0] == e ? E[1] : E[1] == e ? E[2] : E[0];
  };
  const Edge& a = m_edges[after(f1).i];
  const Edge& b = m_edges[after(f2).i];
  return ! (a.v[0] == b.v[0] || a.v[0] == b.v[1] || a.v[1] == b.v[0] || a.v[1] == b.v[1]);
}

//==============================================================================
// Geometry of elements
//==============================================================================

Triangle Geometrium::triangle(FaceId f) const
{
  const auto v = vertices_of(f);
  return { m_vertices[v[0].i].p, m_vertices[v[1].i].p, m_vertices[v[2].i].p };
}

// gts_triangle_normal chooses its vertices by its own rule, which gives the
// same orientation as gts_triangle_vertices but not always the same first
// vertex. It is kept, so that normals and areas agree with GTS to the bit.
Vec3 Geometrium::normal(FaceId f) const
{
  const Edge& e1 = m_edges[m_faces[f.i].e[0].i];
  const Edge& e2 = m_edges[m_faces[f.i].e[1].i];
  const VertexId v1 = e1.v[0];
  VertexId v2, v3;
  if      (e1.v[0] == e2.v[0]) { v2 = e2.v[1]; v3 = e1.v[1]; }
  else if (e1.v[1] == e2.v[1]) { v2 = e1.v[1]; v3 = e2.v[0]; }
  else if (e1.v[0] == e2.v[1]) { v2 = e2.v[0]; v3 = e1.v[1]; }
  else if (e1.v[1] == e2.v[0]) { v2 = e1.v[1]; v3 = e2.v[1]; }
  else throw std::logic_error("gts::Geometrium::normal: the first two edges of the face share no vertex");
  const Vec3 p1 = m_vertices[v1.i].p;
  return cross(m_vertices[v2.i].p - p1, m_vertices[v3.i].p - p1);
}

// gts_triangle_orientation chooses its vertices as gts_triangle_normal does.
double Geometrium::orientation(FaceId f) const
{
  const Edge& e1 = m_edges[m_faces[f.i].e[0].i];
  const Edge& e2 = m_edges[m_faces[f.i].e[1].i];
  const VertexId v1 = e1.v[0];
  VertexId v2, v3;
  if      (e1.v[0] == e2.v[0]) { v2 = e2.v[1]; v3 = e1.v[1]; }
  else if (e1.v[1] == e2.v[1]) { v2 = e1.v[1]; v3 = e2.v[0]; }
  else if (e1.v[0] == e2.v[1]) { v2 = e2.v[0]; v3 = e1.v[1]; }
  else                         { v2 = e1.v[1]; v3 = e2.v[1]; }
  const Vec3 a = m_vertices[v1.i].p, b = m_vertices[v2.i].p, c = m_vertices[v3.i].p;
  const double pa[2] = { a.x, a.y }, pb[2] = { b.x, b.y }, pc[2] = { c.x, c.y };
  return orient2d(pa, pb, pc);
}

double Geometrium::perimeter(FaceId f) const
{
  const Edge& e1 = m_edges[m_faces[f.i].e[0].i];
  const Vec3 a = m_vertices[e1.v[0].i].p, b = m_vertices[e1.v[1].i].p;
  const Vec3 c = m_vertices[opposite_vertex(f).i].p;
  return distance(a, b) + distance(a, c) + distance(b, c);
}

double Geometrium::quality(FaceId f) const
{
  const double p = perimeter(f);
  return p > 0.0 ? kGoldenPerimeter * std::sqrt(area(f)) / p : 0.0;
}

double Geometrium::angle(FaceId f1, FaceId f2) const
{
  const Vec3 n1 = normal(f1), n2 = normal(f2);
  const double theta = std::atan2(norm(cross(n1, n2)), dot(n1, n2)) - std::numbers::pi;
  return theta < -std::numbers::pi ? theta + 2. * std::numbers::pi : theta;
}

//==============================================================================
// Surfaces
//==============================================================================

Surface& Geometrium::new_surface(std::string name)
{
  m_surfaces.push_back(std::unique_ptr<Surface>(new Surface(*this, std::move(name))));
  return *m_surfaces.back();
}

void Geometrium::remove_surface(Surface& s)
{
  const std::vector<FaceId> faces(s.m_faces.begin(), s.m_faces.end());
  for (FaceId f : faces) s.remove_face(f);
  std::erase_if(m_surfaces, [&](const auto& p) { return p.get() == &s; });
}

Surface& Geometrium::import(const Geometrium& g, const Surface& src)
{
  Surface& dst = new_surface(src.name());
  import(g, src, dst);
  return dst;
}

void Geometrium::import(const Geometrium& g, const Surface& src, Surface& dst)
{
  if (&dst == &src)
    throw std::invalid_argument("gts::Geometrium::import: copying a surface into itself");
  if (&dst.geometrium() != this || &src.geometrium() != &g)
    throw std::invalid_argument("gts::Geometrium::import: a surface of another Geometrium");
  std::vector<Index> map_v(g.vertex_capacity(), kNone), map_e(g.edge_capacity(), kNone);

  // Attributes present on both sides under the same name and type.
  using Pairs = std::vector<std::pair<detail::Column*, const detail::Column*>>;
  auto match = [](std::vector<detail::Column>& mine, const std::vector<detail::Column>& theirs) {
    Pairs p;
    for (auto& c : mine)
      for (const auto& s : theirs)
        if (c.name == s.name && c.type == s.type) p.emplace_back(&c, &s);
    return p;
  };
  const Pairs vcols = match(m_vertex_columns, g.m_vertex_columns);
  const Pairs ecols = match(m_edge_columns,   g.m_edge_columns);
  const Pairs fcols = match(m_face_columns,   g.m_face_columns);

  for (FaceId f : src.faces())
  {
    std::array<EdgeId, 3> es;
    for (int k = 0; k < 3; ++k)
    {
      const EdgeId se = g.face(f).e[k];
      if (map_e[se.i] == kNone)
      {
        std::array<VertexId, 2> vs;
        for (int j = 0; j < 2; ++j)
        {
          const VertexId sv = g.edge(se).v[j];
          if (map_v[sv.i] == kNone)
          {
            map_v[sv.i] = add_vertex(g.vertex(sv).p).i;
            m_vertices[map_v[sv.i]].flags = g.vertex(sv).flags & ~kDead;
            for (auto [c, s] : vcols) c->copy_from(map_v[sv.i], *s, sv.i);
          }
          vs[j] = VertexId(map_v[sv.i]);
        }
        map_e[se.i] = add_edge(vs[0], vs[1]).i;
        m_edges[map_e[se.i]].flags = g.edge(se).flags & ~kDead;
        for (auto [c, s] : ecols) c->copy_from(map_e[se.i], *s, se.i);
      }
      es[k] = EdgeId(map_e[se.i]);
    }
    const FaceId nf = add_face(es[0], es[1], es[2]);
    m_faces[nf.i].flags = g.face(f).flags & ~kDead;
    for (auto [c, s] : fcols) c->copy_from(nf.i, *s, f.i);
    dst.add_face(nf);
  }
}

//==============================================================================
// Collection
//==============================================================================

void Geometrium::collect(CollectOptions o)
{
  const Index nv0 = Index(m_vertices.size()), ne0 = Index(m_edges.size()), nf0 = Index(m_faces.size());

  // Which elements stay.
  std::vector<bool> keep_f(nf0), keep_e(ne0), keep_v(nv0);
  for (const auto& s : m_surfaces)
    for (FaceId f : s->m_faces) keep_f[f.i] = true;
  for (Index i = 0; i < nf0; ++i)
    keep_f[i] = (m_faces[i].flags & kParked) || ( ! (m_faces[i].flags & kDead) && (keep_f[i] || ! o.faces));
  for (Index i = 0; i < nf0; ++i)
    if (keep_f[i] && ! (m_faces[i].flags & kParked)) for (EdgeId e : m_faces[i].e) keep_e[e.i] = true;
  for (Index i = 0; i < ne0; ++i)
    keep_e[i] = ! (m_edges[i].flags & kDead) && (keep_e[i] || ! o.edges);
  for (Index i = 0; i < ne0; ++i)
    if (keep_e[i]) for (VertexId v : m_edges[i].v) keep_v[v.i] = true;
  for (Index i = 0; i < nv0; ++i)
    keep_v[i] = (m_vertices[i].flags & kParked) || ( ! (m_vertices[i].flags & kDead) && (keep_v[i] || ! o.vertices));

  // New numbers.
  Renumbering r;
  auto number = [](const std::vector<bool>& keep, std::vector<Index>& map) {
    map.assign(keep.size(), kNone);
    Index n = 0;
    for (std::size_t i = 0; i < keep.size(); ++i) if (keep[i]) map[i] = n++;
    return n;
  };
  const Index nv = number(keep_v, r.vertices), ne = number(keep_e, r.edges), nf = number(keep_f, r.faces);
  const auto& map_v = r.vertices;
  const auto& map_e = r.edges;
  const auto& map_f = r.faces;

  // The new elements, with empty rings.
  std::vector<Vertex> V(nv);
  std::vector<Edge>   E(ne);
  std::vector<Face>   F(nf);
  for (Index i = 0; i < nv0; ++i)
    if (map_v[i] != kNone) { Vertex& x = V[map_v[i]]; x.p = m_vertices[i].p; x.flags = m_vertices[i].flags; }
  for (Index i = 0; i < ne0; ++i)
    if (map_e[i] != kNone)
    {
      Edge& x = E[map_e[i]];
      for (int k = 0; k < 2; ++k) x.v[k] = VertexId(map_v[m_edges[i].v[k].i]);
      x.flags = m_edges[i].flags;
    }
  for (Index i = 0; i < nf0; ++i)
    if (map_f[i] != kNone)
    {
      Face& x = F[map_f[i]];
      for (int k = 0; k < 3; ++k)
        if (m_faces[i].e[k].valid()) x.e[k] = EdgeId(map_e[m_faces[i].e[k].i]);
      x.flags = m_faces[i].flags;
    }

  // Every ring in its old order, without the removed elements.
  std::vector<EdgeId> ring_e;
  for (Index i = 0; i < nv0; ++i)
  {
    if (map_v[i] == kNone) continue;
    const VertexId nvid(map_v[i]);
    ring_e.clear();
    for (EdgeId e : edges_of(VertexId(i))) if (map_e[e.i] != kNone) ring_e.push_back(EdgeId(map_e[e.i]));
    for (std::size_t j = 0; j < ring_e.size(); ++j)
      E[ring_e[j].i].next[vertex_slot(E[ring_e[j].i], nvid)] = ring_e[(j + 1) % ring_e.size()];
    if ( ! ring_e.empty()) V[nvid.i].edge = ring_e.back();
  }
  std::vector<FaceId> ring_f;
  for (Index i = 0; i < ne0; ++i)
  {
    if (map_e[i] == kNone) continue;
    const EdgeId neid(map_e[i]);
    ring_f.clear();
    for (FaceId f : faces_of(EdgeId(i))) if (map_f[f.i] != kNone) ring_f.push_back(FaceId(map_f[f.i]));
    for (std::size_t j = 0; j < ring_f.size(); ++j)
      F[ring_f[j].i].next[edge_slot(F[ring_f[j].i], neid)] = ring_f[(j + 1) % ring_f.size()];
    if ( ! ring_f.empty()) E[neid.i].face = ring_f.back();
  }

  m_vertices.swap(V);
  m_edges.swap(E);
  m_faces.swap(F);
  for (auto& c : m_vertex_columns) c.renumber(map_v, nv);
  for (auto& c : m_edge_columns)   c.renumber(map_e, ne);
  for (auto& c : m_face_columns)   c.renumber(map_f, nf);

  for (auto& s : m_surfaces)
  {
    std::flat_set<FaceId> faces;
    for (FaceId f : s->m_faces)
      if (map_f[f.i] != kNone)
      {
        faces.insert(FaceId(map_f[f.i]));
        ++m_faces[map_f[f.i]].surfaces;
      }
    s->m_faces.swap(faces);
  }

  call(m_renumber_callouts, r);
}

void Geometrium::truncate(std::size_t nv, std::size_t ne, std::size_t nf)
{
  m_vertices.resize(nv);
  m_edges.resize(ne);
  m_faces.resize(nf);
  for (auto& c : m_vertex_columns) c.resize(nv);
  for (auto& c : m_edge_columns)   c.resize(ne);
  for (auto& c : m_face_columns)   c.resize(nf);
}

//==============================================================================
// Edge swap
//==============================================================================

void swap_edge(Surface& s, EdgeId e)
{
  Geometrium& g = s.geometrium();
  FaceId t1, t2;
  for (FaceId t : g.faces_of(e))
  {
    if ( ! s.contains(t)) continue;
    if      ( ! t1.valid()) t1 = t;
    else if ( ! t2.valid()) t2 = t;
    else throw std::invalid_argument("gts::swap_edge: the edge has more than two faces in the surface");
  }
  if ( ! t2.valid())
    throw std::invalid_argument("gts::swap_edge: the edge has fewer than two faces in the surface");

  // GTS's v1 v2 v3 and e e1 e2 of t1, v4 v5 v6 and e e3 e4 of t2.
  const Geometrium::VerticesEdges a = g.vertices_edges(t1, e);
  const Geometrium::VerticesEdges b = g.vertices_edges(t2, e);
  if (a.v[1] != b.v[0] || a.v[0] != b.v[1])
    throw std::invalid_argument("gts::swap_edge: the faces are not oriented consistently");

  EdgeId v3v6 = g.find_edge(a.v[2], b.v[2]);
  if ( ! v3v6.valid()) v3v6 = g.add_edge(a.v[2], b.v[2]);
  auto add = [&](EdgeId x, EdgeId y, EdgeId z, FaceId parent) {
    FaceId f = g.add_face(x, y, z);
    if (const FaceId d = g.duplicate_of(f); d.valid()) { g.kill_face(f); f = d; }
    else g.split_attributes(f, parent);
    s.add_face(f);
  };
  add(a.e[1], v3v6, b.e[2], t1);
  add(v3v6, a.e[2], b.e[1], t2);
  s.remove_face(t1);
  s.remove_face(t2);
}

} // namespace gts
