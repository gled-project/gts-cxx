// Geometrium: where vertices, edges, faces and surfaces live.
//
// A Geometrium owns the element arrays and their attributes. Surfaces are sets
// of faces inside it; a face can belong to several surfaces of the same
// Geometrium. When a face leaves its last surface it is unlinked at once, and
// so is any edge left without faces, and any vertex left without edges, as GTS
// destroys them; KeepOrphans keeps them, as GTS's gts_allow_floating_* flags
// did. Unlinked elements are marked dead, and collect() reclaims their storage
// and renumbers everything, telling subscribers how.
//
// A Geometrium has no locks. One thread may write to it at a time, or any
// number of threads may read it while nobody writes. Separate Geometria are
// independent; import() copies a surface from one to another.
//
// split_edge() and split_face4() are the primitives of refinement and
// tessellation, collapse_edge() that of coarsening. They set the attributes
// of what they make with the hooks, and report to the on_split() and
// on_collapse() subscribers.
//
// Copyright (C) 1999 Stéphane Popinet (GTS: triangle.c, edge.c, segment.c, vertex.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "gts/Attribute.h"
#include "gts/Elements.h"
#include "gts/Geometry.h"

#include <array>
#include <cstdint>
#include <flat_set>
#include <functional>
#include <generator>
#include <iterator>
#include <memory>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace gts {

class Geometrium;

//==============================================================================
// Surface
//==============================================================================

class Surface
{
public:
  Surface(const Surface&) = delete;
  Surface& operator=(const Surface&) = delete;

  Geometrium&       geometrium()       { return *m_geo; }
  const Geometrium& geometrium() const { return *m_geo; }

  const std::string& name() const { return m_name; }
  void set_name(std::string n)    { m_name = std::move(n); }

  // Adding a face that is already there, or removing one that is not, does
  // nothing. A face removed from its last surface is unlinked; see
  // Geometrium.
  void add_face(FaceId f);
  void remove_face(FaceId f);
  bool contains(FaceId f) const { return m_faces.contains(f); }
  // gts_surface_merge: adds the faces of with, which must be a surface of
  // the same Geometrium.
  void merge(const Surface& with);

  // Faces in index order.
  const std::flat_set<FaceId>& faces() const { return m_faces; }

  // Each edge and vertex once, in the order GTS's foreach functions visit
  // them for the same order of faces: e1, e2, e3 of each face, and for
  // vertices the two of e1 followed by the opposite one.
  std::generator<EdgeId>   edges() const;
  std::generator<VertexId> vertices() const;

  std::size_t face_number() const { return m_faces.size(); }
  std::size_t edge_number() const;
  std::size_t vertex_number() const;

private:
  friend class Geometrium;
  Surface(Geometrium& g, std::string n) : m_geo(&g), m_name(std::move(n)) {}

  Geometrium*           m_geo;
  std::string           m_name;
  std::flat_set<FaceId> m_faces;
};

//==============================================================================
// Renumbering
//==============================================================================

// The new numbers after collect(): old -> new. A removed element maps to an
// invalid index.
struct Renumbering
{
  std::vector<Index> vertices, edges, faces;

  VertexId operator()(VertexId v) const { return VertexId(vertices[v.i]); }
  EdgeId   operator()(EdgeId e) const   { return EdgeId(edges[e.i]); }
  FaceId   operator()(FaceId f) const   { return FaceId(faces[f.i]); }
};

using RenumberCallout = std::function<void(const Renumbering&)>;

//==============================================================================
// Splits
//==============================================================================

// Makes the vertex that splits an edge, usually with
// Geometrium::add_interpolated_vertex() so that vertex attributes are
// interpolated.
using Splitter = std::function<VertexId(Geometrium&, EdgeId)>;

// An edge split at a new vertex into two halves, the one from the edge's
// first vertex first. The parent may be dead by the time it is reported;
// its attributes can still be read.
struct EdgeSplitEvent
{
  EdgeId                parent;
  VertexId              vertex;
  std::array<EdgeId, 2> halves;
};

// A face split into pieces by new edges, the cuts. The parent keeps its
// index and is one of the pieces.
struct FaceSplitEvent
{
  FaceId                     parent;
  std::span<const FaceId>    pieces;
  std::span<const EdgeId>    cuts;
};

// What one split_edge() or split_face4() made.
struct SplitEvent
{
  std::span<const EdgeSplitEvent> edges;
  std::span<const FaceSplitEvent> faces;
};

using SplitCallout = std::function<void(const SplitEvent&)>;

// The halves of edges that split_face4() has split, so that the faces on
// both sides of an edge share its halves. Valid until collect().
class SplitMemo
{
private:
  friend class Geometrium;
  std::vector<std::array<EdgeId, 2>> m_halves; // by edge index, in the order they were made
};

//==============================================================================
// Collapses
//==============================================================================

// An edge collapsed into a new vertex. The faces on the edge are removed.
// Each pair of edges that the collapse makes duplicate becomes one: the kept
// edge takes the faces of the other. An edge left without faces is removed,
// and so is the new vertex when it is left without edges; vertex is then
// invalid. Removed elements are dead; their attributes can still be read.
struct CollapseEvent
{
  EdgeId                                edge;
  std::array<VertexId, 2>               vertices;
  VertexId                              vertex;
  std::span<const FaceId>               faces;
  std::span<const std::array<EdgeId, 2>> merged_edges; // kept, removed
  std::span<const EdgeId>               removed_edges;
};

using CollapseCallout = std::function<void(const CollapseEvent&)>;

//==============================================================================
// Subscriptions
//==============================================================================

namespace detail {
template <class F>
struct CalloutRegistry
{
  std::uint64_t next_id = 0;
  std::vector<std::pair<std::uint64_t, F>> callouts;
};
} // namespace detail

// Keeps a callout registered; cancels it when destroyed or reset. May
// outlive its Geometrium.
class Subscription
{
public:
  Subscription() = default;
  Subscription(Subscription&& o) noexcept : m_cancel(std::exchange(o.m_cancel, nullptr)) {}
  Subscription& operator=(Subscription&& o) noexcept { reset(); m_cancel = std::exchange(o.m_cancel, nullptr); return *this; }
  ~Subscription() { reset(); }
  void reset() { if (auto c = std::exchange(m_cancel, nullptr)) c(); }

private:
  friend class Geometrium;
  explicit Subscription(std::function<void()> cancel) : m_cancel(std::move(cancel)) {}

  std::function<void()> m_cancel;
};

//==============================================================================
// Rings
//==============================================================================

// The elements of one ring, newest first, as a range without allocation:
// for (EdgeId e : g.edges_of(v)). The ring must not change during the walk.
template <class ElemId, class CenterId>
class Ring
{
public:
  struct iterator
  {
    const Geometrium* g;
    CenterId          center;
    ElemId            cur, last;
    bool              done;

    ElemId    operator*() const { return cur; }
    iterator& operator++();
    bool operator==(std::default_sentinel_t) const { return done; }
  };

  Ring(const Geometrium* g, CenterId c, ElemId last) : m_g(g), m_center(c), m_last(last) {}
  iterator begin() const;
  std::default_sentinel_t end() const { return {}; }

private:
  const Geometrium* m_g;
  CenterId          m_center;
  ElemId            m_last;
};

using EdgeRing = Ring<EdgeId, VertexId>; // the edges around a vertex
using FaceRing = Ring<FaceId, EdgeId>;   // the faces around an edge

//==============================================================================
// Geometrium
//==============================================================================

// Which unused elements collect() removes, besides the dead ones: faces in no
// surface, edges in no face, vertices in no edge.
struct CollectOptions
{
  bool faces = true, edges = true, vertices = true;
};

// Which orphans survive when a face leaves its last surface.
struct OrphanPolicy
{
  bool faces = false, edges = false, vertices = false;
};

class Geometrium
{
public:
  Geometrium() = default;
  Geometrium(const Geometrium&) = delete;
  Geometrium& operator=(const Geometrium&) = delete;

  //--- Elements ---------------------------------------------------------------

  // add_edge() throws std::invalid_argument when a == b, and add_face() when
  // the three edges do not join three distinct vertices into a triangle.
  VertexId add_vertex(Vec3 p);
  EdgeId   add_edge(VertexId a, VertexId b);
  FaceId   add_face(EdgeId a, EdgeId b, EdgeId c);
  bool     forms_triangle(EdgeId a, EdgeId b, EdgeId c) const;

  // A vertex at a + t (b - a), its attributes interpolated by the hooks.
  VertexId interpolate_vertex(VertexId a, VertexId b, double t);
  // A vertex at p, its attributes interpolated by the hooks as at a + t (b - a).
  // For splitting functions that compute the position their own way.
  VertexId add_interpolated_vertex(Vec3 p, VertexId a, VertexId b, double t);

  // Sets the attributes of piece from those of parent with the split hooks.
  // piece may be parent.
  void split_attributes(EdgeId piece, EdgeId parent);
  void split_attributes(FaceId piece, FaceId parent);

  //--- Splits -----------------------------------------------------------------

  // Splits e at m, a vertex the caller made, as GTS's refinement does: two
  // halves replace e, and each face on e is cut in two by an edge from m to
  // its opposite vertex. The face keeps one half and a new face takes the
  // other; the new face joins the surfaces of the face it comes from.
  // on_new_edge is called with each new edge as it is made, before it has
  // faces, which is when refine() computes its cost. Returns the halves, the
  // one from e's first vertex first.
  std::array<EdgeId, 2> split_edge(EdgeId e, VertexId m, const std::function<void(EdgeId)>& on_new_edge = {});

  // Splits f into four, as GTS's tessellation does: each edge is split at a
  // vertex made by split, unless memo has it split already, and three cuts
  // join the new vertices. f keeps the middle triangle; three new faces take
  // the corners and join the surfaces of f. An edge left without faces is
  // removed; an edge that other faces still use is kept, beside its halves.
  void split_face4(FaceId f, const Splitter& split, SplitMemo& memo);
  void split_face4(FaceId f, const Splitter& split) { SplitMemo m; split_face4(f, split, m); }

  //--- Collapses --------------------------------------------------------------

  // gts_edge_collapse_is_valid: whether collapsing e keeps the surface's
  // topology, by the constraints of Hoppe et al., "Mesh Optimization"
  // (1993). Edges with a duplicate are not collapsible.
  bool collapsible(EdgeId e) const;
  // gts_edge_collapse_creates_fold: whether collapsing e to a vertex at p
  // would make two neighbouring faces meet at an angle whose cosine squared
  // is above max_cos2, or make one face 1e4 times larger than its
  // neighbour.
  bool collapse_folds(EdgeId e, Vec3 p, double max_cos2) const;
  // Collapses e into a new vertex at p, as GTS's coarsening does, without
  // checking collapsible() or collapse_folds(). The new vertex's attributes
  // are merged from e's two vertices, and those of each kept edge from the
  // pair it was made from. on_removed_edge is called with each edge removed
  // besides e. Returns the new vertex, or an invalid index if it was
  // removed as well.
  VertexId collapse_edge(EdgeId e, Vec3 p, const std::function<void(EdgeId)>& on_removed_edge = {});

  // gts_vertex_replace: moves every edge of v to with, each becoming the
  // newest edge of with, and marks v dead. Throws std::logic_error when an
  // edge joins v and with, which would join with to itself.
  void replace_vertex(VertexId v, VertexId with);

  // gts_edge_replace: moves every face of e to with, which must join the same
  // two vertices. Each face keeps its place in its other rings and becomes the
  // newest face of with. e is left without faces.
  void replace_edge(EdgeId e, EdgeId with);

  // gts_vertices_merge: each vertex of vertices, in turn, takes the place
  // of the others within epsilon of it in x, y and z (and for which check,
  // when given, returns true), by replace_vertex(). Returns the vertices
  // left; the others are dead. An edge between two merged vertices is
  // removed with its faces; see QUIRKS.md. The edges that become duplicates
  // stay; merge_edges() merges them.
  std::vector<VertexId> merge_vertices(const std::vector<VertexId>& vertices, double epsilon,
                                       const std::function<bool(VertexId, VertexId)>& check = {});

  // gts_edges_merge: each edge of edges that has a duplicate gives it its
  // faces, by replace_edge(), and is removed. Returns the edges left.
  std::vector<EdgeId> merge_edges(const std::vector<EdgeId>& edges);

  // gts_triangle_revert: swaps the first two edges of f, which turns its
  // orientation over.
  void revert_face(FaceId f);

  // A vertex at p, its attributes merged by the hooks from those of a and b.
  VertexId add_merged_vertex(Vec3 p, VertexId a, VertexId b);

  // Gives f the edges a, b, c in place, as GTS code assigns t->e1, t->e2 and
  // t->e3. f leaves the rings of the edges it loses, keeps its place in the
  // rings of the edges it keeps, and becomes the newest face of the edges it
  // gains. Throws std::invalid_argument when a, b, c do not form a triangle.
  void set_face_edges(FaceId f, EdgeId a, EdgeId b, EdgeId c);

  // Unlinks from the rings and marks dead. kill_face() leaves the edges
  // linked; kill_edge() requires an edge without faces, and kill_vertex()
  // a vertex without edges.
  void kill_face(FaceId f);
  void kill_edge(EdgeId e);
  void kill_vertex(VertexId v);

  const Vertex& vertex(VertexId v) const { return m_vertices[v.i]; }
  const Edge&   edge(EdgeId e) const     { return m_edges[e.i]; }
  const Face&   face(FaceId f) const     { return m_faces[f.i]; }

  Vec3& position(VertexId v) { return m_vertices[v.i].p; }
  std::uint32_t& flags(VertexId v) { return m_vertices[v.i].flags; }
  std::uint32_t& flags(EdgeId e)   { return m_edges[e.i].flags; }
  std::uint32_t& flags(FaceId f)   { return m_faces[f.i].flags; }

  bool dead(VertexId v) const { return m_vertices[v.i].flags & kDead; }
  bool dead(EdgeId e) const   { return m_edges[e.i].flags & kDead; }
  bool dead(FaceId f) const   { return m_faces[f.i].flags & kDead; }

  // Sizes of the arrays, dead elements included.
  std::size_t vertex_capacity() const { return m_vertices.size(); }
  std::size_t edge_capacity() const   { return m_edges.size(); }
  std::size_t face_capacity() const   { return m_faces.size(); }

  //--- Adjacency --------------------------------------------------------------

  // Newest first, as GTS's lists.
  EdgeRing edges_of(VertexId v) const { return EdgeRing(this, v, m_vertices[v.i].edge); }
  FaceRing faces_of(EdgeId e) const   { return FaceRing(this, e, m_edges[e.i].face); }

  // The successor of e in the ring around v, and of f in the ring around e.
  EdgeId next_around(EdgeId e, VertexId v) const { return m_edges[e.i].next[vertex_slot(m_edges[e.i], v)]; }
  FaceId next_around(FaceId f, EdgeId e) const   { return m_faces[f.i].next[edge_slot(m_faces[f.i], e)]; }

  // gts_triangle_vertices: the vertices of f in its orientation.
  std::array<VertexId, 3> vertices_of(FaceId f) const;
  // gts_triangle_vertex: the vertex of f that is not on its first edge.
  VertexId opposite_vertex(FaceId f) const;
  // gts_triangle_vertex_opposite: the vertex of f that is not on e.
  VertexId opposite_vertex(FaceId f, EdgeId e) const;
  // gts_triangle_edge_opposite: the edge of f that does not touch v, or an
  // invalid index if v is not a vertex of f.
  EdgeId opposite_edge(FaceId f, VertexId v) const;
  // gts_triangles_common_edge: an edge of both faces, the first of f1's.
  EdgeId common_edge(FaceId f1, FaceId f2) const;
  // gts_triangle_vertices_edges: the edges of f in its order starting from e
  // (from its first edge if e is invalid), v[0] and v[1] the ends of e[0] in
  // f's orientation, v[2] the vertex opposite.
  struct VerticesEdges { std::array<VertexId, 3> v; std::array<EdgeId, 3> e; };
  VerticesEdges vertices_edges(FaceId f, EdgeId e = {}) const;
  // The other end of e.
  VertexId other_vertex(EdgeId e, VertexId v) const { return m_edges[e.i].v[0] == v ? m_edges[e.i].v[1] : m_edges[e.i].v[0]; }

  // gts_vertices_are_connected: the newest edge between a and b, if any.
  EdgeId find_edge(VertexId a, VertexId b) const;
  // gts_segment_is_duplicate: another edge with the same two vertices.
  EdgeId duplicate_of(EdgeId e) const;
  // gts_triangle_is_duplicate: another face with the same three edges.
  FaceId duplicate_of(FaceId f) const;
  // gts_triangles_are_compatible: whether f1 and f2, sharing e, are
  // oriented consistently.
  bool compatible(FaceId f1, FaceId f2, EdgeId e) const;

  //--- Geometry of elements, computed as GTS does -----------------------------

  // gts_triangle_orientation: orient2d of f's vertices in the xy plane,
  // positive when they run counterclockwise seen from +z.
  double orientation(FaceId f) const;

  Segment  segment(EdgeId e) const { return { m_vertices[m_edges[e.i].v[0].i].p, m_vertices[m_edges[e.i].v[1].i].p }; }
  // The triangle in the vertex order of vertices_of().
  Triangle triangle(FaceId f) const;
  // gts_triangle_normal, gts_triangle_area, gts_triangle_perimeter,
  // gts_triangle_quality, gts_triangles_angle.
  Vec3   normal(FaceId f) const;
  double area(FaceId f) const { return norm(normal(f)) / 2.; }
  double perimeter(FaceId f) const;
  double quality(FaceId f) const;
  double angle(FaceId f1, FaceId f2) const;

  //--- Surfaces ---------------------------------------------------------------

  Surface& new_surface(std::string name = {});
  // Removes s from its faces first, which may unlink them.
  void     remove_surface(Surface& s);
  std::size_t surface_number() const { return m_surfaces.size(); }

  // Copies the faces of src, a surface of src_geo, into a new surface here,
  // with new edges and vertices, and the attributes that exist here under
  // the same name and type. src_geo may be this Geometrium.
  Surface& import(const Geometrium& src_geo, const Surface& src);
  // The same into an existing surface dst. With src_geo this Geometrium, it
  // is gts_surface_copy (dst, src).
  void import(const Geometrium& src_geo, const Surface& src, Surface& dst);

  //--- Orphans and collection -------------------------------------------------

  const OrphanPolicy& orphan_policy() const { return m_orphans; }
  void set_orphan_policy(OrphanPolicy p)    { m_orphans = p; }

  // Removes dead elements and the unused ones selected by o, renumbers
  // everything that remains, attributes and surfaces included, keeping the
  // order of every ring, then calls the renumbering callouts.
  void collect(CollectOptions o = {});

  // Registers a function that collect() calls after renumbering. Callouts
  // run in the writing thread, in the order of registration, when the
  // Geometrium is already consistent under the new numbers.
  [[nodiscard]] Subscription on_renumber(RenumberCallout f) { return subscribe(m_renumber_callouts, std::move(f)); }

  // Registers a function that split_edge() and split_face4() call after each
  // split, when the elements are linked and their attributes set. The
  // callout may set attributes, for instance of the cuts, which have no
  // parent of their own kind; it must not change the topology.
  [[nodiscard]] Subscription on_split(SplitCallout f) { return subscribe(m_split_callouts, std::move(f)); }

  // Registers a function that collapse_edge() calls after each collapse,
  // when the elements are linked and their attributes set.
  [[nodiscard]] Subscription on_collapse(CollapseCallout f) { return subscribe(m_collapse_callouts, std::move(f)); }

  //--- Attributes -------------------------------------------------------------

  template <class T> Attribute<VertexId, T> add_vertex_attribute(std::string name, Hooks<T> h = {})
  { return Attribute<VertexId, T>(add_column<T>(m_vertex_columns, m_vertices.size(), std::move(name), std::move(h))); }
  template <class T> Attribute<EdgeId, T> add_edge_attribute(std::string name, Hooks<T> h = {})
  { return Attribute<EdgeId, T>(add_column<T>(m_edge_columns, m_edges.size(), std::move(name), std::move(h))); }
  template <class T> Attribute<FaceId, T> add_face_attribute(std::string name, Hooks<T> h = {})
  { return Attribute<FaceId, T>(add_column<T>(m_face_columns, m_faces.size(), std::move(name), std::move(h))); }

  // Truncates the arrays back to the given sizes. Only valid when nothing
  // older refers to the removed elements, as after a failed read.
  void truncate(std::size_t nv, std::size_t ne, std::size_t nf);

private:
  friend class Surface;

  template <class T>
  std::vector<T>* add_column(std::vector<detail::Column>& cols, std::size_t n, std::string name, Hooks<T> h)
  {
    cols.push_back(detail::Column::make<T>(std::move(name), n, std::move(h)));
    return cols.back().vector<T>();
  }

  static int vertex_slot(const Edge& e, VertexId v) { return e.v[0] == v ? 0 : 1; }
  static int edge_slot(const Face& f, EdgeId e)     { return f.e[0] == e ? 0 : f.e[1] == e ? 1 : 2; }

  template <class F>
  static Subscription subscribe(const std::shared_ptr<detail::CalloutRegistry<F>>& r, F f)
  {
    const std::uint64_t id = ++r->next_id;
    r->callouts.emplace_back(id, std::move(f));
    return Subscription([w = std::weak_ptr(r), id] {
      if (auto r = w.lock()) std::erase_if(r->callouts, [&](const auto& c) { return c.first == id; });
    });
  }
  // Calls the callouts on a copy of the list, so that one may subscribe or cancel.
  template <class F, class... A>
  static void call(const std::shared_ptr<detail::CalloutRegistry<F>>& r, const A&... a)
  {
    if (r->callouts.empty()) return;
    const auto callouts = r->callouts;
    for (const auto& c : callouts) c.second(a...);
  }

  // Adds piece to every surface that has parent.
  void join_surfaces(FaceId piece, FaceId parent);
  // split_face4's handling of one edge, running from a to b.
  struct EdgeHalves { VertexId v; EdgeId from_a, from_b; };
  EdgeHalves split_or_find(EdgeId e, VertexId a, VertexId b, const Splitter& split, SplitMemo& memo,
                           std::vector<EdgeSplitEvent>& events);

  // Removes f from every surface and unlinks it, without the release cascade.
  void destroy_face(FaceId f);
  // Unlinks e, which must have no faces, and marks dead its vertices left
  // without edges, as GTS's segment_destroy.
  void destroy_edge(EdgeId e);

  void link_edge(EdgeId e);
  void link_edge_slot(EdgeId e, int k);   // into the ring of its k-th vertex, as the newest
  void link_face(FaceId f);
  void link_face_slot(FaceId f, int k);   // into the ring of its k-th edge, as the newest
  void unlink_face_slot(FaceId f, int k); // out of the ring of its k-th edge
  void release_face(FaceId f); // called when a face leaves its last surface

  std::vector<Vertex> m_vertices;
  std::vector<Edge>   m_edges;
  std::vector<Face>   m_faces;
  std::vector<detail::Column> m_vertex_columns, m_edge_columns, m_face_columns;
  std::vector<std::unique_ptr<Surface>> m_surfaces;
  OrphanPolicy m_orphans;
  std::shared_ptr<detail::CalloutRegistry<RenumberCallout>> m_renumber_callouts = std::make_shared<detail::CalloutRegistry<RenumberCallout>>();
  std::shared_ptr<detail::CalloutRegistry<SplitCallout>>    m_split_callouts    = std::make_shared<detail::CalloutRegistry<SplitCallout>>();
  std::shared_ptr<detail::CalloutRegistry<CollapseCallout>> m_collapse_callouts = std::make_shared<detail::CalloutRegistry<CollapseCallout>>();
};

template <class ElemId, class CenterId>
typename Ring<ElemId, CenterId>::iterator Ring<ElemId, CenterId>::begin() const
{
  if ( ! m_last.valid()) return { m_g, m_center, m_last, m_last, true };
  return { m_g, m_center, m_g->next_around(m_last, m_center), m_last, false };
}

template <class ElemId, class CenterId>
typename Ring<ElemId, CenterId>::iterator& Ring<ElemId, CenterId>::iterator::operator++()
{
  if (cur == last) done = true;
  else             cur = g->next_around(cur, center);
  return *this;
}

// gts_edge_swap: replaces e, the edge of two faces of s, by the edge that
// joins their opposite vertices, using an existing one if there is. Two new
// faces replace the old ones in s, oriented as they were, and take their
// attributes with the split hooks; where a face with the same edges exists,
// it is added to s instead. Throws std::invalid_argument when e does not
// have exactly two faces in s or when they are not oriented consistently.
void swap_edge(Surface& s, EdgeId e);

// Sets an orphan policy for a scope, as GTS code set gts_allow_floating_*.
class KeepOrphans
{
public:
  KeepOrphans(Geometrium& g, OrphanPolicy p) : m_geo(g), m_saved(g.orphan_policy()) { g.set_orphan_policy(p); }
  ~KeepOrphans() { m_geo.set_orphan_policy(m_saved); }
  KeepOrphans(const KeepOrphans&) = delete;
  KeepOrphans& operator=(const KeepOrphans&) = delete;

private:
  Geometrium&  m_geo;
  OrphanPolicy m_saved;
};

} // namespace gts
