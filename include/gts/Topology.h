// Neighbourhoods: the faces, edges and vertices around an element, within a
// surface or in any; the fan of a vertex; contact edges; tetrahedra; the
// breadth-first traversal of a surface and its connected components.
//
// Lists are in the order of GTS's lists, which put the last element found
// first. Where a function takes a surface pointer, null means any face.
//
// Copyright (C) 1999 Stéphane Popinet (GTS: vertex.c, edge.c, face.c, triangle.c, segment.c, surface.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-11.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "gts/Geometrium.h"

#include <array>
#include <deque>
#include <functional>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace gts {

// gts_vertex_triangles and gts_vertex_faces: the faces around v, of s when
// given.
std::vector<FaceId> vertex_faces(const Geometrium& g, VertexId v, const Surface* s = nullptr);
inline std::vector<FaceId> vertex_faces(const Surface& s, VertexId v) { return vertex_faces(s.geometrium(), v, &s); }

// gts_vertex_neighbors: the vertices joined to v by an edge, by an edge
// with a face of s when s is given.
std::vector<VertexId> vertex_neighbors(const Geometrium& g, VertexId v, const Surface* s = nullptr);

// gts_vertex_fan_oriented: the edges opposite v in its faces of s, in the
// orientation of the faces. Empty where v has no face in s or an edge of v
// has more than two faces of s. Where the fan meets a non-manifold edge or
// a face oriented the other way, the list ends there.
std::vector<EdgeId> fan_oriented(const Surface& s, VertexId v);

// gts_face_neighbors and gts_triangle_neighbors: the other faces on the
// edges of f, of s when given, once for each edge they share with f. The
// number of neighbours (gts_face_neighbor_number) is the size of the list.
std::vector<FaceId> face_neighbors(const Geometrium& g, FaceId f, const Surface* s = nullptr);

// gts_face_is_compatible: whether f is oriented consistently with each of
// its neighbours in s.
bool is_compatible(const Surface& s, FaceId f);

// gts_edge_manifold_faces: the two faces of s on e, when e has exactly two.
std::optional<std::array<FaceId, 2>> manifold_faces(const Surface& s, EdgeId e);

// gts_edge_belongs_to_tetrahedron: whether two faces on e and two more
// faces close a tetrahedron around it.
bool belongs_to_tetrahedron(const Geometrium& g, EdgeId e);

// gts_edge_is_contact: the number of sets of faces around the vertices of
// e that meet only at e.
unsigned edge_contacts(const Geometrium& g, EdgeId e);

// gts_segments_from_vertices and gts_edges_from_vertices: the edges with a
// vertex in vertices, with a face of s when s is given.
std::vector<EdgeId> edges_from_vertices(const Geometrium& g, std::span<const VertexId> vertices, const Surface* s = nullptr);
// gts_triangles_from_edges and gts_faces_from_edges: the faces with an edge
// in edges, of s when given.
std::vector<FaceId> faces_from_edges(const Geometrium& g, std::span<const EdgeId> edges, const Surface* s = nullptr);
// gts_vertices_from_segments: the vertices of edges.
std::vector<VertexId> vertices_from_edges(const Geometrium& g, std::span<const EdgeId> edges);

// gts_surface_foreach_face_remove: removes from s the faces for which pred
// returns true, and returns their number.
std::size_t remove_faces_if(Surface& s, const std::function<bool(FaceId)>& pred);

// gts_surface_traverse_*: the faces of s in breadth-first order from start,
// each with its level, 0 for start, 1 for its neighbours and so on. s must
// not change during the traversal.
class SurfaceTraverse
{
public:
  SurfaceTraverse(const Surface& s, FaceId start);
  // The next face and its level; nullopt when no faces are left.
  std::optional<std::pair<FaceId, unsigned>> next();

private:
  const Surface&        m_s;
  std::deque<FaceId>    m_queue;
  std::vector<unsigned> m_level; // by face index, 0 for not reached, else level + 1
};

// gts_surface_split: the connected and manifold components of s, as new
// surfaces of its Geometrium, in the order of GTS's list. A component is
// grown from a face across edges with exactly two faces; the faces of a
// component on an edge with more than two of its faces then leave it and
// make components of their own. Faces of other surfaces take part as in
// GTS; see QUIRKS.md.
std::vector<Surface*> split_components(Surface& s);

} // namespace gts
