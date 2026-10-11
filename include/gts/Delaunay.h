// Delaunay triangulation in the xy plane: point location, insertion and
// removal of vertices, constrained edges, the convex hull, and conforming and
// refined triangulations after Ruppert (1995) and Shewchuk (1996); and the
// planar tests and constructions they use.
//
// A surface is a triangulation of the projections of its vertices on the xy
// plane; z is carried along but not used, except by interpolate_height().
// A triangulation starts as one face enclosing all the points, usually from
// enclosing_triangle(), to which vertices are added, then constraints. A
// constraint is an edge with the flag kConstraint, which swaps and
// insertions keep.
//
// Copyright (C) 1999 Stéphane Popinet (GTS: cdt.c, refine.c, point.c, triangle.c, segment.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-11.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "gts/Geometrium.h"

#include <functional>
#include <optional>
#include <span>
#include <vector>

namespace gts {

//==============================================================================
// Planar tests and constructions
//==============================================================================

// GtsIntersect.
enum class Intersect { Out = -1, On = 0, In = 1 };

// gts_point_is_in_triangle: whether the projection of p lies inside, on the
// boundary of, or outside that of t, whose vertices run counterclockwise.
Intersect point_in_triangle(Vec3 p, const Triangle& t);
// gts_point_in_triangle_circle: positive when the projection of p lies
// inside the circle through the projections of t's vertices, which run
// counterclockwise.
double in_triangle_circle(Vec3 p, const Triangle& t);
// gts_segments_are_intersecting: In when the projections of a and b
// cross, On when an end of one lies on the other.
Intersect segments_intersect(const Segment& a, const Segment& b);
// gts_triangle_circumcircle_center: the centre of the circle through the
// projections of t's vertices, at z = 0; none for collinear vertices.
std::optional<Vec3> circumcircle_center(const Triangle& t);
// gts_triangle_interpolate_height: the z of the plane of f at the x and y
// of p; the mean z of f's vertices when f is degenerate in the xy plane.
double interpolate_height(const Geometrium& g, FaceId f, Vec3 p);
// gts_triangle_enclosing: a new equilateral face in the plane z = 0, on new
// vertices and edges of g, enclosing the rectangle round the projections of
// points scaled by scale about its centre. Invalid for no points.
FaceId enclosing_triangle(Geometrium& g, std::span<const Vec3> points, double scale);

//==============================================================================
// Triangulation
//==============================================================================

inline bool is_constraint(const Geometrium& g, EdgeId e) { return g.edge(e).flags & kConstraint; }
inline void set_constraint(Geometrium& g, EdgeId e, bool on = true)
{
  if (on) g.flags(e) |= kConstraint; else g.flags(e) &= ~kConstraint;
}

// gts_point_locate: the face of s whose projection contains that of p, by a
// walk from guess, or from the nearest of n^1/3 faces of s for n faces.
// The projection of s must be connected, without holes, and convex. Invalid
// when p lies outside.
FaceId locate(const Surface& s, Vec3 p, FaceId guess = {});

// gts_delaunay_add_vertex_to_face: adds v, which lies in f, to the
// triangulation s, and swaps edges until s is Delaunay again. Returns an
// invalid index when v was added or is a vertex of f already, or the vertex
// of f with v's x and y, when there is one; v is not added then.
VertexId delaunay_add_vertex_to_face(Surface& s, VertexId v, FaceId f);
// gts_delaunay_add_vertex: delaunay_add_vertex_to_face() for the face that
// locate() finds; returns v when v lies outside s.
VertexId delaunay_add_vertex(Surface& s, VertexId v, FaceId guess = {});
// gts_delaunay_remove_vertex: removes the faces of s around v, then fills
// the hole with a Delaunay triangulation. v must not be on a constraint.
void delaunay_remove_vertex(Surface& s, VertexId v);
// gts_delaunay_add_constraint: makes c, an edge between two vertices of s,
// a constraint and an edge of s, removing the edges it crosses and
// triangulating the two sides. Returns the constraints that crossed c,
// which are no longer edges of s. Where an edge of s joins the two
// vertices already, c stays out of s; see QUIRKS.md.
std::vector<EdgeId> delaunay_add_constraint(Surface& s, EdgeId c);
// gts_delaunay_check: a face of s whose circumcircle has a neighbouring
// vertex inside, or an invalid index when s is Delaunay. Constraints are
// not taken into account.
FaceId delaunay_check(const Surface& s);
// gts_delaunay_remove_hull: removes from s the faces on its boundary
// edges that are not constraints, until only constraints are left on the
// boundary.
void delaunay_remove_hull(Surface& s);

// Whether v encroaches on e, a constraint of s.
using Encroaches = std::function<bool(const Surface& s, VertexId v, EdgeId e)>;
// gts_vertex_encroaches_edge: whether v lies strictly inside the circle
// with e as diameter, in the xy plane.
bool vertex_encroaches_edge(const Geometrium& g, VertexId v, EdgeId e);
// gts_edge_is_encroached: a vertex opposite e in a face of s that
// encroaches on e, by encroaches (vertex_encroaches_edge() by default), or
// an invalid index.
VertexId edge_is_encroached(const Surface& s, EdgeId e, const Encroaches& encroaches = {});

// gts_delaunay_conform: splits the constraints of s that vertices encroach
// on, until none is, or until steiner_max vertices have been added when it
// is not negative. The pieces of a constraint are constraints. Returns the
// number of encroached constraints left.
unsigned delaunay_conform(Surface& s, int steiner_max, const Encroaches& encroaches = {});

// The cost of a face for refinement; 0 means the face is good.
using FaceCost = std::function<double(const Geometrium&, FaceId)>;
// gts_delaunay_refine: Ruppert's refinement of s, a conforming Delaunay
// triangulation: inserts the circumcentre of the face of smallest non-zero
// cost, or splits the constraints it would encroach on instead, until no
// face has a cost or steiner_max vertices were added when it is not
// negative. Returns the number of faces with a cost left. Throws
// std::runtime_error when a circumcentre falls outside s or on a vertex,
// where GTS failed an assertion.
unsigned delaunay_refine(Surface& s, int steiner_max, const FaceCost& cost, const Encroaches& encroaches = {});

} // namespace gts
