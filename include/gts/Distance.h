// Distances: from a point to an edge or a face, and the closest point on
// it; where a segment crosses a face; whether a point is inside a surface;
// and the distance between two surfaces, sampled on the faces and boundary
// edges of one and measured to the closest face or edge of the other.
//
// Copyright (C) 1999 Stéphane Popinet (GTS: point.c, bbtree.c, surface.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "gts/BBTree.h"
#include "gts/Stats.h"

#include <optional>

namespace gts {

// gts_point_segment_distance2, gts_point_segment_distance and
// gts_point_segment_closest, along e from its first vertex.
double distance2(const Geometrium& g, Vec3 p, EdgeId e);
double distance(const Geometrium& g, Vec3 p, EdgeId e);
Vec3   closest_point(const Geometrium& g, Vec3 p, EdgeId e);

// gts_point_triangle_distance2, gts_point_triangle_distance and
// gts_point_triangle_closest.
double distance2(const Geometrium& g, Vec3 p, FaceId f);
double distance(const Geometrium& g, Vec3 p, FaceId f);
Vec3   closest_point(const Geometrium& g, Vec3 p, FaceId f);

// gts_segment_triangle_intersection: where the segment a-b crosses f, if it
// does and does not lie in f's plane. With boundary, a crossing on f's
// edges counts.
std::optional<Vec3> intersection(const Geometrium& g, Vec3 a, Vec3 b, FaceId f, bool boundary);

// gts_point_is_inside_surface: whether p is inside the closed surface whose
// faces make tree, by the parity of the faces crossed by a ray from p
// towards +x. With is_open, the answer is reversed. Ties are broken by
// Simulation of Simplicity, which ranks p and the ray's end after all
// vertices.
bool is_inside(const BBTree& tree, Vec3 p, bool is_open = false);

// gts_bb_tree_triangle_distance and gts_bb_tree_segment_distance: the
// statistics of the distances, by d, from points sampled every delta on f
// or e of g to the closest leaf of tree.
Range face_distance(const BBTree& tree, const Geometrium& g, FaceId f, double delta,
                    const BBTree::DistanceFunc& d = {});
Range edge_distance(const BBTree& tree, const Geometrium& g, EdgeId e, double delta,
                    const BBTree::DistanceFunc& d = {});

// gts_bb_tree_surface_distance and gts_bb_tree_surface_boundary_distance:
// face_distance() for each face of s, edge_distance() for each boundary
// edge, with delta a fraction of the diagonal of tree's box, between 0 and
// 1. The mean and standard deviation are weighted by area, or by length;
// n counts the sampled points. A boundary edge is one with a single face
// in any surface, as in GTS.
Range surface_distance(const BBTree& tree, const Surface& s, double delta,
                       const BBTree::DistanceFunc& d = {});
Range boundary_distance(const BBTree& tree, const Surface& s, double delta,
                        const BBTree::DistanceFunc& d = {});

// gts_surface_distance: the distances from the faces of s1 to those of s2,
// and from the boundary edges of s1 to those of s2, with delta a fraction
// of the diagonal of s2's box, between 0 and 1. A range with nothing to
// measure is all zero.
struct SurfaceDistance { Range faces, boundary; };
SurfaceDistance distance(const Surface& s1, const Surface& s2, double delta);

} // namespace gts
