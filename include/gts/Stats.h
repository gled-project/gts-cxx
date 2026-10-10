// Measures and statistics of a surface: area, volume, whether it is closed,
// manifold and orientable, and the statistics GTS prints for a surface.
//
// Copyright (C) 1999 Stéphane Popinet (GTS: surface.c, misc.c, face.c, edge.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "gts/Geometrium.h"

#include <cfloat>
#include <ostream>
#include <string>
#include <vector>

namespace gts {

// GtsRange: count, extremes, mean and standard deviation of added values.
struct Range
{
  double   min = DBL_MAX, max = -DBL_MAX;
  double   sum = 0, sum2 = 0;
  double   mean = 0, stddev = 0;
  unsigned n = 0;

  void add(double v)
  {
    if (v < min) min = v;
    if (v > max) max = v;
    sum += v;
    sum2 += v * v;
    ++n;
  }

  // Computes mean and stddev; must be called after the last add().
  void update();
};

// gts_range_print: "min: %g mean: %g | %g max: %g".
std::string to_string(const Range& r);

double area(const Surface& s);
// Positive for a closed surface whose faces are oriented outwards.
double volume(const Surface& s);

// gts_surface_is_manifold: no edge has more than two faces of s.
bool is_manifold(const Surface& s);
// gts_surface_is_closed: every edge has exactly two faces of s.
bool is_closed(const Surface& s);
// gts_surface_is_orientable: every edge has at most two faces of s, and they
// are oriented consistently.
bool is_orientable(const Surface& s);

// Number of faces of s that have e as an edge.
unsigned face_number(const Surface& s, EdgeId e);

// gts_edge_is_boundary: the only face on e in s, or in any surface when s
// is null; invalid when e has no face or several.
FaceId boundary_face(const Geometrium& g, EdgeId e, const Surface* s = nullptr);
// gts_vertex_is_boundary: whether an edge of v is a boundary edge.
bool is_boundary(const Geometrium& g, VertexId v, const Surface* s = nullptr);
// gts_surface_boundary: the boundary edges of s, in the order of GTS's list.
std::vector<EdgeId> boundary(const Surface& s);

// The triangle cleanup of GTS's example cleanup.c: removes from s each face
// that has a duplicate in s, a face with the same three edges, in the order
// of GTS's list, so that the last of each set stays. Returns the number
// removed.
std::size_t remove_duplicate_faces(Surface& s);

// gts_surface_center_of_mass: the centre of mass of the volume bounded by
// s, and the signed volume.
struct CenterOfMass { Vec3 center; double volume; };
CenterOfMass center_of_mass(const Surface& s);
// gts_surface_center_of_area: the centre of the area of s, and the area.
struct CenterOfArea { Vec3 center; double area; };
CenterOfArea center_of_area(const Surface& s);

struct SurfaceStats
{
  unsigned n_faces = 0;
  unsigned n_incompatible_faces = 0;
  unsigned n_duplicate_faces = 0;
  unsigned n_duplicate_edges = 0;
  unsigned n_boundary_edges = 0;
  unsigned n_non_manifold_edges = 0;
  Range    edges_per_vertex, faces_per_edge;
};

struct QualityStats
{
  Range face_quality, face_area, edge_length, edge_angle;
};

SurfaceStats stats(const Surface& s);
QualityStats quality_stats(const Surface& s);

// gts_surface_print_stats: the connectivity and geometric statistics of s,
// in GTS's text.
void print_stats(const Surface& s, std::ostream& out);

} // namespace gts
