// Coarsening: collapsing edges in the order of a cost, and the
// volume-optimized vertex and cost of Lindstrom and Turk. The loop is over
// Geometrium::collapsible(), collapse_folds() and collapse_edge(), which
// merge attributes with the hooks and report to on_collapse() subscribers.
//
// Copyright (C) 1999 Stéphane Popinet (GTS: surface.c, vopt.c, matrix.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "gts/Geometrium.h"
#include "gts/Refine.h"

#include <cstddef>
#include <functional>

namespace gts {

// Gives the position of the vertex that replaces a collapsed edge.
using CollapsePosition = std::function<Vec3(const Geometrium&, EdgeId)>;
// Called with the cost of the next edge and the number of edges there would
// be after its collapse; returns true to stop.
using CoarsenStop = std::function<bool(double cost, std::size_t n_edges)>;

// gts_coarsen_stop_number and gts_coarsen_stop_cost.
CoarsenStop stop_below_edges(std::size_t min_edges);
CoarsenStop stop_above_cost(double max_cost);

// GtsVolumeOptimizedParams: the weights of the volume, boundary and shape
// terms of the cost.
struct VolumeOptimizedParams
{
  double volume_weight   = 0.5;
  double boundary_weight = 0.5;
  double shape_weight    = 0.0;
};

// gts_volume_optimized_vertex: the position that keeps the volume enclosed
// by the faces around e, and then the boundary and the shape of the faces,
// as well as it can. From Lindstrom and Turk, "Fast and memory efficient
// polygonal simplification" (1998) and "Evaluation of memoryless
// simplification" (1999).
Vec3 volume_optimized_vertex(const Geometrium& g, EdgeId e, const VolumeOptimizedParams& p = {});
// gts_volume_optimized_cost: the cost of collapsing e to that position.
double volume_optimized_cost(const Geometrium& g, EdgeId e, const VolumeOptimizedParams& p = {});

// The middle of e, as gts_segment_midvertex.
Vec3 midpoint(const Geometrium& g, EdgeId e);

// gts_surface_coarsen: collapses the edges of s in the order of cost,
// smallest first, until stop returns true. An edge that is not collapsible,
// or whose collapse would fold the surface by more than min_angle (in
// radians), goes to the end of the order. After a collapse, the costs of the
// edges within two edges of the new vertex are recomputed. Edges of equal
// cost are taken in the order of edge_order_key(). The default cost is the
// squared length, the default position the midpoint.
void coarsen(Surface& s, const CoarsenStop& stop, const EdgeCost& cost = {},
             const CollapsePosition& position = {}, double min_angle = 0.);

} // namespace gts
