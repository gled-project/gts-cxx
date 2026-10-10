// Refinement: splitting edges at new vertices. Tessellation splits every
// face into four, refinement splits the edges in the order of a cost, and a
// sphere is an icosahedron tessellated with vertices on the unit sphere. All
// three are loops over Geometrium::split_edge() and split_face4(), which set
// the attributes of what they make and report it to on_split() subscribers.
//
// Copyright (C) 1999 Stéphane Popinet (GTS: surface.c, segment.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "gts/Geometrium.h"

#include <cstddef>
#include <functional>

namespace gts {

using EdgeCost  = std::function<double(const Geometrium&, EdgeId)>;
// Called with the cost of the next edge and an estimate of the number of
// edges after it is split; returns true to stop.
using RefineStop = std::function<bool(double cost, std::size_t n_edges)>;

// gts_segment_midvertex: the vertex at the middle of e.
VertexId midvertex(Geometrium& g, EdgeId e);
// The middle of e, projected onto the unit sphere.
VertexId unit_sphere_arc_midvertex(Geometrium& g, EdgeId e);

// gts_surface_tessellate: splits every face of s into four, splitting each
// edge once with split (by default unit_sphere_arc_midvertex). Faces of other
// surfaces that share an edge with s keep the edge.
void tessellate(Surface& s, const Splitter& split = {});

// gts_surface_refine: splits the edges of s in the order of cost, smallest
// first, until stop returns true. Splitting an edge splits each face on it
// into two; the new edges join the order. The default cost is minus the
// squared length, which splits the longest edges first; the default split
// is midvertex.
void refine(Surface& s, const RefineStop& stop, const EdgeCost& cost = {}, const Splitter& split = {});

// gts_surface_generate_sphere: adds to s a unit sphere made from an
// icosahedron (order 1), tessellated order - 1 times. Throws
// std::invalid_argument for order 0.
Surface& generate_sphere(Surface& s, unsigned order);

} // namespace gts
