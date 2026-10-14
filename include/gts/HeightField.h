// Height fields: a triangulation with few vertices that approximates a field
// of heights, by greedy insertion, algorithm III of Garland and Heckbert,
// "Fast polygonal approximation of terrains and height fields" (1995), as
// GTS's example happrox.c.
//
// Each face of a Delaunay triangulation keeps the points its projection
// holds (FacePoints) and the one farthest from its plane, by cost. The
// point of largest cost is inserted as a vertex, the faces around it are
// updated, and so on until stop returns true.
//
// Copyright (C) 1999 Stéphane Popinet (GTS: examples/happrox.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-13.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "gts/Geometrium.h"

#include <functional>
#include <span>

namespace gts {

// The cost of point p in a face whose plane is z = a x + b y + c, given as
// plane = (a, b, c).
using HeightCost = std::function<double(Vec3 plane, Vec3 p)>;
// The default: |a x + b y + c - z|.
double height_cost(Vec3 plane, Vec3 p);
// That divided by |z|, or by z0 where |z| < z0.
HeightCost relative_height_cost(double z0);

// Stops the insertion when it returns true, given the largest cost left and
// the number of vertices, counted from 4 as GTS did.
using HeightStop = std::function<bool(double cost, std::size_t n_vertices)>;
// When there are n vertices.
HeightStop stop_at_vertices(std::size_t n);
// When the largest cost is below c.
HeightStop stop_below_cost(double c);

// happrox on a grid: z[j * width + i] is the height at (i, j), for a width
// and height of at least 2. The four corners are inserted first. Throws
// std::invalid_argument for a grid of the wrong size.
void approximate_grid(Surface& s, std::size_t width, std::size_t height, std::span<const double> z,
                      const HeightStop& stop, const HeightCost& cost = height_cost);

// happrox_list on scattered points, whose order is that of GTS's list. The
// triangle that encloses them is removed, at z = -1e30 until then, or kept
// at z = -10, or kept and closed by a face on its boundary, which makes s
// closed. A point with the x and y of a vertex of its face is not
// inserted, and its face is left out from then on, as in GTS. Throws
// std::invalid_argument for no points.
enum class Enclosing { Remove, Keep, Close };
void approximate_points(Surface& s, std::span<const Vec3> points, const HeightStop& stop,
                        const HeightCost& cost = height_cost, Enclosing enclosing = Enclosing::Remove);

} // namespace gts
