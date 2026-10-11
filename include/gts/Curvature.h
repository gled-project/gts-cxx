// Discrete curvature at a vertex of a surface, by the operators of Meyer,
// Desbrun, Schröder and Barr, "Discrete differential-geometry operators for
// triangulated 2-manifolds" (2002), as GTS has them. None is defined at a
// boundary vertex.
//
// Copyright (C) 1999 Stéphane Popinet, 2002 Ray Jones (GTS: curvature.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "gts/Geometrium.h"

#include <array>
#include <optional>
#include <utility>

namespace gts {

// gts_vertex_mean_curvature_normal: the mean curvature normal Kh at v, a
// vector of length twice the mean curvature. None at a boundary vertex, or
// where v's faces have no area or do not form an oriented fan.
std::optional<Vec3> mean_curvature_normal(const Surface& s, VertexId v);

// gts_vertex_gaussian_curvature: the Gaussian curvature Kg at v. None at a
// boundary vertex, or where v's faces have no area or do not form an
// oriented fan.
std::optional<double> gaussian_curvature(const Surface& s, VertexId v);

// gts_vertex_principal_curvatures: K1 and K2 from the mean curvature Kh, a
// scalar (half the length of the mean curvature normal), and Kg.
std::pair<double, double> principal_curvatures(double kh, double kg);

// gts_vertex_principal_directions: the principal directions e1 and e2 at v,
// from the mean curvature normal and the Gaussian curvature. Throws
// std::invalid_argument when an edge of v does not have exactly two faces
// in s, or when no normal can be found.
std::array<Vec3, 2> principal_directions(const Surface& s, VertexId v, Vec3 kh, double kg);

} // namespace gts
