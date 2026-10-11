// Matrix4: 4 x 4 matrices for transforming points, as GtsMatrix: identity,
// scale, translation, rotation, the projection onto a face's plane, the
// product, transpose, determinant and inverse.
//
// Copyright (C) 1999 Stéphane Popinet (GTS: matrix.c, point.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-11.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "gts/Geometrium.h"

#include <array>
#include <optional>
#include <string>

namespace gts {

struct Matrix4
{
  std::array<std::array<double, 4>, 4> m{}; // m[row][column]

  std::array<double, 4>&       operator[](int i)       { return m[i]; }
  const std::array<double, 4>& operator[](int i) const { return m[i]; }
  friend bool operator==(const Matrix4&, const Matrix4&) = default;

  // gts_matrix_identity, gts_matrix_scale and gts_matrix_translate.
  static Matrix4 identity();
  static Matrix4 scale(Vec3 s);
  static Matrix4 translate(Vec3 t);
  // gts_matrix_rotate: by angle, in radians, round axis, which need not be
  // of unit length.
  static Matrix4 rotate(Vec3 axis, double angle);
  // gts_matrix_projection: the rotation whose columns are the unit vectors
  // along f's first edge, along the normal of f crossed with it, and along
  // the normal, from f's first vertex. Throws std::invalid_argument for a
  // degenerate face.
  static Matrix4 projection(const Geometrium& g, FaceId f);
};

// gts_matrix_product
Matrix4 operator*(const Matrix4& a, const Matrix4& b);
// gts_matrix_transpose
Matrix4 transpose(const Matrix4& m);
// gts_matrix_determinant, by cofactors along the first row.
double determinant(const Matrix4& m);
// gts_matrix_inverse: the adjoint divided by the determinant; none when the
// determinant is zero.
std::optional<Matrix4> inverse(const Matrix4& m);
// gts_matrix3_inverse: the inverse of the upper left 3 x 3 block, in a
// matrix whose fourth row and column are zero; none when singular.
std::optional<Matrix4> inverse3(const Matrix4& m);

// gts_point_transform: m p, with p as (x, y, z, 1).
Vec3 transform(const Matrix4& m, Vec3 p);
// Moves every vertex of s by m, as GTS's example transform.c.
void transform(Surface& s, const Matrix4& m);

// gts_matrix_print, gts_vector4_print and gts_vector_print: "%15.7g" for
// each element.
std::string to_string(const Matrix4& m);

} // namespace gts
