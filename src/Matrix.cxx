// Copyright (C) 1999 Stéphane Popinet (GTS: matrix.c, point.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-11.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/Matrix.h"

#include <cmath>
#include <format>
#include <stdexcept>

namespace gts {

Matrix4 Matrix4::identity()
{
  Matrix4 r;
  r[0][0] = r[1][1] = r[2][2] = r[3][3] = 1.;
  return r;
}

Matrix4 Matrix4::scale(Vec3 s)
{
  Matrix4 r;
  r[0][0] = s.x;
  r[1][1] = s.y;
  r[2][2] = s.z;
  r[3][3] = 1.;
  return r;
}

Matrix4 Matrix4::translate(Vec3 t)
{
  Matrix4 r = identity();
  r[0][3] = t.x;
  r[1][3] = t.y;
  r[2][3] = t.z;
  return r;
}

Matrix4 Matrix4::rotate(Vec3 a, double angle)
{
  // gts_vector_normalize
  const double n = norm(a);
  if (n > 0.) a = { a.x/n, a.y/n, a.z/n };
  const double c = std::cos(angle), c1 = 1. - c, s = std::sin(angle);
  Matrix4 m;
  m[0][0] = a.x*a.x*c1 + c;
  m[0][1] = a.x*a.y*c1 - a.z*s;
  m[0][2] = a.x*a.z*c1 + a.y*s;
  m[1][0] = a.y*a.x*c1 + a.z*s;
  m[1][1] = a.y*a.y*c1 + c;
  m[1][2] = a.y*a.z*c1 - a.x*s;
  m[2][0] = a.z*a.x*c1 - a.y*s;
  m[2][1] = a.z*a.y*c1 + a.x*s;
  m[2][2] = a.z*a.z*c1 + c;
  m[3][3] = 1.;
  return m;
}

Matrix4 Matrix4::projection(const Geometrium& g, FaceId f)
{
  const auto ve = g.vertices_edges(f);
  const Vec3 p1 = g.vertex(ve.v[0]).p, p2 = g.vertex(ve.v[1]).p, p3 = g.vertex(ve.v[2]).p;
  double x1 = p2.x - p1.x, y1 = p2.y - p1.y, z1 = p2.z - p1.z;
  double x2 = p3.x - p1.x, y2 = p3.y - p1.y, z2 = p3.z - p1.z;
  const double x3 = y1*z2 - z1*y2, y3 = z1*x2 - x1*z2, z3 = x1*y2 - y1*x2;
  x2 = y3*z1 - z3*y1; y2 = z3*x1 - x3*z1; z2 = x3*y1 - y3*x1;

  Matrix4 m;
  auto column = [&](int j, double x, double y, double z) {
    const double l = std::sqrt(x*x + y*y + z*z);
    if ( ! (l > 0.0)) throw std::invalid_argument("gts::Matrix4::projection: a degenerate face");
    m[0][j] = x/l; m[1][j] = y/l; m[2][j] = z/l;
  };
  column(0, x1, y1, z1);
  column(1, x2, y2, z2);
  column(2, x3, y3, z3);
  m[3][3] = 1.;
  return m;
}

Matrix4 operator*(const Matrix4& a, const Matrix4& b)
{
  Matrix4 m;
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      m[i][j] = a[i][0]*b[0][j] + a[i][1]*b[1][j] + a[i][2]*b[2][j] + a[i][3]*b[3][j];
  return m;
}

Matrix4 transpose(const Matrix4& m)
{
  Matrix4 t;
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      t[j][i] = m[i][j];
  return t;
}

namespace {

double det2x2(double a, double b, double c, double d) { return a*d - b*c; }

double det3x3(double a1, double a2, double a3, double b1, double b2, double b3,
              double c1, double c2, double c3)
{
  return a1*det2x2(b2, b3, c2, c3) - b1*det2x2(a2, a3, c2, c3) + c1*det2x2(a2, a3, b2, b3);
}

} // namespace

double determinant(const Matrix4& m)
{
  const double a1 = m[0][0], b1 = m[0][1], c1 = m[0][2], d1 = m[0][3];
  const double a2 = m[1][0], b2 = m[1][1], c2 = m[1][2], d2 = m[1][3];
  const double a3 = m[2][0], b3 = m[2][1], c3 = m[2][2], d3 = m[2][3];
  const double a4 = m[3][0], b4 = m[3][1], c4 = m[3][2], d4 = m[3][3];
  return a1*det3x3(b2, b3, b4, c2, c3, c4, d2, d3, d4)
       - b1*det3x3(a2, a3, a4, c2, c3, c4, d2, d3, d4)
       + c1*det3x3(a2, a3, a4, b2, b3, b4, d2, d3, d4)
       - d1*det3x3(a2, a3, a4, b2, b3, b4, c2, c3, c4);
}

std::optional<Matrix4> inverse(const Matrix4& m)
{
  const double det = determinant(m);
  if (det == 0.) return std::nullopt;
  const double a1 = m[0][0], b1 = m[0][1], c1 = m[0][2], d1 = m[0][3];
  const double a2 = m[1][0], b2 = m[1][1], c2 = m[1][2], d2 = m[1][3];
  const double a3 = m[2][0], b3 = m[2][1], c3 = m[2][2], d3 = m[2][3];
  const double a4 = m[3][0], b4 = m[3][1], c4 = m[3][2], d4 = m[3][3];
  // The adjoint; row and column labels reversed, as it is transposed.
  Matrix4 a;
  a[0][0] =   det3x3(b2, b3, b4, c2, c3, c4, d2, d3, d4);
  a[1][0] = - det3x3(a2, a3, a4, c2, c3, c4, d2, d3, d4);
  a[2][0] =   det3x3(a2, a3, a4, b2, b3, b4, d2, d3, d4);
  a[3][0] = - det3x3(a2, a3, a4, b2, b3, b4, c2, c3, c4);
  a[0][1] = - det3x3(b1, b3, b4, c1, c3, c4, d1, d3, d4);
  a[1][1] =   det3x3(a1, a3, a4, c1, c3, c4, d1, d3, d4);
  a[2][1] = - det3x3(a1, a3, a4, b1, b3, b4, d1, d3, d4);
  a[3][1] =   det3x3(a1, a3, a4, b1, b3, b4, c1, c3, c4);
  a[0][2] =   det3x3(b1, b2, b4, c1, c2, c4, d1, d2, d4);
  a[1][2] = - det3x3(a1, a2, a4, c1, c2, c4, d1, d2, d4);
  a[2][2] =   det3x3(a1, a2, a4, b1, b2, b4, d1, d2, d4);
  a[3][2] = - det3x3(a1, a2, a4, b1, b2, b4, c1, c2, c4);
  a[0][3] = - det3x3(b1, b2, b3, c1, c2, c3, d1, d2, d3);
  a[1][3] =   det3x3(a1, a2, a3, c1, c2, c3, d1, d2, d3);
  a[2][3] = - det3x3(a1, a2, a3, b1, b2, b3, d1, d2, d3);
  a[3][3] =   det3x3(a1, a2, a3, b1, b2, b3, c1, c2, c3);
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      a[i][j] /= det;
  return a;
}

std::optional<Matrix4> inverse3(const Matrix4& m)
{
  const double det = m[0][0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2]) -
                     m[0][1]*(m[1][0]*m[2][2] - m[2][0]*m[1][2]) +
                     m[0][2]*(m[1][0]*m[2][1] - m[2][0]*m[1][1]);
  if (det == 0.0) return std::nullopt;
  Matrix4 mi;
  mi[0][0] = (m[1][1]*m[2][2] - m[1][2]*m[2][1])/det;
  mi[0][1] = (m[2][1]*m[0][2] - m[0][1]*m[2][2])/det;
  mi[0][2] = (m[0][1]*m[1][2] - m[1][1]*m[0][2])/det;
  mi[1][0] = (m[1][2]*m[2][0] - m[1][0]*m[2][2])/det;
  mi[1][1] = (m[0][0]*m[2][2] - m[2][0]*m[0][2])/det;
  mi[1][2] = (m[1][0]*m[0][2] - m[0][0]*m[1][2])/det;
  mi[2][0] = (m[1][0]*m[2][1] - m[2][0]*m[1][1])/det;
  mi[2][1] = (m[2][0]*m[0][1] - m[0][0]*m[2][1])/det;
  mi[2][2] = (m[0][0]*m[1][1] - m[0][1]*m[1][0])/det;
  return mi;
}

Vec3 transform(const Matrix4& m, Vec3 p)
{
  return { m[0][0]*p.x + m[0][1]*p.y + m[0][2]*p.z + m[0][3],
           m[1][0]*p.x + m[1][1]*p.y + m[1][2]*p.z + m[1][3],
           m[2][0]*p.x + m[2][1]*p.y + m[2][2]*p.z + m[2][3] };
}

void transform(Surface& s, const Matrix4& m)
{
  Geometrium& g = s.geometrium();
  for (VertexId v : s.vertices()) g.position(v) = transform(m, g.vertex(v).p);
}

std::string to_string(const Matrix4& m)
{
  std::string r;
  for (int i = 0; i < 4; i++)
    r += std::format("{}[{:15.7g} {:15.7g} {:15.7g} {:15.7g}]{}", i ? " " : "[",
                     m[i][0], m[i][1], m[i][2], m[i][3], i == 3 ? "]\n" : "\n");
  return r;
}

} // namespace gts
