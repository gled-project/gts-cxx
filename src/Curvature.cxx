// Copyright (C) 1999 Stéphane Popinet, 2002 Ray Jones (GTS: curvature.c, vertex.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/Curvature.h"
#include "gts/Stats.h"
#include "gts/Topology.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace gts {

namespace {

Vec3 pos(const Geometrium& g, VertexId v) { return g.vertex(v).p; }

// gts_vector_normalize
Vec3 normalized(Vec3 v)
{
  const double n = norm(v);
  return n > 0. ? Vec3{ v.x / n, v.y / n, v.z / n } : v;
}

bool angle_obtuse(const Geometrium& g, VertexId v, FaceId f)
{
  const Edge& e = g.edge(g.opposite_edge(f, v));
  const Vec3 vec1 = pos(g, e.v[0]) - pos(g, v);
  const Vec3 vec2 = pos(g, e.v[1]) - pos(g, v);
  return dot(vec1, vec2) < 0.0;
}

bool triangle_obtuse(const Geometrium& g, VertexId v, FaceId f)
{
  const Edge& e = g.edge(g.opposite_edge(f, v));
  return angle_obtuse(g, v, f) || angle_obtuse(g, e.v[0], f) || angle_obtuse(g, e.v[1], f);
}

// Appendix B of Meyer et al.
double cotan(const Geometrium& g, VertexId vo, VertexId v1, VertexId v2)
{
  const Vec3 u = pos(g, v1) - pos(g, vo);
  const Vec3 v = pos(g, v2) - pos(g, vo);
  const double udotv = dot(u, v);
  const double denom = std::sqrt(dot(u, u)*dot(v, v) - udotv*udotv);
  // denom is zero when u == v; 0 suits the callers.
  if (denom == 0.0) return 0.0;
  return udotv/denom;
}

// Appendix B and the caption of Table 1 of Meyer et al.
double angle_from_cotan(const Geometrium& g, VertexId vo, VertexId v1, VertexId v2)
{
  const Vec3 u = pos(g, v1) - pos(g, vo);
  const Vec3 v = pos(g, v2) - pos(g, vo);
  const double udotv = dot(u, v);
  const double denom = std::sqrt(dot(u, u)*dot(v, v) - udotv*udotv);
  return std::fabs(std::atan2(denom, udotv));
}

// Section 3.3 of Meyer et al.
double region_area(const Geometrium& g, VertexId v, FaceId f)
{
  const double area = g.area(f);
  if (area == 0.0) return 0.0;
  if (triangle_obtuse(g, v, f))
    return angle_obtuse(g, v, f) ? area/2.0 : area/4.0;
  const Edge& e = g.edge(g.opposite_edge(f, v));
  return (cotan(g, e.v[0], v, e.v[1]) * distance2(pos(g, v), pos(g, e.v[1])) +
          cotan(g, e.v[1], v, e.v[0]) * distance2(pos(g, v), pos(g, e.v[0])))/8.0;
}

// Sum of region_area() over v's faces; the faces must be there.
double mixed_area(const Surface& s, VertexId v, const std::vector<FaceId>& faces)
{
  double area = 0.0;
  for (FaceId f : faces) area += region_area(s.geometrium(), v, f);
  return area;
}

} // namespace

std::optional<Vec3> mean_curvature_normal(const Surface& s, VertexId v)
{
  const Geometrium& g = s.geometrium();
  // The operator is not defined at a boundary.
  if (is_boundary(g, v, &s)) return std::nullopt;
  const std::vector<FaceId> faces = vertex_faces(s, v);
  if (faces.empty()) return std::nullopt;
  const std::vector<EdgeId> edges = fan_oriented(s, v);
  if (edges.empty()) return std::nullopt;

  const double area = mixed_area(s, v, faces);
  const Vec3 p = pos(g, v);
  double kh[3] = { 0.0, 0.0, 0.0 };
  for (EdgeId e : edges)
  {
    const VertexId v1 = g.edge(e).v[0], v2 = g.edge(e).v[1];
    const Vec3 p1 = pos(g, v1), p2 = pos(g, v2);
    double temp = cotan(g, v1, v, v2);
    kh[0] += temp*(p2.x - p.x);
    kh[1] += temp*(p2.y - p.y);
    kh[2] += temp*(p2.z - p.z);
    temp = cotan(g, v2, v, v1);
    kh[0] += temp*(p1.x - p.x);
    kh[1] += temp*(p1.y - p.y);
    kh[2] += temp*(p1.z - p.z);
  }
  if ( ! (area > 0.0)) return std::nullopt;
  return Vec3{ kh[0]/(2*area), kh[1]/(2*area), kh[2]/(2*area) };
}

std::optional<double> mean_curvature(const Surface& s, VertexId v)
{
  const std::optional<Vec3> kh = mean_curvature_normal(s, v);
  if ( ! kh) return std::nullopt;
  const Geometrium& g = s.geometrium();
  Vec3 n = {};
  for (FaceId f : vertex_faces(s, v)) n = n + g.normal(f);
  const double h = norm(*kh) / 2;
  return dot(*kh, n) > 0 ? - h : h;
}

std::optional<double> gaussian_curvature(const Surface& s, VertexId v)
{
  const Geometrium& g = s.geometrium();
  if (is_boundary(g, v, &s)) return std::nullopt;
  const std::vector<FaceId> faces = vertex_faces(s, v);
  if (faces.empty()) return std::nullopt;
  const std::vector<EdgeId> edges = fan_oriented(s, v);
  if (edges.empty()) return std::nullopt;

  const double area = mixed_area(s, v, faces);
  double angle_sum = 0.0;
  for (EdgeId e : edges)
    angle_sum += angle_from_cotan(g, v, g.edge(e).v[0], g.edge(e).v[1]);
  if ( ! (area > 0.0)) return std::nullopt;
  return (2.0*std::numbers::pi - angle_sum)/area;
}

std::pair<double, double> principal_curvatures(double kh, double kg)
{
  double temp = kh*kh - kg;
  if (temp < 0.0) temp = 0.0;
  temp = std::sqrt(temp);
  return { kh + temp, kh - temp };
}

std::array<Vec3, 2> principal_directions(const Surface& s, VertexId v, Vec3 kh, double kg)
{
  const Geometrium& g = s.geometrium();

  // the unit normal
  const double norm_kh = std::sqrt(dot(kh, kh));
  Vec3 n;
  if (norm_kh > 0.0)
    n = { kh.x / norm_kh, kh.y / norm_kh, kh.z / norm_kh };
  else
  {
    // A point of zero mean curvature, flat or a saddle: the normal is the
    // sum of those of the faces.
    for (FaceId f : vertex_faces(s, v))
      n = n + g.normal(f);
    if ( ! (norm(n) > 0.0))
      throw std::invalid_argument("gts::principal_directions: no normal at the vertex");
    n = normalized(n);
  }

  // A basis from n: basis1 along a component that is not n's largest.
  Vec3 basis1;
  if (std::fabs(n.x) > std::fabs(n.y)) basis1.y = 1.0;
  else                                  basis1.x = 1.0;
  const Vec3 basis2 = normalized(cross(n, basis1));
  basis1 = normalized(cross(n, basis2));

  double aterm_da = 0., bterm_da = 0., cterm_da = 0., const_da = 0.;
  double aterm_db = 0., bterm_db = 0., cterm_db = 0., const_db = 0.;
  struct Sample { double weight, kappa, d1, d2; };
  std::vector<Sample> samples;

  for (EdgeId e : g.edges_of(v))
  {
    // the two faces of e in s
    FaceId f1, f2;
    for (FaceId t : g.faces_of(e))
    {
      if ( ! s.contains(t)) continue;
      if ( ! f1.valid()) f1 = t;
      else if ( ! f2.valid()) f2 = t;
      else { f2 = FaceId(); break; }
    }
    if ( ! f2.valid())
      throw std::invalid_argument("gts::principal_directions: an edge of the vertex does not have two faces");

    // Section 5 of Meyer et al., solving for the curvature tensor
    // B = [a b; b c]; the comments of GTS's curvature.c explain the steps.
    const VertexId other = g.other_vertex(e, v);
    const Vec3 vec_edge = pos(g, other) - pos(g, v);
    const double ve2 = dot(vec_edge, vec_edge);
    const double vdotn = dot(vec_edge, n);
    const double kappa = 2.0 * vdotn / ve2;

    double weight = 0.0;
    for (FaceId f : { f1, f2 })
    {
      if ( ! triangle_obtuse(g, v, f))
        weight += ve2 * cotan(g, g.opposite_vertex(f, e), g.edge(e).v[0], g.edge(e).v[1]) / 8.0;
      else if (angle_obtuse(g, v, f))
        weight += ve2 * g.area(f) / 4.0;
      else
        weight += ve2 * g.area(f) / 8.0;
    }

    // The edge projected perpendicular to n, in the basis.
    const Vec3 d = normalized({ vec_edge.x - vdotn * n.x, vec_edge.y - vdotn * n.y, vec_edge.z - vdotn * n.z });
    const double d1 = dot(d, basis1);
    const double d2 = dot(d, basis2);
    samples.push_back({ weight, kappa, d1, d2 });

    aterm_da += weight * d1 * d1 * d1 * d1;
    bterm_da += weight * d1 * d1 * 2 * d1 * d2;
    cterm_da += weight * d1 * d1 * d2 * d2;
    const_da += weight * d1 * d1 * (- kappa);

    aterm_db += weight * d1 * d2 * d1 * d1;
    bterm_db += weight * d1 * d2 * 2 * d1 * d2;
    cterm_db += weight * d1 * d2 * d2 * d2;
    const_db += weight * d1 * d2 * (- kappa);
  }

  // a + c = |Kh|
  aterm_da -= cterm_da;
  const_da += cterm_da * norm_kh;
  aterm_db -= cterm_db;
  const_db += cterm_db * norm_kh;

  double eig[2];
  if ((aterm_da * bterm_db - aterm_db * bterm_da) != 0.0 && (const_da != 0.0 || const_db != 0.0))
  {
    // linsolve, from Maple
    const double m11 = aterm_da, m12 = bterm_da, b1 = -const_da;
    const double m21 = aterm_db, m22 = bterm_db, b2 = -const_db;
    const double temp = 1.0 / (m21*m12 - m11*m22);
    const double a = (m12*b2 - m22*b1)*temp;
    const double b = (m11*b2 - m21*b1)*temp;
    const double c = norm_kh - a;
    // eigenvector, from Maple: the largest eigenvector of [a b; b c]
    eig[0] = b == 0.0 ? 0.0 : -(c - a - std::sqrt(c*c - 2*a*c + a*a + 4*b*b))/(2*b);
    eig[1] = 1.0;
  }
  else
  {
    // the region of v is planar
    eig[0] = 1.0;
    eig[1] = 0.0;
  }

  // Which eigenvector goes with which curvature: the assignment that
  // matches the curvatures along the edges better.
  const auto [k1, k2] = principal_curvatures(0.5 * norm_kh, kg);
  double err_e1 = 0.0, err_e2 = 0.0;
  for (const Sample& smp : samples)
  {
    double temp1 = std::fabs(eig[0] * smp.d1 + eig[1] * smp.d2);
    temp1 = temp1 * temp1;
    double temp2 = std::fabs(eig[1] * smp.d1 - eig[0] * smp.d2);
    temp2 = temp2 * temp2;
    double delta = k1 * temp1 + k2 * temp2 - smp.kappa;
    err_e1 += smp.weight * delta * delta;
    delta = k2 * temp1 + k1 * temp2 - smp.kappa;
    err_e2 += smp.weight * delta * delta;
  }
  if (err_e2 < err_e1)
  {
    const double temp = eig[0];
    eig[0] = eig[1];
    eig[1] = -temp;
  }

  const Vec3 e1 = normalized({ eig[0] * basis1.x + eig[1] * basis2.x,
                               eig[0] * basis1.y + eig[1] * basis2.y,
                               eig[0] * basis1.z + eig[1] * basis2.z });
  const Vec3 e2 = normalized(cross(n, e1));
  return { e1, e2 };
}

} // namespace gts
