// Copyright (C) 1999 Stéphane Popinet (GTS: point.c, bbtree.c, surface.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/Distance.h"
#include "gts/Predicates.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <stdexcept>

namespace gts {

//==============================================================================
// Point and edge
//==============================================================================

double distance2(const Geometrium& g, Vec3 p, EdgeId e)
{
  const Segment s = g.segment(e);
  const Vec3 p1 = s.a, p2 = s.b;
  const double ns2 = gts::distance2(p1, p2);
  if (ns2 == 0.0)
    return gts::distance2(p, p1);
  const double t = ((p2.x - p1.x)*(p.x - p1.x) +
                    (p2.y - p1.y)*(p.y - p1.y) +
                    (p2.z - p1.z)*(p.z - p1.z))/ns2;
  if (t > 1.0)
    return gts::distance2(p, p2);
  if (t < 0.0)
    return gts::distance2(p, p1);
  const double x = (1. - t)*p1.x + t*p2.x - p.x;
  const double y = (1. - t)*p1.y + t*p2.y - p.y;
  const double z = (1. - t)*p1.z + t*p2.z - p.z;
  return x*x + y*y + z*z;
}

double distance(const Geometrium& g, Vec3 p, EdgeId e)
{
  return std::sqrt(distance2(g, p, e));
}

Vec3 closest_point(const Geometrium& g, Vec3 p, EdgeId e)
{
  const Segment s = g.segment(e);
  const Vec3 p1 = s.a, p2 = s.b;
  const double ns2 = gts::distance2(p1, p2);
  if (ns2 == 0.0)
    return p1;
  const double t = ((p2.x - p1.x)*(p.x - p1.x) +
                    (p2.y - p1.y)*(p.y - p1.y) +
                    (p2.z - p1.z)*(p.z - p1.z))/ns2;
  if (t > 1.0) return p2;
  if (t < 0.0) return p1;
  return { (1. - t)*p1.x + t*p2.x, (1. - t)*p1.y + t*p2.y, (1. - t)*p1.z + t*p2.z };
}

//==============================================================================
// Point and face
//==============================================================================

namespace {

// The parameters of the point of f's plane closest to p, as GTS computes
// them; det is zero for a degenerate face.
struct FacePlane
{
  Geometrium::VerticesEdges ve;
  Vec3   p1, p1p2, p1p3, pp1;
  double det, t1, t2;
};

FacePlane face_plane(const Geometrium& g, Vec3 p, FaceId f)
{
  FacePlane fp;
  fp.ve = g.vertices_edges(f);
  fp.p1 = g.vertex(fp.ve.v[0]).p;
  const Vec3 p2 = g.vertex(fp.ve.v[1]).p, p3 = g.vertex(fp.ve.v[2]).p;
  fp.p1p2 = p2 - fp.p1;
  fp.p1p3 = p3 - fp.p1;
  fp.pp1  = fp.p1 - p;

  const double B = dot(fp.p1p3, fp.p1p2);
  const double E = dot(fp.p1p2, fp.p1p2);
  const double C = dot(fp.p1p3, fp.p1p3);
  fp.det = B*B - E*C;
  if (fp.det == 0.) return fp;

  const double A = dot(fp.p1p3, fp.pp1);
  const double D = dot(fp.p1p2, fp.pp1);
  fp.t1 = (D*C - A*B)/fp.det;
  fp.t2 = (A*E - D*B)/fp.det;
  return fp;
}

} // namespace

double distance2(const Geometrium& g, Vec3 p, FaceId f)
{
  const FacePlane fp = face_plane(g, p, f);
  const auto& e = fp.ve.e;
  if (fp.det == 0.) // p1p2 and p1p3 are collinear
  {
    const double d1 = distance2(g, p, e[0]);
    const double d2 = distance2(g, p, e[2]);
    return d1 < d2 ? d1 : d2;
  }
  if (fp.t1 < 0.) return distance2(g, p, e[2]);
  if (fp.t2 < 0.) return distance2(g, p, e[0]);
  if (fp.t1 + fp.t2 > 1.) return distance2(g, p, e[1]);

  const double x = fp.pp1.x + fp.t1*fp.p1p2.x + fp.t2*fp.p1p3.x;
  const double y = fp.pp1.y + fp.t1*fp.p1p2.y + fp.t2*fp.p1p3.y;
  const double z = fp.pp1.z + fp.t1*fp.p1p2.z + fp.t2*fp.p1p3.z;
  return x*x + y*y + z*z;
}

double distance(const Geometrium& g, Vec3 p, FaceId f)
{
  return std::sqrt(distance2(g, p, f));
}

Vec3 closest_point(const Geometrium& g, Vec3 p, FaceId f)
{
  const FacePlane fp = face_plane(g, p, f);
  const auto& e = fp.ve.e;
  if (fp.det == 0.)
  {
    const Vec3 cp = closest_point(g, p, e[0]);
    const Vec3 closest = closest_point(g, p, e[2]);
    return gts::distance2(cp, p) < gts::distance2(closest, p) ? cp : closest;
  }
  if (fp.t1 < 0.) return closest_point(g, p, e[2]);
  if (fp.t2 < 0.) return closest_point(g, p, e[0]);
  if (fp.t1 + fp.t2 > 1.) return closest_point(g, p, e[1]);
  return { fp.p1.x + fp.t1*fp.p1p2.x + fp.t2*fp.p1p3.x,
           fp.p1.y + fp.t1*fp.p1p2.y + fp.t2*fp.p1p3.y,
           fp.p1.z + fp.t1*fp.p1p2.z + fp.t2*fp.p1p3.z };
}

//==============================================================================
// Segment and face
//==============================================================================

std::optional<Vec3> intersection(const Geometrium& g, Vec3 d, Vec3 e, FaceId f, bool boundary)
{
  const Edge& e1 = g.edge(g.face(f).e[0]);
  const Vec3 A = g.vertex(e1.v[0]).p;
  const Vec3 B = g.vertex(e1.v[1]).p;
  const Vec3 C = g.vertex(g.opposite_vertex(f)).p;
  Vec3 D = d, E = e;

  double ABCE = orient3d(A, B, C, E);
  double ABCD = orient3d(A, B, C, D);
  if (ABCE < 0.0 || ABCD > 0.0)
  {
    std::swap(D, E);
    std::swap(ABCE, ABCD);
  }
  if (ABCE < 0.0 || ABCD > 0.0)
    return std::nullopt;
  const double ADCE = orient3d(A, D, C, E);
  if ((boundary && ADCE < 0.) || ( ! boundary && ADCE <= 0.))
    return std::nullopt;
  const double ABDE = orient3d(A, B, D, E);
  if ((boundary && ABDE < 0.) || ( ! boundary && ABDE <= 0.))
    return std::nullopt;
  const double BCDE = orient3d(B, C, D, E);
  if ((boundary && BCDE < 0.) || ( ! boundary && BCDE <= 0.))
    return std::nullopt;
  if (ABCE == 0.0)
  {
    if (ABCD == 0.0) // the segment lies in f's plane
      return std::nullopt;
    return E;
  }
  if (ABCD == 0.0)
    return D;
  if (boundary) // corners of f
  {
    if (ABDE == 0.)
    {
      if (ADCE == 0.) return A;
      if (BCDE == 0.) return B;
    }
    else if (BCDE == 0. && ADCE == 0.)
      return C;
  }
  const double c = ABCE/(ABCE - ABCD);
  return Vec3{ E.x + c*(D.x - E.x), E.y + c*(D.y - E.y), E.z + c*(D.z - E.z) };
}

//==============================================================================
// Inside
//==============================================================================

namespace {

struct RankedPoint { Vec3 p; unsigned rank; };

int orient_sos(const RankedPoint& a, const RankedPoint& b, const RankedPoint& c, const RankedPoint& d)
{
  const Vec3* p[4] = { &a.p, &b.p, &c.p, &d.p };
  const unsigned r[4] = { a.rank, b.rank, c.rank, d.rank };
  return orient3d_sos(p, r);
}

// ray_intersects_triangle
bool ray_crosses(const Geometrium& g, RankedPoint D, RankedPoint E, FaceId f)
{
  const auto v = g.vertices_of(f);
  const RankedPoint A{ g.vertex(v[0]).p, v[0].i }, B{ g.vertex(v[1]).p, v[1].i }, C{ g.vertex(v[2]).p, v[2].i };

  int ABCE = orient_sos(A, B, C, E);
  int ABCD = orient_sos(A, B, C, D);
  if (ABCE < 0 || ABCD > 0)
  {
    std::swap(D, E);
    std::swap(ABCE, ABCD);
  }
  if (ABCE < 0 || ABCD > 0) return false;
  if (orient_sos(A, D, C, E) < 0) return false;
  if (orient_sos(A, B, D, E) < 0) return false;
  if (orient_sos(B, C, D, E) < 0) return false;
  return true;
}

} // namespace

bool is_inside(const BBTree& tree, Vec3 p, bool is_open)
{
  if (tree.empty())
    throw std::invalid_argument("gts::is_inside: the tree is empty");
  const BBox& bb = tree.box();
  const RankedPoint rp{ p, UINT_MAX - 1 };
  // The ray ends beyond the box, also when max.x is 0.
  const double step = std::max(std::fabs(bb.max.x), bb.max.x - bb.min.x);
  const RankedPoint rp1{ { bb.max.x + (step > 0. ? step/10. : 1.), p.y, p.z }, UINT_MAX };
  unsigned nc = 0;
  for (const BBLeaf& l : tree.stabbed(p))
    if (ray_crosses(tree.geometrium(), rp, rp1, l.face))
      nc++;
  return is_open ? (nc % 2 == 0) : (nc % 2 != 0);
}

//==============================================================================
// Sampled distances
//==============================================================================

Range face_distance(const BBTree& tree, const Geometrium& g, FaceId f, double delta, const BBTree::DistanceFunc& d)
{
  if ( ! (delta > 0.))
    throw std::invalid_argument("gts::face_distance: delta must be positive");
  const Triangle tr = g.triangle(f);
  const Vec3 p1 = tr.a, p2 = tr.b;
  const Vec3 p1p2 = tr.b - tr.a, p1p3 = tr.c - tr.a;
  Range range;

  const double l1 = std::sqrt(dot(p1p2, p1p2));
  const unsigned n1 = l1/delta + 1;
  const double dt1 = 1.0/double(n1);
  double t1 = 0.0;
  for (unsigned i = 0; i <= n1; i++, t1 += dt1)
  {
    double t2 = 1. - t1;
    double x = t2*p1p3.x;
    double y = t2*p1p3.y;
    double z = t2*p1p3.z;
    const double l2 = std::sqrt(x*x + y*y + z*z);
    const unsigned n2 = unsigned(l2/delta + 1);
    const double dt2 = t2/double(n2);

    x = t2*p1.x + t1*p2.x;
    y = t2*p1.y + t1*p2.y;
    z = t2*p1.z + t1*p2.z;

    t2 = 0.0;
    for (unsigned j = 0; j <= n2; j++, t2 += dt2)
      range.add(tree.point_distance({ x + t2*p1p3.x, y + t2*p1p3.y, z + t2*p1p3.z }, d));
  }
  range.update();
  return range;
}

Range edge_distance(const BBTree& tree, const Geometrium& g, EdgeId e, double delta, const BBTree::DistanceFunc& d)
{
  if ( ! (delta > 0.))
    throw std::invalid_argument("gts::edge_distance: delta must be positive");
  const Segment s = g.segment(e);
  const Vec3 p1 = s.a, p1p2 = s.b - s.a;
  Range range;

  const double l = std::sqrt(dot(p1p2, p1p2));
  const unsigned n = unsigned(l/delta + 1);
  const double dt = 1.0/double(n);
  double t = 0.0;
  for (unsigned i = 0; i <= n; i++, t += dt)
    range.add(tree.point_distance({ p1.x + t*p1p2.x, p1.y + t*p1p2.y, p1.z + t*p1p2.z }, d));
  range.update();
  return range;
}

namespace {

// The weighted mean and standard deviation of a surface's distances.
void finish(Range& range, double total)
{
  if (total > 0.)
  {
    const double v = range.sum2 - range.sum*range.sum/total;
    range.stddev = v >= 0. ? std::sqrt(v/total) : 0.;
    range.mean = range.sum/total;
  }
  else
    range.min = range.max = range.mean = range.stddev = 0.;
}

void check_fraction(double delta, const char* who)
{
  if ( ! (delta > 0. && delta < 1.))
    throw std::invalid_argument(std::string("gts::") + who + ": delta must be between 0 and 1");
}

} // namespace

Range surface_distance(const BBTree& tree, const Surface& s, double delta, const BBTree::DistanceFunc& d)
{
  check_fraction(delta, "surface_distance");
  const Geometrium& g = s.geometrium();
  Range range;
  delta *= std::sqrt(diagonal2(tree.box()));
  double total_area = 0.;
  for (FaceId f : s.faces())
  {
    const Range rt = face_distance(tree, g, f, delta, d);
    if (rt.min < range.min) range.min = rt.min;
    if (rt.max > range.max) range.max = rt.max;
    range.n += rt.n;
    const double area = g.area(f);
    total_area += area;
    range.sum += area*rt.mean;
    range.sum2 += area*rt.mean*rt.mean;
  }
  finish(range, total_area);
  return range;
}

Range boundary_distance(const BBTree& tree, const Surface& s, double delta, const BBTree::DistanceFunc& d)
{
  check_fraction(delta, "boundary_distance");
  const Geometrium& g = s.geometrium();
  Range range;
  delta *= std::sqrt(diagonal2(tree.box()));
  double total_length = 0.;
  for (EdgeId e : s.edges())
  {
    if ( ! boundary_face(g, e, &s).valid()) continue;
    const Range re = edge_distance(tree, g, e, delta, d);
    if (re.min < range.min) range.min = re.min;
    if (re.max > range.max) range.max = re.max;
    range.n += re.n;
    const Segment sg = g.segment(e);
    const double length = gts::distance(sg.a, sg.b);
    total_length += length;
    range.sum += length*re.mean;
    range.sum2 += length*re.mean*re.mean;
  }
  finish(range, total_length);
  return range;
}

SurfaceDistance distance(const Surface& s1, const Surface& s2, double delta)
{
  check_fraction(delta, "distance");
  const Range zero{ 0., 0., 0., 0., 0., 0., 0 };
  SurfaceDistance r{ zero, zero };
  if (s2.face_number() == 0) return r;

  const BBTree faces(s2);
  r.faces = surface_distance(faces, s1, delta);

  const Geometrium& g2 = s2.geometrium();
  std::vector<EdgeId> edges;
  for (EdgeId e : s2.edges())
    if (boundary_face(g2, e, &s2).valid()) edges.push_back(e);
  if ( ! edges.empty())
  {
    const BBTree boundary(g2, edges);
    r.boundary = boundary_distance(boundary, s1, delta);
  }
  return r;
}

} // namespace gts
