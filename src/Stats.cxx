// Copyright (C) 1999 Stéphane Popinet (GTS: surface.c, misc.c, face.c, edge.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/Stats.h"

#include <algorithm>
#include <cmath>
#include <format>

namespace gts {

void Range::update()
{
  if (n > 0)
  {
    const double v = sum2 - sum * sum / double(n);
    stddev = v >= 0. ? std::sqrt(v / double(n)) : 0.;
    mean = sum / double(n);
  }
  else
  {
    min = max = mean = stddev = 0.;
  }
}

double area(const Surface& s)
{
  const Geometrium& g = s.geometrium();
  double a = 0.0;
  for (FaceId f : s.faces()) a += g.area(f);
  return a;
}

double volume(const Surface& s)
{
  const Geometrium& g = s.geometrium();
  double v = 0.0;
  for (FaceId f : s.faces())
  {
    const Triangle t = g.triangle(f);
    const Vec3 &pa = t.a, &pb = t.b, &pc = t.c;
    v += (pa.x * (pb.y * pc.z - pb.z * pc.y) +
          pb.x * (pc.y * pa.z - pc.z * pa.y) +
          pc.x * (pa.y * pb.z - pa.z * pb.y));
  }
  return v / 6.;
}

unsigned face_number(const Surface& s, EdgeId e)
{
  unsigned n = 0;
  for (FaceId f : s.geometrium().faces_of(e))
    if (s.contains(f)) ++n;
  return n;
}

FaceId boundary_face(const Geometrium& g, EdgeId e, const Surface* s)
{
  FaceId f;
  for (FaceId t : g.faces_of(e))
  {
    if (s && ! s->contains(t)) continue;
    if (f.valid()) return FaceId();
    f = t;
  }
  return f;
}

bool is_boundary(const Geometrium& g, VertexId v, const Surface* s)
{
  for (EdgeId e : g.edges_of(v))
    if (boundary_face(g, e, s).valid()) return true;
  return false;
}

std::vector<EdgeId> boundary(const Surface& s)
{
  std::vector<EdgeId> list;
  for (EdgeId e : s.edges())
    if (boundary_face(s.geometrium(), e, &s).valid()) list.push_back(e);
  std::reverse(list.begin(), list.end());
  return list;
}

std::size_t remove_duplicate_faces(Surface& s)
{
  const Geometrium& g = s.geometrium();
  // build_list prepended the faces in the order of GTS's foreach.
  std::vector<FaceId> faces(s.faces().begin(), s.faces().end());
  std::reverse(faces.begin(), faces.end());
  std::size_t n = 0;
  for (FaceId f : faces)
  {
    const auto& e = g.face(f).e;
    for (FaceId t : g.faces_of(e[0]))
    {
      if (t == f || ! s.contains(t)) continue;
      const auto& x = g.face(t).e;
      auto has = [&](EdgeId y) { return x[0] == y || x[1] == y || x[2] == y; };
      if (has(e[1]) && has(e[2]))
      {
        s.remove_face(f);
        ++n;
        break;
      }
    }
  }
  return n;
}

CenterOfMass center_of_mass(const Surface& s)
{
  const Geometrium& g = s.geometrium();
  double volume = 0.;
  double cm[3] = { 0., 0., 0. };
  for (FaceId f : s.faces())
  {
    const Triangle t = g.triangle(f);
    const Vec3 &p1 = t.a, &p2 = t.b, &p3 = t.c;
    const double x1 = p2.x - p1.x, y1 = p2.y - p1.y, z1 = p2.z - p1.z;
    const double x2 = p3.x - p1.x, y2 = p3.y - p1.y, z2 = p3.z - p1.z;
    const double nx = y1*z2 - z1*y2, ny = z1*x2 - x1*z2, nz = x1*y2 - y1*x2;
    cm[0] += nx*(p1.x*p1.x + p2.x*p2.x + p3.x*p3.x + p1.x*p2.x + p1.x*p3.x + p2.x*p3.x);
    cm[1] += ny*(p1.y*p1.y + p2.y*p2.y + p3.y*p3.y + p1.y*p2.y + p1.y*p3.y + p2.y*p3.y);
    cm[2] += nz*(p1.z*p1.z + p2.z*p2.z + p3.z*p3.z + p1.z*p2.z + p1.z*p3.z + p2.z*p3.z);
    volume += nx*(p1.x + p2.x + p3.x);
  }
  if (volume != 0.)
  {
    cm[0] /= 4.*volume;
    cm[1] /= 4.*volume;
    cm[2] /= 4.*volume;
  }
  return { { cm[0], cm[1], cm[2] }, volume/6. };
}

CenterOfArea center_of_area(const Surface& s)
{
  const Geometrium& g = s.geometrium();
  double area = 0.;
  double cm[3] = { 0., 0., 0. };
  for (FaceId f : s.faces())
  {
    const Triangle t = g.triangle(f);
    const double a = g.area(f);
    cm[0] += a*(t.a.x + t.b.x + t.c.x);
    cm[1] += a*(t.a.y + t.b.y + t.c.y);
    cm[2] += a*(t.a.z + t.b.z + t.c.z);
    area += a;
  }
  if (area != 0.)
  {
    cm[0] /= 3.*area;
    cm[1] /= 3.*area;
    cm[2] /= 3.*area;
  }
  return { { cm[0], cm[1], cm[2] }, area };
}

bool is_manifold(const Surface& s)
{
  for (EdgeId e : s.edges())
    if (face_number(s, e) > 2) return false;
  return true;
}

bool is_closed(const Surface& s)
{
  for (EdgeId e : s.edges())
    if (face_number(s, e) != 2) return false;
  return true;
}

bool is_orientable(const Surface& s)
{
  const Geometrium& g = s.geometrium();
  for (EdgeId e : s.edges())
  {
    FaceId f1, f2;
    for (FaceId f : g.faces_of(e))
    {
      if ( ! s.contains(f)) continue;
      if      ( ! f1.valid()) f1 = f;
      else if ( ! f2.valid()) f2 = f;
      else return false;
    }
    if (f1.valid() && f2.valid() && ! g.compatible(f1, f2, e))
      return false;
  }
  return true;
}

namespace {

// triangle_is_incompatible in face.c
bool incompatible_on(const Surface& s, FaceId t, EdgeId e)
{
  const Geometrium& g = s.geometrium();
  for (FaceId f : g.faces_of(e))
    if (f != t && s.contains(f) && ! g.compatible(t, f, e))
      return true;
  return false;
}

} // namespace

SurfaceStats stats(const Surface& s)
{
  const Geometrium& g = s.geometrium();
  SurfaceStats st;

  for (VertexId v : s.vertices())
  {
    unsigned n = 0;
    for (EdgeId e : g.edges_of(v))
      if (face_number(s, e) > 0) ++n; // gts_edge_has_parent_surface
    st.edges_per_vertex.add(n);
  }
  for (EdgeId e : s.edges())
  {
    const unsigned nt = face_number(s, e);
    if (g.duplicate_of(e).valid()) ++st.n_duplicate_edges;
    if (nt == 1)      ++st.n_boundary_edges;
    else if (nt > 2)  ++st.n_non_manifold_edges;
    st.faces_per_edge.add(nt);
  }
  for (FaceId f : s.faces())
  {
    const auto& e = g.face(f).e;
    if (incompatible_on(s, f, e[0]) || incompatible_on(s, f, e[1]) || incompatible_on(s, f, e[2]))
      ++st.n_incompatible_faces;
    if (g.duplicate_of(f).valid()) ++st.n_duplicate_faces;
    ++st.n_faces;
  }

  st.edges_per_vertex.update();
  st.faces_per_edge.update();
  return st;
}

QualityStats quality_stats(const Surface& s)
{
  const Geometrium& g = s.geometrium();
  QualityStats q;

  for (EdgeId e : s.edges())
  {
    const Segment sg = g.segment(e);
    q.edge_length.add(distance(sg.a, sg.b));
    // All faces on the edge, as GTS takes e->triangles, not only those of s.
    std::vector<FaceId> faces;
    for (FaceId f : g.faces_of(e)) faces.push_back(f);
    for (std::size_t i = 0; i < faces.size(); ++i)
      for (std::size_t j = i + 1; j < faces.size(); ++j)
        q.edge_angle.add(std::fabs(g.angle(faces[i], faces[j])));
  }
  for (FaceId f : s.faces())
  {
    q.face_quality.add(g.quality(f));
    q.face_area.add(g.area(f));
  }

  q.face_quality.update();
  q.face_area.update();
  q.edge_length.update();
  q.edge_angle.update();
  return q;
}

std::string to_string(const Range& r)
{
  return std::format("min: {:g} mean: {:g} | {:g} max: {:g}", r.min, r.mean, r.stddev, r.max);
}

void print_stats(const Surface& s, std::ostream& out)
{
  const SurfaceStats st = stats(s);
  const QualityStats qs = quality_stats(s);
  out << std::format("# vertices: {} edges: {} faces: {}\n"
                     "# Connectivity statistics\n"
                     "#   incompatible faces: {}\n"
                     "#   duplicate faces: {}\n"
                     "#   boundary edges: {}\n"
                     "#   duplicate edges: {}\n"
                     "#   non-manifold edges: {}\n",
                     st.edges_per_vertex.n, st.faces_per_edge.n, st.n_faces,
                     st.n_incompatible_faces, st.n_duplicate_faces, st.n_boundary_edges,
                     st.n_duplicate_edges, st.n_non_manifold_edges);
  out << "#   edges per vertex: " << to_string(st.edges_per_vertex)
      << "\n#   faces per edge: " << to_string(st.faces_per_edge)
      << "\n# Geometric statistics\n#   face quality: " << to_string(qs.face_quality)
      << "\n#   face area  : " << to_string(qs.face_area)
      << "\n#   edge length : " << to_string(qs.edge_length) << '\n';
}

} // namespace gts
