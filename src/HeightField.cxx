// Copyright (C) 1999 Stéphane Popinet (GTS: examples/happrox.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-13.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/HeightField.h"

#include "gts/Delaunay.h"
#include "gts/EHeap.h"
#include "gts/FacePoints.h"
#include "gts/Stats.h"
#include "gts/Topology.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace gts {

double height_cost(Vec3 plane, Vec3 p)
{
  return std::fabs(p.x*plane.x + p.y*plane.y + plane.z - p.z);
}

HeightCost relative_height_cost(double z0)
{
  return [z0](Vec3 plane, Vec3 p) {
    const double dz = p.x*plane.x + p.y*plane.y + plane.z - p.z;
    return std::fabs(p.z) < z0 ? std::fabs(dz/z0) : std::fabs(dz/p.z);
  };
}

HeightStop stop_at_vertices(std::size_t n)
{
  return [n](double, std::size_t number) { return number >= n; };
}

HeightStop stop_below_cost(double c)
{
  return [c](double cost, std::size_t) { return cost < c; };
}

namespace {

// triangle_plane: (a, b, c) of the plane z = a x + b y + c through the
// vertices of f's first edge and the vertex opposite it.
Vec3 face_plane(const Geometrium& g, FaceId f, const char* who)
{
  const Edge& e1 = g.edge(g.face(f).e[0]);
  const Vec3 v1 = g.vertex(e1.v[0]).p, v2 = g.vertex(e1.v[1]).p, v3 = g.vertex(g.opposite_vertex(f)).p;
  const double x1 = v2.x - v1.x, y1 = v2.y - v1.y;
  const double x2 = v3.x - v1.x, y2 = v3.y - v1.y;
  const double det = x1*y2 - x2*y1;
  if (det == 0.) throw std::runtime_error(std::string("gts::") + who + ": a face is degenerate in the xy plane");
  return { (y2*(v2.z - v1.z) - y1*(v3.z - v1.z))/det,
           (-x2*(v2.z - v1.z) + x1*(v3.z - v1.z))/det,
           ((- v1.x*y2 + v1.y*x2)*(v2.z - v1.z) + (- v1.y*x1 + v1.x*y1)*(v3.z - v1.z))/det + v1.z };
}

// surface_hf_refine.
void refine(Surface& s, FacePoints& fp, const HeightCost& cost, const HeightStop& stop, const char* who)
{
  Geometrium& g = s.geometrium();
  using Heap = EHeap<FaceId>;
  Heap heap;
  constexpr Heap::Handle kNoHandle = ~Heap::Handle(0);
  std::vector<Heap::Handle> handle;      // by face index
  std::vector<FacePoints::PointId> best; // by face index
  auto handle_of = [&](FaceId f) -> Heap::Handle& {
    if (handle.size() <= f.i) handle.resize(f.i + 1, kNoHandle);
    return handle[f.i];
  };
  auto best_of = [&](FaceId f) -> FacePoints::PointId& {
    if (best.size() <= f.i) best.resize(f.i + 1);
    return best[f.i];
  };
  // A face that leaves s leaves the heap, as list_face_destroy.
  const Subscription sub = s.on_face([&](const FaceEvent& ev) {
    if (ev.added) return;
    Heap::Handle& h = handle_of(ev.face);
    if (h != kNoHandle) { heap.remove(h); h = kNoHandle; }
  });
  // list_face_update: a face without points is left as it is.
  auto update = [&](FaceId f) {
    bool any = false;
    Vec3 plane;
    double maxerr = 0.;
    for (FacePoints::PointId i : fp.points(f))
    {
      if ( ! any) { plane = face_plane(g, f, who); any = true; }
      const double e = cost(plane, fp.point(i));
      if (e > maxerr) { maxerr = e; best_of(f) = i; }
    }
    if ( ! any) return;
    Heap::Handle& h = handle_of(f);
    if (h != kNoHandle) { heap.remove(h); h = kNoHandle; }
    if (maxerr > 0.) h = heap.insert(f, - maxerr);
  };

  for (FaceId f : s.faces()) update(f);
  std::size_t nv = 4;
  while (const auto top = heap.remove_top())
  {
    const FaceId f = top->first;
    handle_of(f) = kNoHandle;
    if (stop(- top->second, nv)) break;
    const FacePoints::PointId b = best_of(f);
    const Vec3 p = fp.point(b);
    // gts_delaunay_add_vertex_to_face refuses a point on a vertex of f; the
    // point stays, and f is out of the heap.
    bool duplicate = false;
    for (VertexId w : g.vertices_of(f))
      if (g.vertex(w).p.x == p.x && g.vertex(w).p.y == p.y) duplicate = true;
    if ( ! duplicate)
    {
      fp.remove(b);
      const VertexId v = g.add_vertex(p);
      delaunay_add_vertex_to_face(s, v, f);
      for (FaceId t : vertex_faces(g, v)) update(t);
    }
    ++nv;
  }
}

void destroy_vertices(Geometrium& g, const std::array<VertexId, 3>& vs)
{
  KeepOrphans keep(g, { .vertices = true });
  for (VertexId v : vs) g.destroy(v);
}

} // namespace

void approximate_grid(Surface& s, std::size_t width, std::size_t height, std::span<const double> z,
                      const HeightStop& stop, const HeightCost& cost)
{
  if (width < 2 || height < 2 || z.size() != width * height)
    throw std::invalid_argument("gts::approximate_grid: the grid needs at least 2 x 2 heights, width * height of them");
  Geometrium& g = s.geometrium();
  auto point = [&](std::size_t i, std::size_t j) { return Vec3{ double(i), double(j), z[j * width + i] }; };
  const std::array<Vec3, 4> corners = { point(0, 0), point(0, height - 1), point(width - 1, 0),
                                        point(width - 1, height - 1) };
  std::array<VertexId, 4> cv;
  for (int k = 0; k < 4; ++k) cv[k] = g.add_vertex(corners[k]);

  // gts_triangle_enclosing takes the corners from a list built by prepending.
  const std::array<Vec3, 4> prepended = { corners[3], corners[2], corners[1], corners[0] };
  const FaceId t = enclosing_triangle(g, prepended, 100.);
  const auto hull = g.vertices_of(t);
  s.add_face(t);

  // The points but the corners, in the order of GTS's list, which was built
  // by prepending.
  std::vector<Vec3> points;
  for (std::size_t i = 1; i + 1 < width; i++)
  {
    for (std::size_t j = 1; j + 1 < height; j++) points.push_back(point(i, j));
    points.push_back(point(i, 0));
    points.push_back(point(i, height - 1));
  }
  for (std::size_t j = 1; j + 1 < height; j++)
  {
    points.push_back(point(0, j));
    points.push_back(point(width - 1, j));
  }
  FacePoints fp(s);
  for (auto p = points.rbegin(); p != points.rend(); ++p) fp.add(*p, t);

  for (int k = 0; k < 4; ++k)
  {
    const FaceId f = k == 0 ? t : locate(s, corners[k]);
    if ( ! f.valid() || delaunay_add_vertex_to_face(s, cv[k], f).valid())
      throw std::logic_error("gts::approximate_grid: a corner could not be inserted");
  }
  refine(s, fp, cost, stop, "approximate_grid");
  destroy_vertices(g, hull);
}

void approximate_points(Surface& s, std::span<const Vec3> points, const HeightStop& stop,
                        const HeightCost& cost, Enclosing enclosing)
{
  if (points.empty()) throw std::invalid_argument("gts::approximate_points: no points");
  Geometrium& g = s.geometrium();
  const FaceId t = enclosing_triangle(g, points, 10.);
  const auto hull = g.vertices_of(t);
  for (VertexId w : hull) g.position(w).z = enclosing == Enclosing::Remove ? -1e30 : -10.;
  s.add_face(t);
  FacePoints fp(s);
  for (const Vec3& p : points) fp.add(p, t);

  refine(s, fp, cost, stop, "approximate_points");
  if (enclosing == Enclosing::Remove)
    destroy_vertices(g, hull);
  else if (enclosing == Enclosing::Close)
  {
    const std::vector<EdgeId> l = boundary(s);
    if (l.size() != 3) throw std::logic_error("gts::approximate_points: the boundary is not the enclosing triangle");
    const FaceId f = g.add_face(l[0], l[1], l[2]);
    s.add_face(f);
    if ( ! is_compatible(s, f)) g.revert_face(f);
  }
}

} // namespace gts
