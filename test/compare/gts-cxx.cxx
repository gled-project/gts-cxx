// gts-cxx: the twin of gts-ref, running the same operations through the C++
// library and printing the same summary. Operations not ported yet print
// "status unsupported" and exit with 3.
//
// Copyright (C) 2026 Matevz Tadel
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/all.h"

#include <fstream>
#include <limits>
#include <array>
#include <map>
#include <memory>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <numbers>
#include <stdexcept>
#include <print>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

using namespace gts;

namespace {

[[noreturn]] void fail(std::string_view msg)
{
  std::println("status error {}", msg);
  std::exit(1);
}

[[noreturn]] void unsupported(std::string_view cmd)
{
  std::println("status unsupported {}", cmd);
  std::exit(3);
}

Surface& read_surface(Geometrium& g, const char* path)
{
  auto r = read_gts(g, std::filesystem::path(path));
  if ( ! r) fail(std::format("cannot read {}: {}", path, r.error().message));
  return **r;
}

void print_range(const char* name, const Range& r)
{
  std::println("{} n={} min={:.12e} mean={:.12e} max={:.12e}", name, r.n, r.min, r.mean, r.max);
}

// A dot file of a graph: its nodes and edges, and the range of the degrees
// of its nodes, numbered "nK" from 1.
void print_dot(std::string_view key, const std::string& dot)
{
  std::vector<unsigned> degree;
  std::size_t edges = 0;
  for (std::size_t b = 0, e; b < dot.size(); b = e + 1)
  {
    e = dot.find('\n', b);
    if (e == std::string::npos) e = dot.size();
    const std::string_view line(dot.data() + b, e - b);
    if ( ! line.starts_with("  n")) continue;
    const unsigned n1 = unsigned(std::strtoul(line.data() + 3, nullptr, 10));
    const std::size_t arrow = line.find(" -> n");
    if (arrow == std::string_view::npos)
    {
      if (degree.size() <= n1) degree.resize(n1 + 1);
      continue;
    }
    const unsigned n2 = unsigned(std::strtoul(line.data() + arrow + 5, nullptr, 10));
    if (degree.size() <= std::max(n1, n2)) degree.resize(std::max(n1, n2) + 1);
    ++degree[n1];
    ++degree[n2];
    ++edges;
  }
  unsigned dmin = degree.size() > 1 ? ~0u : 0, dmax = 0, sum = 0;
  for (std::size_t i = 1; i < degree.size(); ++i)
  {
    dmin = std::min(dmin, degree[i]);
    dmax = std::max(dmax, degree[i]);
    sum += degree[i];
  }
  std::println("{}_dot nodes={} edges={} degree_min={} degree_max={} degree_sum={} closed={}", key,
               degree.empty() ? 0 : degree.size() - 1, edges, dmin, dmax, sum, dot.ends_with("}\n") ? 1 : 0);
}

void print_summary(const Surface& s)
{
  const Geometrium& g = s.geometrium();
  std::println("vertices {}", s.vertex_number());
  std::println("edges {}",    s.edge_number());
  std::println("faces {}",    s.face_number());
  if (s.face_number() == 0)
    return;

  std::println("area {:.12e}",   area(s));
  std::println("volume {:.12e}", volume(s));
  std::println("manifold {}",    is_manifold(s)   ? 1 : 0);
  std::println("closed {}",      is_closed(s)     ? 1 : 0);
  std::println("orientable {}",  is_orientable(s) ? 1 : 0);
  std::size_t self_ties = 0;
  std::println("self_intersecting {}", is_self_intersecting(s, &self_ties) ? 1 : 0);
  // Only from gts-cxx: compare.py treats a differing self_intersecting with
  // ties as degenerate.
  std::println("self_sos_ties {}", self_ties);

  BBox b;
  for (VertexId v : s.vertices()) b.add(g.vertex(v).p);
  std::println("bbox {:.12e} {:.12e} {:.12e} {:.12e} {:.12e} {:.12e}", b.min.x, b.min.y, b.min.z, b.max.x, b.max.y, b.max.z);

  const SurfaceStats st = stats(s);
  std::println("incompatible_faces {}", st.n_incompatible_faces);
  std::println("duplicate_faces {}",    st.n_duplicate_faces);
  std::println("duplicate_edges {}",    st.n_duplicate_edges);
  std::println("boundary_edges {}",     st.n_boundary_edges);
  std::println("non_manifold_edges {}", st.n_non_manifold_edges);
  print_range("edges_per_vertex", st.edges_per_vertex);
  print_range("faces_per_edge",   st.faces_per_edge);

  const QualityStats qs = quality_stats(s);
  print_range("face_quality", qs.face_quality);
  print_range("face_area",    qs.face_area);
  print_range("edge_length",  qs.edge_length);
  print_range("edge_angle",   qs.edge_angle);
}

//==============================================================================
// Refinement, as in GledGTS's GTSRetriangulator
//==============================================================================

// Cost from the angle opposite to the edge, from GTSRetriangulator.cxx:
// 1 for 60 degrees, falling linearly to 0 at 0 and 180 degrees. The cost of
// an edge is the smallest over its triangles.
double edge_angle(const Geometrium& g, FaceId t, EdgeId e)
{
  static const double r60 = std::numbers::pi / 3;
  const VertexId v = g.opposite_vertex(t, e);
  if ( ! v.valid()) return 1.0;
  const Vec3 p = g.vertex(v).p;
  const Vec3 a = g.vertex(g.edge(e).v[0]).p - p, b = g.vertex(g.edge(e).v[1]).p - p;
  const double phi = std::atan2(norm(cross(a, b)), dot(a, b));
  return phi <= r60 ? phi / r60 : 1.5 - 0.5 * phi / r60;
}

double cost_angle(const Geometrium& g, EdgeId e)
{
  double cost = 1.0;
  for (FaceId t : g.faces_of(e))
    cost = std::fmin(edge_angle(g, t, e), cost);
  return cost;
}

// refine IN COST STOP VALUE
void refine_cmd(Surface& s, std::span<char*> a)
{
  const std::string_view c = a[0], st = a[1];
  EdgeCost cost;
  if      (c == "length") cost = {};
  else if (c == "angle")  cost = cost_angle;
  else if (c == "volume") unsupported("refine volume");
  else fail("COST must be length, volume or angle");

  RefineStop stop;
  if (st == "number")
  {
    const std::size_t max = std::strtoul(a[2], nullptr, 10);
    stop = [max](double, std::size_t n) { return n > max; };
  }
  else if (st == "cost")
  {
    double min = std::strtod(a[2], nullptr);
    // The length cost works on the square of the edge length.
    if (c == "length") min *= min;
    stop = [min](double cost, std::size_t) { return cost < min; };
  }
  else fail("STOP must be number or cost");

  refine(s, stop, cost);
}

// coarsen IN COST STOP VALUE [MIDVERT [MINANGLE_DEG]], as gts-ref's
// parse_retriangulate and coarsen, with GledGTS's defaults.
void coarsen_cmd(Surface& s, std::span<char*> a)
{
  const std::string_view c = a[0], st = a[1];
  const VolumeOptimizedParams vo = { 0.5, 0.5, 0 };
  EdgeCost cost;
  if      (c == "length") cost = {};
  else if (c == "volume") cost = [vo](const Geometrium& g, EdgeId e) { return volume_optimized_cost(g, e, vo); };
  else if (c == "angle")  cost = cost_angle;
  else fail("COST must be length, volume or angle");

  CoarsenStop stop;
  if (st == "number")
    stop = stop_below_edges(std::strtoul(a[2], nullptr, 10));
  else if (st == "cost")
  {
    double max = std::strtod(a[2], nullptr);
    // The length cost works on the square of the edge length.
    if (c == "length") max *= max;
    stop = stop_above_cost(max);
  }
  else fail("STOP must be number or cost");

  bool midvert_volume = true;
  if (a.size() > 3)
  {
    const std::string_view m = a[3];
    if      (m == "mid")    midvert_volume = false;
    else if (m == "volume") midvert_volume = true;
    else fail("MIDVERT must be mid or volume");
  }
  const double min_angle_deg = a.size() > 4 ? std::strtod(a[4], nullptr) : 1.0;

  CollapsePosition position;
  if (midvert_volume)
    position = [vo](const Geometrium& g, EdgeId e) { return volume_optimized_vertex(g, e, vo); };
  coarsen(s, stop, cost, position, min_angle_deg * std::numbers::pi / 180);
}

//==============================================================================
// Iso-surfaces, as in GTSIsoMaker::MakeSurface and gts-ref
//==============================================================================

double field_sphere(Vec3 p) { return p.x*p.x + p.y*p.y + p.z*p.z; }

double field_torus(Vec3 p)
{
  // Major radius 0.7, minor radius 0.3; the surface is at 0.09.
  const double q = std::sqrt(p.x*p.x + p.y*p.y) - 0.7;
  return q*q + p.z*p.z;
}

double field_blobs(Vec3 p)
{
  // Two blobs, 1/r^2 potentials summed; the surface is at 4.
  auto f = [](double dx, double dy, double dz) { return 1.0 / (dx*dx + dy*dy + dz*dz + 1e-12); };
  return f(p.x - 0.4, p.y, p.z) + f(p.x + 0.4, p.y - 0.1, p.z + 0.05);
}

// iso ALGO FUNC N
void iso_cmd(Surface& s, std::span<char*> a)
{
  const std::string_view algo = a[0], func = a[1];
  double (*f)(Vec3);
  double iso;
  Vec3 inside; // across the surface from (1.2, 0.1, 0.05), for seeded
  if      (func == "sphere") { f = field_sphere; iso = 1;    inside = { 0, 0, 0 };   }
  else if (func == "torus")  { f = field_torus;  iso = 0.09; inside = { 0.7, 0, 0 }; }
  else if (func == "blobs")  { f = field_blobs;  iso = 4;    inside = { 0.4, 0, 0 }; }
  else fail("FUNC must be sphere, torus or blobs");

  const unsigned n = unsigned(std::strtoul(a[2], nullptr, 10));
  const double lo = -1.25, d = 2.5 / n;
  const CartesianGrid g = { n + 1, n + 1, n + 1, lo, d, lo, d, lo, d };
  const PlaneFunc pf = sample(f);
  if      (algo == "cartesian")     isosurface_cartesian    (s, g, pf, iso);
  else if (algo == "seeded")        isosurface_cartesian    (s, g, f, iso, cells_crossing(g, f, iso, inside, { 1.2, 0.1, 0.05 }));
  else if (algo == "tetra")         isosurface_tetra        (s, g, pf, iso);
  else if (algo == "tetra_bounded") isosurface_tetra_bounded(s, g, pf, iso);
  else if (algo == "tetra_bcl")     isosurface_tetra_bcl    (s, g, pf, iso);
  else fail("ALGO must be cartesian, seeded, tetra, tetra_bounded or tetra_bcl");
}

// isoslice FUNC N, as gts-ref.
Vec3 warp_point(unsigned i, unsigned j, unsigned k, double lo, double d)
{
  return { lo + d*i + 0.15*d*std::sin(1.7*k + 0.3*j),
           lo + d*j + 0.15*d*std::cos(2.3*i),
           lo + d*k + 0.15*d*std::sin(1.1*j + 0.5*i) };
}

void isoslice_cmd(Surface& s, std::span<char*> a)
{
  const std::string_view func = a[0];
  double (*f)(Vec3);
  double iso;
  if      (func == "sphere") { f = field_sphere; iso = 1;    }
  else if (func == "torus")  { f = field_torus;  iso = 0.09; }
  else if (func == "blobs")  { f = field_blobs;  iso = 4;    }
  else fail("FUNC must be sphere, torus or blobs");

  const unsigned n = unsigned(std::strtoul(a[1], nullptr, 10)), m = n + 1;
  const double lo = -1.25, d = 2.5 / n;
  GridPlane planes[2] = { GridPlane(m, m), GridPlane(m, m) };
  PlaneValues values[2] = { PlaneValues(m, m), PlaneValues(m, m) };
  auto fill_plane = [&](int w, unsigned k) {
    for (unsigned i = 0; i < m; ++i)
      for (unsigned j = 0; j < m; ++j)
      {
        const Vec3 p = warp_point(i, j, k, lo, d);
        planes[w](i, j) = p;
        values[w](i, j) = f(p);
      }
  };
  IsoSlice sa(m, m), sb(m, m);
  IsoSlice* slices[2] = { &sa, &sb };
  int x = 0, y = 1;
  fill_plane(x, 0);
  fill_plane(y, 1);
  slices[0]->fill(planes[x], &planes[y], values[x], &values[y], iso);
  for (unsigned k = 2; k < m; ++k)
  {
    std::swap(x, y);
    fill_plane(y, k);
    slices[1]->fill(planes[x], &planes[y], values[x], &values[y], iso);
    isosurface_slice(s, *slices[0], *slices[1]);
    std::swap(slices[0], slices[1]);
  }
  slices[1]->fill(planes[y], nullptr, values[y], nullptr, iso);
  isosurface_slice(s, *slices[0], *slices[1]);
}

//==============================================================================
// Booleans, as in GTSBoolOpHelper with post-processing off, and gts-ref
//==============================================================================

// bool OP A B: the result in a new surface of g; a and b are removed.
Surface& make_boolean(Geometrium& g, std::string_view op, Surface& a, Surface& b)
{
  if (op != "union" && op != "inter" && op != "diff")
    fail("OP must be union, inter or diff");

  Surface& r = g.new_surface();
  {
    SurfaceInter si(a, b);
    const SurfaceInter::Check c = si.check();
    std::println("inter_edges {}",      si.edges().size());
    std::println("inter_orientable {}", c.ok     ? 1 : 0);
    std::println("inter_closed {}",     c.closed ? 1 : 0);
    // Only from gts-cxx: compare.py treats a case with ties as degenerate.
    std::println("sos_ties {}",         si.sos_ties());
    if ( ! c.ok || ! c.closed)
      fail("intersection curve not orientable or not closed");

    if (op == "union")
    {
      si.boolean(r, BooleanOp::A_out_B);
      si.boolean(r, BooleanOp::B_out_A);
    }
    else if (op == "inter")
    {
      si.boolean(r, BooleanOp::A_in_B);
      si.boolean(r, BooleanOp::B_in_A);
    }
    else
    {
      si.boolean(r, BooleanOp::A_out_B);
      si.boolean(r, BooleanOp::B_in_A);
      for (FaceId f : si.pieces_b().faces()) g.revert_face(f);
      for (FaceId f : b.faces())             g.revert_face(f);
    }
  }
  g.remove_surface(a);
  g.remove_surface(b);
  return r;
}

//==============================================================================
// Cluster grid, as in GTSRetriangulator::OutOfCoreSimplification and gts-ref
//==============================================================================

// cluster IN DELTA: the result in a new surface of g; in is removed.
Surface& make_cluster(Geometrium& g, Surface& in, double delta)
{
  Surface& d = g.new_surface();
  BBox b;
  for (VertexId v : in.vertices()) b.add(g.vertex(v).p);
  {
    ClusterGrid grid(d, b, delta);
    for (FaceId f : in.faces())
    {
      const auto v = g.vertices_of(f);
      grid.add_triangle(g.vertex(v[0]).p, g.vertex(v[1]).p, g.vertex(v[2]).p);
    }
    print_range("clusters", grid.update());
  }
  g.remove_surface(in);
  return d;
}

//==============================================================================
// Distance, curvature, inside, merging and centres
//==============================================================================

std::string key(const char* what, double x, double y, double z)
{
  return std::format("{}@{:.17g},{:.17g},{:.17g}", what, x, y, z);
}

// Lines by key; vertices at the same position share a key.
using Lines = std::map<std::string, std::vector<std::string>>;

// Prints the lines sorted, numbering the values of a shared key in their
// sorted order, which does not depend on the order of the vertices.
void print_lines(Lines& lines)
{
  for (auto& [k, vs] : lines)
  {
    std::sort(vs.begin(), vs.end());
    for (std::size_t i = 0; i < vs.size(); ++i)
      std::println("{}{} {}", k, i ? std::format("#{}", i) : std::string(), vs[i]);
  }
}

void curvature_cmd(const Surface& s)
{
  const Geometrium& g = s.geometrium();
  Lines lines;
  for (VertexId v : s.vertices())
  {
    const Vec3 p = g.vertex(v).p;
    const std::string k = key("curv", p.x, p.y, p.z);
    const auto kh = mean_curvature_normal(s, v);
    const auto kg = gaussian_curvature(s, v);
    if ( ! kh || ! kg)
    {
      lines[k].push_back(std::format("mean={} gauss={}", int(kh.has_value()), int(kg.has_value())));
      continue;
    }
    const auto [k1, k2] = principal_curvatures(0.5*norm(*kh), *kg);
    const auto [e1, e2] = principal_directions(s, v, *kh, *kg);
    lines[k].push_back(std::format(
      "khx={:.12e} khy={:.12e} khz={:.12e} kg={:.12e} k1={:.12e} k2={:.12e} "
      "e1x={:.12e} e1y={:.12e} e1z={:.12e} e2x={:.12e} e2y={:.12e} e2z={:.12e}",
      kh->x, kh->y, kh->z, *kg, k1, k2, e1.x, e1.y, e1.z, e2.x, e2.y, e2.z));
  }
  print_lines(lines);
}

void inside_cmd(const Surface& s, unsigned n)
{
  const Geometrium& g = s.geometrium();
  BBox bb;
  for (VertexId v : s.vertices()) bb.add(g.vertex(v).p);
  const double b[6] = { bb.min.x, bb.min.y, bb.min.z, bb.max.x, bb.max.y, bb.max.z };
  const BBTree tree(s);
  unsigned n_in = 0;
  for (unsigned i = 0; i < n; ++i)
    for (unsigned j = 0; j < n; ++j)
      for (unsigned k = 0; k < n; ++k)
      {
        double q[3];
        const unsigned ijk[3] = { i, j, k };
        for (int c = 0; c < 3; ++c)
        {
          const double w = b[c + 3] - b[c];
          q[c] = b[c] - 0.1*w + 1.2*w*(ijk[c] + 0.5)/n;
        }
        const Vec3 p{ q[0], q[1], q[2] };
        const bool in = is_inside(tree, p);
        n_in += in;
        const double d = tree.point_distance(p);
        double d2;
        tree.point_closest(p, {}, &d2);
        // Closest points can tie, so only their distance is printed.
        std::println("{} in={} d={:.12e} d2={:.12e}", key("pt", q[0], q[1], q[2]), int(in), d, d2);
      }
  std::println("inside n={} in={}", n*n*n, n_in);
}

// Every face of in on vertices of its own, then merged back together.
Surface& merge_cmd(Geometrium& g, const Surface& in, double eps)
{
  Surface& s = g.new_surface();
  std::vector<VertexId> vertices;
  std::vector<EdgeId> edges;
  const std::vector<FaceId> faces(in.faces().begin(), in.faces().end());
  for (FaceId f : faces)
  {
    const auto o = g.vertices_of(f);
    VertexId v[3];
    for (int k = 0; k < 3; ++k) vertices.push_back(v[k] = g.add_vertex(g.vertex(o[k]).p));
    EdgeId e[3];
    for (int k = 0; k < 3; ++k) edges.push_back(e[k] = g.add_edge(v[k], v[(k + 1) % 3]));
    s.add_face(g.add_face(e[0], e[1], e[2]));
  }
  std::println("unwelded vertices={} edges={}", vertices.size(), edges.size());
  vertices = g.merge_vertices(vertices, eps);
  edges = g.merge_edges(edges);
  std::println("merged vertices={} edges={}", vertices.size(), edges.size());
  std::println("duplicate_faces_removed {}", remove_duplicate_faces(s));
  return s;
}

void centers_cmd(const Surface& s)
{
  const CenterOfMass cm = center_of_mass(s);
  const CenterOfArea ca = center_of_area(s);
  std::println("center_of_mass volume={:.12e} x={:.12e} y={:.12e} z={:.12e}", cm.volume, cm.center.x, cm.center.y, cm.center.z);
  std::println("center_of_area area={:.12e} x={:.12e} y={:.12e} z={:.12e}", ca.area, ca.center.x, ca.center.y, ca.center.z);
  std::println("boundary edges={}", boundary(s).size());
}

//==============================================================================
// Delaunay triangulation, as gts-ref
//==============================================================================

// delaunay IN [noconst] [nohull] [holes] [hullconst] [conform] [steiner=N] [quality=Q] [area=A]
Surface& delaunay_cmd(Geometrium& g, const char* path, std::span<char*> a)
{
  bool constraints = true, keep_hull = true, holes = false, conform = false, refine = false, hull_constraints = false;
  int steiner_max = -1;
  double quality = 0., max_area = std::numeric_limits<double>::max();
  for (std::string_view x : a)
  {
    if      (x == "noconst") constraints = false;
    else if (x == "nohull")  keep_hull = false;
    else if (x == "holes")   holes = true;
    else if (x == "conform") conform = true;
    else if (x == "hullconst") hull_constraints = true;
    else if (x.starts_with("steiner=")) steiner_max = std::atoi(x.data() + 8);
    else if (x.starts_with("quality=")) { conform = refine = true; quality = std::atof(x.data() + 8); }
    else if (x.starts_with("area="))    { conform = refine = true; max_area = std::atof(x.data() + 5); }
    else fail(std::string("unknown delaunay option ") + std::string(x));
  }

  std::ifstream in(path);
  if ( ! in) fail(std::string("cannot open ") + path);
  unsigned nv, ne, nf;
  if ( ! (in >> nv >> ne >> nf)) fail("bad header");
  std::vector<VertexId> vertices(nv);
  std::vector<Vec3> points(nv);
  for (unsigned i = 0; i < nv; ++i)
  {
    if ( ! (in >> points[i].x >> points[i].y >> points[i].z)) fail("bad vertex");
    vertices[i] = g.add_vertex(points[i]);
  }
  std::vector<EdgeId> edges(ne);
  for (auto& e : edges)
  {
    unsigned i, j;
    if ( ! (in >> i >> j) || i < 1 || i > nv || j < 1 || j > nv) fail("bad edge");
    e = g.add_edge(vertices[i - 1], vertices[j - 1]);
    if (constraints) set_constraint(g, e);
  }

  // gts_triangle_enclosing takes the points from a list built by prepending.
  std::vector<Vec3> reversed(points.rbegin(), points.rend());
  const FaceId t = enclosing_triangle(g, reversed, 100.);
  const auto hull = g.vertices_of(t);
  Surface& s = g.new_surface();
  s.add_face(t);

  unsigned n_duplicate = 0;
  for (VertexId v : vertices)
    if (const VertexId d = delaunay_add_vertex(s, v); d.valid())
    {
      ++n_duplicate;
      g.replace_vertex(v, d);
    }
  std::println("delaunay_duplicates {}", n_duplicate);

  if (constraints)
  {
    unsigned n_conflicts = 0;
    for (EdgeId e : edges) n_conflicts += unsigned(delaunay_add_constraint(s, e).size());
    std::println("delaunay_conflicts {}", n_conflicts);
  }

  {
    KeepOrphans keep(g, { .vertices = true });
    for (VertexId v : hull) g.destroy(v);
  }

  if ( ! keep_hull) delaunay_remove_hull(s);
  if (hull_constraints)
    for (EdgeId e : boundary(s)) set_constraint(g, e);
  if (holes)
    std::println("delaunay_holes {}", remove_faces_if(s, [&](FaceId f) {
      const auto ve = g.vertices_edges(f);
      for (int k = 0; k < 3; ++k)
        if (is_constraint(g, ve.e[k]) && g.edge(ve.e[k]).v[0] != ve.v[k]) return true;
      return false;
    }));

  if (conform)
  {
    const unsigned encroached = delaunay_conform(s, steiner_max);
    std::println("delaunay_encroached {}", encroached);
    if (encroached == 0 && refine)
    {
      const FaceCost cost = [&](const Geometrium& g, FaceId f) {
        const double q = g.quality(f), ar = g.area(f);
        return q < quality || ar > max_area ? q : 0.;
      };
      std::println("delaunay_unrefined {}", delaunay_refine(s, steiner_max, cost));
    }
  }
  std::println("delaunay_check {}", delaunay_check(s).valid() ? 1 : 0);
  return s;
}

//==============================================================================
// Height fields, as gts-ref
//==============================================================================

// Heights from a PGM file, P2 or P5, without comments; row j at y = j.
std::vector<double> read_pgm(const char* path, std::size_t& w, std::size_t& h)
{
  std::ifstream in(path, std::ios::binary);
  std::string magic;
  unsigned maxval;
  if ( ! (in >> magic >> w >> h >> maxval) || (magic != "P2" && magic != "P5")) fail("bad PGM header");
  in.get();
  std::vector<double> z(w * h);
  for (double& v : z)
  {
    unsigned x = 0;
    if (magic == "P2") { if ( ! (in >> x)) fail("bad PGM data"); }
    else if (maxval < 256) x = unsigned(in.get());
    else { x = unsigned(in.get()) << 8; x |= unsigned(in.get()); }
    v = x;
  }
  return z;
}

// happrox IN number=N|cost=C [relative=Z] [keep|closed]
void happrox_cmd(Surface& s, const char* path, std::span<char*> a)
{
  std::size_t number = 0;
  double cmin = 0., relative = 0.;
  bool by_cost = false;
  Enclosing enclosing = Enclosing::Remove;
  for (std::string_view x : a)
  {
    if      (x.starts_with("number="))   number = std::size_t(std::atoi(x.data() + 7));
    else if (x.starts_with("cost="))     { by_cost = true; cmin = std::atof(x.data() + 5); }
    else if (x.starts_with("relative=")) relative = std::atof(x.data() + 9);
    else if (x == "keep")   enclosing = Enclosing::Keep;
    else if (x == "closed") enclosing = Enclosing::Close;
    else fail(std::string("unknown happrox option ") + std::string(x));
  }
  const HeightCost cost = relative > 0. ? relative_height_cost(relative) : HeightCost(height_cost);
  const HeightStop stop = by_cost ? stop_below_cost(cmin) : stop_at_vertices(number);
  if (std::string_view(path).ends_with(".xyz"))
  {
    std::ifstream in(path);
    std::vector<Vec3> points;
    double x, y, z;
    while (in >> x >> y >> z) points.push_back({ x, y, z });
    // happrox -f prepends each point to its list.
    std::ranges::reverse(points);
    approximate_points(s, points, stop, cost, enclosing);
    return;
  }
  std::size_t w, h;
  const std::vector<double> z = read_pgm(path, w, h);
  approximate_grid(s, w, h, z, stop, cost);
}

// As gts-ref's strip_chain.
bool strip_chain(const std::vector<std::array<VertexId, 3>>& f)
{
  if (f.size() < 2) return true;
  auto in = [](VertexId v, const std::array<VertexId, 3>& t) { return v == t[0] || v == t[1] || v == t[2]; };
  std::vector<VertexId> shared, sq;
  for (VertexId v : f[0]) if (in(v, f[1])) shared.push_back(v); else sq.push_back(v);
  if (shared.size() != 2 || sq.size() != 1) return false;
  if (f.size() > 2 && in(shared[0], f[2])) std::swap(shared[0], shared[1]);
  sq.push_back(shared[0]);
  sq.push_back(shared[1]);
  for (std::size_t k = 1; k < f.size(); ++k)
  {
    if ( ! in(sq[k], f[k]) || ! in(sq[k + 1], f[k])) return false;
    VertexId next;
    for (VertexId v : f[k]) if (v != sq[k] && v != sq[k + 1]) next = v;
    sq.push_back(next);
  }
  return true;
}

// strip IN, as gts-ref.
void strip_cmd(const Surface& s)
{
  const Geometrium& g = s.geometrium();
  const auto list = strips(s);
  std::vector<bool> seen(g.face_capacity());
  Range len;
  bool cover = true, adjacent = true, chain = true;
  std::size_t total = 0;
  for (const auto& strip : list)
  {
    len.add(double(strip.size()));
    std::vector<std::array<VertexId, 3>> faces;
    for (FaceId f : strip) faces.push_back(g.vertices_of(f));
    if ( ! strip_chain(faces)) chain = false;
    for (std::size_t j = 0; j < strip.size(); ++j, ++total)
    {
      if (seen[strip[j].i]) cover = false;
      seen[strip[j].i] = true;
      if (j + 1 < strip.size() && ! g.common_edge(strip[j], strip[j + 1]).valid()) adjacent = false;
    }
  }
  len.update();
  if (total != s.face_number()) cover = false;
  std::println("strips {}", list.size());
  std::println("strip_faces {}", total);
  std::println("strip_cover {}", cover ? 1 : 0);
  std::println("strip_adjacent {}", adjacent ? 1 : 0);
  std::println("strip_chain {}", chain ? 1 : 0);
  print_range("strip_length", len);
}

//==============================================================================
// Progressive and hierarchical surfaces, as gts-ref
//==============================================================================

void print_level(const char* what, std::size_t level, const Surface& s)
{
  std::println("{}-{} vertices={} edges={} faces={} area={:.12e} volume={:.12e} closed={} manifold={} orientable={}",
               what, level, s.vertex_number(), s.edge_number(), s.face_number(), area(s), volume(s),
               is_closed(s) ? 1 : 0, is_manifold(s) ? 1 : 0, is_orientable(s) ? 1 : 0);
}

std::unique_ptr<PSurface> make_psurface(Surface& s, std::size_t min_edges)
{
  return std::make_unique<PSurface>(s, stop_below_edges(min_edges), EdgeCost{}, CollapsePosition{},
                                    std::numbers::pi / 180);
}

void psurface_levels(PSurface& ps, std::span<char*> levels)
{
  std::println("psurface min={} max={}", ps.min_vertex_number(), ps.max_vertex_number());
  for (char* l : levels)
  {
    const std::size_t n = std::strtoul(l, nullptr, 10);
    ps.set_vertex_number(n);
    print_level("level", n, ps.surface());
  }
  ps.set_vertex_number(ps.max_vertex_number());
}

void hsurface_cmd(Surface& s, std::size_t min_edges, unsigned n)
{
  const Geometrium& g = s.geometrium();
  auto ps = make_psurface(s, min_edges);
  auto x = [&](const VertexSplit& vs) { return g.vertex(vs.v).p.x; };
  HSurface hs(*ps, x, x);
  std::println("hsurface height={} vertices={}", hs.height(), hs.vertex_number());
  // GTS took the top of each heap. gts-cxx refuses a collapse out of order
  // that would break the surface, so the split of the smallest key that it
  // allows is taken.
  const auto& sp = hs.splits();
  auto best = [&](bool expand) {
    int b = -1;
    for (int i = 0; i < int(sp.size()); ++i)
      if ((expand ? hs.can_expand(i) : hs.can_collapse(i)) && (b < 0 || x(sp[i]) < x(sp[b]))) b = i;
    return b;
  };
  for (int i; (i = best(false)) >= 0;) hs.collapse(i);
  print_level("coarsest", hs.vertex_number(), s);
  for (unsigned k = 0; k < n; ++k)
  {
    const int i = best(true);
    if (i < 0) break;
    hs.expand(i);
  }
  print_level("expanded", hs.vertex_number(), s);
  for (unsigned k = 0; k < n / 2; ++k)
  {
    const int i = best(false);
    if (i < 0) break;
    hs.collapse(i);
  }
  print_level("collapsed", hs.vertex_number(), s);
  while (auto top = hs.top_expandable()) hs.expand(top->first);
}

//==============================================================================
// Graphs and partitions, as gts-ref
//==============================================================================

void graph_stats(const NodeSet& g, std::string_view key = "graph")
{
  Range degree;
  for (GNodeId n : g.nodes()) degree.add(gts::degree(g.graph(), n, &g));
  degree.update();
  std::println("{} nodes={} edges={} weight={:.6e}", key, g.size(), edges_of(g).size(), g.weight());
  print_range(std::format("{}_degree", key).c_str(), degree);
}

// The lines of text, each as "key_N line".
void print_lines(std::string_view key, const std::string& text)
{
  std::istringstream in(text);
  int i = 0;
  for (std::string line; std::getline(in, line); ) std::println("{}_{} {}", key, ++i, line);
}

// The text of print_stats(), and the nodes, edges and degrees of the dot
// file, whose labels and order differ from GTS's.
void graph_print(const NodeSet& g, std::string_view key)
{
  std::ostringstream out;
  print_stats(g, out);
  print_lines(key, out.str());
  std::ostringstream dot;
  write_dot(g, dot);
  print_dot(key, dot.str());
}

void partition_stats(const NodeSet& g, std::span<const NodeSet> parts)
{
  std::size_t n = 0;
  Range size;
  for (const NodeSet& p : parts) { n += p.size(); size.add(double(p.size())); }
  size.update();
  bool ok = true;
  for (GNodeId x : g.nodes())
  {
    unsigned k = 0;
    for (const NodeSet& p : parts) k += p.contains(x);
    if (k != 1) ok = false;
  }
  std::println("parts {}", parts.size());
  std::println("parts_nodes {}", n);
  std::println("parts_cover {}", ok ? 1 : 0);
  std::println("parts_cut {}", partition_edges_cut(parts));
  std::println("parts_balance {:.6e}", partition_balance(parts));
  print_range("parts_size", size);
  std::println("parts_cut_weight {:.6e}", partition_edges_cut_weight(parts));
  Range cut, cut_weight;
  for (const NodeSet& p : parts) { cut.add(edges_cut(p)); cut_weight.add(edges_cut_weight(p)); }
  cut.update();
  cut_weight.update();
  print_range("parts_edges_cut", cut);
  print_range("parts_edges_cut_weight", cut_weight);
  std::ostringstream out;
  print_partition_stats(parts, out);
  print_lines("parts_print", out.str());
}

void graph_cmd(const Surface& s, std::span<char*> a)
{
  Graph gr;
  NodeSet g = surface_graph(gr, s);
  graph_stats(g);
  const std::string_view what = a[0];
  if (what == "bubble" && a.size() == 3)
  {
    const auto parts = bubble_partition(g, unsigned(std::strtoul(a[1], nullptr, 10)), unsigned(std::strtoul(a[2], nullptr, 10)));
    partition_stats(g, parts);
  }
  else if (what == "bisect" && (a.size() == 3 || a.size() == 5))
  {
    const unsigned ntry = unsigned(std::strtoul(a[2], nullptr, 10));
    GraphBisection bg = std::string_view(a[1]) == "ggg" ? ggg_bisection(g, ntry) : bfgg_bisection(g, ntry);
    if (a.size() == 5)
    {
      const unsigned mmax = unsigned(std::strtoul(a[4], nullptr, 10));
      if (std::string_view(a[3]) == "kl") kl_refine(bg, mmax);
      else while (bkl_refine(bg, mmax, 0.1f) != 0.) {}
    }
    std::println("bisection_check {}", check(bg) ? 1 : 0);
    const NodeSet parts[2] = { bg.g1, bg.g2 };
    partition_stats(g, parts);
  }
  else if (what == "recursive" && a.size() == 2)
  {
    NodeSet whole = g;
    const auto parts = recursive_bisection(gr, whole, unsigned(std::strtoul(a[1], nullptr, 10)), 10, 50, 100, 0.1f);
    partition_stats(g, parts);
  }
  else if (what == "pgraph" && a.size() == 2)
  {
    PGraph pg(gr, g, std::strtoul(a[1], nullptr, 10));
    std::println("pgraph min={} max={} levels={}", pg.min_node_number(), pg.max_node_number(), pg.levels());
    graph_stats(g, "pgraph_coarse");
    pg.set_node_number(pg.max_node_number());
    graph_stats(g, "pgraph_fine");
    // Half the splits collapsed one by one, then expanded again.
    const std::size_t half = (pg.min_node_number() + pg.max_node_number()) / 2;
    unsigned n = 0;
    while (pg.node_number() > half && pg.remove_node()) ++n;
    std::println("pgraph_removed {} number={}", n, pg.node_number());
    graph_stats(g, "pgraph_half");
    n = 0;
    while (pg.add_node()) ++n;
    std::println("pgraph_added {} number={}", n, pg.node_number());
    graph_stats(g, "pgraph_back");
  }
  else if (what == "segments")
  {
    Graph sgr;
    std::vector<EdgeId> edges;
    for (EdgeId e : s.edges()) edges.push_back(e);
    NodeSet sg = segments_graph(sgr, s.geometrium(), edges);
    graph_stats(sg, "segments");
    graph_print(sg, "segments_print");
  }
  else if (what == "io")
  {
    std::stringstream io;
    write_graph(g, io);
    Graph gr1;
    auto g1 = read_graph(gr1, io);
    if ( ! g1) fail(std::format("cannot read the graph back: {}:{}: {}", g1.error().line, g1.error().pos, g1.error().message));
    graph_stats(*g1);
  }
  else if (what == "stats")
    graph_print(g, "graph_print");
  else
    fail("graph: WHAT must be stats, bubble, bisect, recursive, pgraph, segments or io");
}

} // namespace

int main(int argc, char* argv[])
{
  std::vector<char*> args(argv + 1, argv + argc);
  const char* out = nullptr;
  if (args.size() >= 2 && std::string_view(args[0]) == "-o")
  {
    out = args[1];
    args.erase(args.begin(), args.begin() + 2);
  }
  if (args.empty())
  {
    std::println(stderr, "usage: gts-cxx [-o OUT.gts] COMMAND ARGS...   (as gts-ref)");
    return 2;
  }

  const std::string_view cmd = args[0];
  std::span<char*> a(args.begin() + 1, args.end());

  std::print("command");
  for (char* x : args) std::print(" {}", x);
  std::println();

  Geometrium g;
  Surface* s = nullptr;

  if (cmd == "stats" && a.size() == 1)
  {
    s = &read_surface(g, a[0]);
  }
  else if (cmd == "coarsen" && a.size() >= 4)
  {
    s = &read_surface(g, a[0]);
    coarsen_cmd(*s, a.subspan(1));
  }
  else if (cmd == "refine" && a.size() == 4)
  {
    s = &read_surface(g, a[0]);
    refine_cmd(*s, a.subspan(1));
  }
  else if (cmd == "tessellate" && (a.size() == 2 || a.size() == 3))
  {
    s = &read_surface(g, a[0]);
    const bool mid = a.size() == 3 && std::string_view(a[2]) == "mid";
    // As GTSurf::Tessellate's mid-edge splitter: a + 0.5 (b - a).
    const Splitter mid_edge = [](Geometrium& g, EdgeId e) { return g.interpolate_vertex(g.edge(e).v[0], g.edge(e).v[1], 0.5); };
    for (unsigned long n = std::strtoul(a[1], nullptr, 10); n > 0; --n)
      tessellate(*s, mid ? Splitter(mid_edge) : Splitter(unit_sphere_arc_midvertex)); // GTS's default
  }
  else if (cmd == "bool" && a.size() == 3)
  {
    Surface& sa = read_surface(g, a[1]);
    Surface& sb = read_surface(g, a[2]);
    try { s = &make_boolean(g, a[0], sa, sb); }
    catch (const std::exception& e) { fail(e.what()); }
  }
  else if (cmd == "cluster" && a.size() == 2)
  {
    Surface& in = read_surface(g, a[0]);
    s = &make_cluster(g, in, std::strtod(a[1], nullptr));
  }
  else if (cmd == "iso" && a.size() == 3)
  {
    s = &g.new_surface();
    iso_cmd(*s, a);
  }
  else if (cmd == "sphere" && a.size() == 1)
  {
    s = &g.new_surface();
    // GTS refuses order 0 with a CRITICAL message and leaves the surface
    // empty; the summary is that of the empty surface.
    try { generate_sphere(*s, unsigned(std::strtoul(a[0], nullptr, 10))); }
    catch (const std::invalid_argument&) {}
  }
  else if (cmd == "distance" && (a.size() == 2 || a.size() == 3))
  {
    Surface& sa = read_surface(g, a[0]);
    Surface& sb = read_surface(g, a[1]);
    const double delta = a.size() == 3 ? std::strtod(a[2], nullptr) : 0.1;
    SurfaceDistance d = distance(sa, sb, delta);
    print_range("distance_a_to_b", d.faces);
    print_range("boundary_distance_a_to_b", d.boundary);
    d = distance(sb, sa, delta);
    print_range("distance_b_to_a", d.faces);
    print_range("boundary_distance_b_to_a", d.boundary);
  }
  else if (cmd == "curvature" && a.size() == 1)
  {
    curvature_cmd(read_surface(g, a[0]));
  }
  else if (cmd == "inside" && a.size() == 2)
  {
    inside_cmd(read_surface(g, a[0]), unsigned(std::strtoul(a[1], nullptr, 10)));
  }
  else if (cmd == "merge" && a.size() == 2)
  {
    Surface& in = read_surface(g, a[0]);
    s = &merge_cmd(g, in, std::strtod(a[1], nullptr));
    g.remove_surface(in);
  }
  else if (cmd == "centers" && a.size() == 1)
  {
    centers_cmd(read_surface(g, a[0]));
  }
  else if (cmd == "psurface" && a.size() >= 2)
  {
    s = &read_surface(g, a[0]);
    auto ps = make_psurface(*s, std::strtoul(a[1], nullptr, 10));
    psurface_levels(*ps, a.subspan(2));
  }
  else if (cmd == "psurface-write" && a.size() == 3)
  {
    s = &read_surface(g, a[0]);
    auto ps = make_psurface(*s, std::strtoul(a[1], nullptr, 10));
    std::ofstream o(a[2]);
    if ( ! o) fail(std::string("cannot write ") + a[2]);
    ps->write(o);
  }
  else if (cmd == "psurface-read" && a.size() >= 1)
  {
    std::ifstream in(a[0]);
    if ( ! in) fail(std::string("cannot open ") + a[0]);
    try
    {
      auto ps = PSurface::read(g, in);
      if ( ! ps) fail(std::format("cannot read {}: {}", a[0], ps.error().message));
      s = &(*ps)->surface();
      psurface_levels(**ps, a.subspan(1));
    }
    catch (const std::exception& e) { fail(e.what()); }
  }
  else if (cmd == "hsurface" && a.size() == 3)
  {
    s = &read_surface(g, a[0]);
    try { hsurface_cmd(*s, std::strtoul(a[1], nullptr, 10), unsigned(std::strtoul(a[2], nullptr, 10))); }
    catch (const std::exception& e) { fail(e.what()); }
  }
  else if (cmd == "jostle" && a.size() == 2)
  {
    std::ifstream in(a[0]);
    if ( ! in) fail(std::string("cannot open ") + a[0]);
    Graph gr;
    auto gs = read_jostle(gr, in);
    if ( ! gs) fail(std::format("cannot read {}: {}:{}: {}", a[0], gs.error().line, gs.error().pos, gs.error().message));
    graph_stats(*gs);
    NodeSet whole = *gs;
    const auto parts = recursive_bisection(gr, whole, unsigned(std::strtoul(a[1], nullptr, 10)), 10, 50, 100, 0.1f);
    partition_stats(*gs, parts);
  }
  else if (cmd == "graph" && a.size() >= 2)
  {
    try { graph_cmd(read_surface(g, a[0]), a.subspan(1)); }
    catch (const std::exception& e) { fail(e.what()); }
  }
  else if (cmd == "strip" && a.size() == 1)
  {
    s = &read_surface(g, a[0]);
    strip_cmd(*s);
  }
  else if (cmd == "isoslice" && a.size() == 2)
  {
    s = &g.new_surface();
    isoslice_cmd(*s, a);
  }
  else if (cmd == "delaunay" && a.size() >= 1)
  {
    try { s = &delaunay_cmd(g, a[0], a.subspan(1)); }
    catch (const std::exception& e) { fail(e.what()); }
  }
  else if (cmd == "happrox" && a.size() >= 2)
  {
    s = &g.new_surface();
    try { happrox_cmd(*s, a[0], a.subspan(1)); }
    catch (const std::exception& e) { fail(e.what()); }
  }
  else
    unsupported(cmd);

  if (s)
  {
    print_summary(*s);
    if (out && ! write_gts(*s, std::filesystem::path(out)))
      fail(std::string("cannot write ") + out);
  }

  std::println("glib_critical 0");
  std::println("glib_warning 0");
  std::println("status ok");
  return 0;
}
