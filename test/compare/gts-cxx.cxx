// gts-cxx: the twin of gts-ref, running the same operations through the C++
// library and printing the same summary. Operations not ported yet print
// "status unsupported" and exit with 3.
//
// Copyright (C) 2026 Matevz Tadel
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/all.h"

#include <map>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <numbers>
#include <stdexcept>
#include <print>
#include <span>
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
  std::println("self_intersecting {}", is_self_intersecting(s) ? 1 : 0);

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
  if      (func == "sphere") { f = field_sphere; iso = 1;    }
  else if (func == "torus")  { f = field_torus;  iso = 0.09; }
  else if (func == "blobs")  { f = field_blobs;  iso = 4;    }
  else fail("FUNC must be sphere, torus or blobs");

  const unsigned n = unsigned(std::strtoul(a[2], nullptr, 10));
  const double lo = -1.25, d = 2.5 / n;
  const CartesianGrid g = { n + 1, n + 1, n + 1, lo, d, lo, d, lo, d };
  const PlaneFunc pf = sample(f);
  if      (algo == "cartesian")     isosurface_cartesian    (s, g, pf, iso);
  else if (algo == "tetra")         isosurface_tetra        (s, g, pf, iso);
  else if (algo == "tetra_bounded") isosurface_tetra_bounded(s, g, pf, iso);
  else if (algo == "tetra_bcl")     isosurface_tetra_bcl    (s, g, pf, iso);
  else fail("ALGO must be cartesian, tetra, tetra_bounded or tetra_bcl");
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
