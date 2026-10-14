// Unit tests of the core: Geometrium, surfaces, rings, collection,
// attributes, renumbering callouts, statistics and .gts input and output;
// default hooks, split, collapse and swap primitives and their events; the edge
// heap, tessellation, refinement, spheres, coarsening, iso-surfaces,
// predicates, bounding-box trees, booleans and the cluster grid; distances,
// the inside test, curvature, the kd-tree and merging; neighbourhoods,
// Delaunay triangulation, matrices, box overlaps, the writers, strips,
// parked elements, progressive and hierarchical surfaces, graphs, the
// iso-surface slices and marching cubes from seeds; face events, points
// kept by faces, and height fields.
//
// Copyright (C) 2026 Matevz Tadel
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/all.h"

#include <atomic>
#include <algorithm>
#include <cmath>
#include <map>
#include <numbers>
#include <print>
#include <sstream>
#include <string>
#include <thread>

using namespace gts;

// A type with its own default hooks: interpolation takes the nearer source.
struct Label { int id = 0; };
template <> struct gts::DefaultHooks<Label> : gts::GenericHooks<Label>
{
  static Label interpolate(const Label& a, const Label& b, double t) { return t < 0.5 ? a : b; }
};

static std::atomic<int> g_failed = 0;
#define CHECK(c) do { if ( ! (c)) { std::println("FAIL {}:{}: {}", __FILE__, __LINE__, #c); ++g_failed; } } while (0)

template <class R> static std::size_t count(R&& r) { std::size_t n = 0; for (auto x : r) { (void) x; ++n; } return n; }

struct Tet { std::array<VertexId, 4> v; std::array<EdgeId, 6> e; std::array<FaceId, 4> f; };

// A tetrahedron with its faces oriented outwards. Edges: 01 02 03 12 13 23.
static Tet tetrahedron(Geometrium& g, Surface& s, Vec3 o = {})
{
  Tet t;
  t.v = { g.add_vertex(o), g.add_vertex(o + Vec3{1, 0, 0}), g.add_vertex(o + Vec3{0, 1, 0}), g.add_vertex(o + Vec3{0, 0, 1}) };
  t.e = { g.add_edge(t.v[0], t.v[1]), g.add_edge(t.v[0], t.v[2]), g.add_edge(t.v[0], t.v[3]),
          g.add_edge(t.v[1], t.v[2]), g.add_edge(t.v[1], t.v[3]), g.add_edge(t.v[2], t.v[3]) };
  // Bottom 0-2-1, front 0-1-3, left 0-3-2, slanted 1-2-3, all seen from outside.
  t.f = { g.add_face(t.e[1], t.e[3], t.e[0]), g.add_face(t.e[0], t.e[4], t.e[2]),
          g.add_face(t.e[2], t.e[5], t.e[1]), g.add_face(t.e[3], t.e[5], t.e[4]) };
  for (FaceId f : t.f) s.add_face(f);
  return t;
}

static void test_build_and_rings()
{
  Geometrium g;
  Surface& s = g.new_surface("tet");
  Tet t = tetrahedron(g, s);
  CHECK(s.vertex_number() == 4 && s.edge_number() == 6 && s.face_number() == 4);
  for (VertexId v : t.v) CHECK(count(g.edges_of(v)) == 3);
  for (EdgeId e : t.e)   CHECK(count(g.faces_of(e)) == 2);
  CHECK(is_closed(s) && is_manifold(s) && is_orientable(s));
  CHECK(std::abs(area(s) - (1.5 + std::sqrt(3.0) / 2)) < 1e-12);
  CHECK(std::abs(volume(s) - 1.0 / 6) < 1e-15);
  CHECK(g.find_edge(t.v[0], t.v[3]) == t.e[2] && ! g.find_edge(t.v[0], t.v[0]).valid());
}

static void test_orientation()
{
  // gts_triangle_vertices: the four ways the first two edges can share a vertex.
  Geometrium g;
  VertexId a = g.add_vertex({0, 0, 0}), b = g.add_vertex({1, 0, 0}), c = g.add_vertex({0, 1, 0});
  EdgeId ab = g.add_edge(a, b), ba = g.add_edge(b, a), bc = g.add_edge(b, c), cb = g.add_edge(c, b);
  EdgeId ca = g.add_edge(c, a);
  using V3 = std::array<VertexId, 3>;
  CHECK(g.vertices_of(g.add_face(ab, bc, ca)) == (V3{ a, b, c })); // e1.v2 == e2.v1
  CHECK(g.vertices_of(g.add_face(ab, cb, ca)) == (V3{ a, b, c })); // e1.v2 == e2.v2
  CHECK(g.vertices_of(g.add_face(ba, ab == ab ? g.add_edge(b, c) : bc, ca)) == (V3{ a, b, c })); // e1.v1 == e2.v1
  CHECK(g.vertices_of(g.add_face(ba, cb, ca)) == (V3{ a, b, c })); // e1.v1 == e2.v2... b is e1.v1 and e2.v2
  // Counter-clockwise seen from +z: the normal points to +z.
  const Vec3 n = g.normal(g.add_face(ab, bc, ca));
  CHECK(n.z > 0 && n.x == 0 && n.y == 0);
}

static void test_release_cascade()
{
  // A face leaving its last surface is unlinked at once, with its orphans.
  Geometrium g;
  Surface& s = g.new_surface();
  Tet t = tetrahedron(g, s);
  s.remove_face(t.f[0]);
  CHECK(g.dead(t.f[0]) && ! g.dead(t.e[0]) && count(g.faces_of(t.e[0])) == 1);
  s.remove_face(t.f[1]);                // edge 01 loses its last face
  CHECK(g.dead(t.e[0]) && count(g.edges_of(t.v[0])) == 2);
  CHECK( ! g.dead(t.v[0]) && ! g.dead(t.v[1]));
  s.remove_face(t.f[2]);
  s.remove_face(t.f[3]);                // everything goes
  CHECK(g.dead(t.v[0]) && g.dead(t.v[3]));
  g.collect();
  CHECK(g.vertex_capacity() == 0 && g.edge_capacity() == 0 && g.face_capacity() == 0);
}

static void test_keep_orphans()
{
  Geometrium g;
  Surface& s = g.new_surface();
  Tet t = tetrahedron(g, s);
  {
    KeepOrphans keep(g, { .edges = true });
    s.remove_face(t.f[0]);
    s.remove_face(t.f[1]);
  }
  CHECK(g.dead(t.f[0]) && ! g.dead(t.e[0]) && count(g.faces_of(t.e[0])) == 0);
  g.collect({ .edges = false });        // keep the faceless edge
  CHECK(g.face_capacity() == 2 && g.edge_capacity() == 6);
  g.collect();
  CHECK(g.edge_capacity() == 5);
}

static void test_shared_faces()
{
  Geometrium g;
  Surface& s1 = g.new_surface("one");
  Surface& s2 = g.new_surface("two");
  Tet t = tetrahedron(g, s1);
  s2.add_face(t.f[0]);
  s1.remove_face(t.f[0]);
  CHECK( ! g.dead(t.f[0]) && g.face(t.f[0]).surfaces == 1);
  g.collect();
  CHECK(g.face_capacity() == 4 && s1.face_number() == 3 && s2.face_number() == 1);
  g.remove_surface(s2);
  g.collect();
  CHECK(g.face_capacity() == 3 && g.surface_number() == 1);
}

static void test_non_manifold_and_duplicates()
{
  Geometrium g;
  Surface& s = g.new_surface();
  VertexId a = g.add_vertex({0, 0, 0}), b = g.add_vertex({1, 0, 0});
  EdgeId ab = g.add_edge(a, b);
  for (int k = 0; k < 3; ++k)
  {
    VertexId c = g.add_vertex({0.5, double(k), 1});
    s.add_face(g.add_face(ab, g.add_edge(b, c), g.add_edge(c, a)));
  }
  CHECK(count(g.faces_of(ab)) == 3 && count(g.edges_of(a)) == 4);
  CHECK( ! is_manifold(s) && ! is_closed(s));
  SurfaceStats st = stats(s);
  CHECK(st.n_non_manifold_edges == 1 && st.n_boundary_edges == 6 && st.n_duplicate_edges == 0);

  EdgeId ab2 = g.add_edge(b, a);
  CHECK(g.duplicate_of(ab2) == ab && g.duplicate_of(ab) == ab2);
  FaceId f0 = *s.faces().begin();
  FaceId dup = g.add_face(g.face(f0).e[0], g.face(f0).e[1], g.face(f0).e[2]);
  CHECK(g.duplicate_of(dup) == f0);
}

static void test_ring_order()
{
  Geometrium g;
  VertexId c = g.add_vertex({0, 0, 0});
  std::array<EdgeId, 4> e;
  for (int k = 0; k < 4; ++k) e[k] = g.add_edge(c, g.add_vertex({double(k), 1, 0}));
  auto ring = [&] { std::vector<Index> r; for (EdgeId x : g.edges_of(c)) r.push_back(x.i); return r; };
  CHECK(ring() == (std::vector<Index>{ e[3].i, e[2].i, e[1].i, e[0].i }));
  g.kill_edge(e[0]);
  g.kill_edge(e[2]);
  CHECK(ring() == (std::vector<Index>{ e[3].i, e[1].i }));
  EdgeId e4 = g.add_edge(c, g.add_vertex({9, 9, 9}));
  CHECK(ring() == (std::vector<Index>{ e4.i, e[3].i, e[1].i }));
  std::vector<EdgeId> kept = { e4, e[3], e[1] };
  Subscription sub = g.on_renumber([&](const Renumbering& r) { for (auto& x : kept) x = r(x); });
  g.collect({ .edges = false, .vertices = false });
  CHECK(ring() == (std::vector<Index>{ kept[0].i, kept[1].i, kept[2].i }));
}

static void test_attributes()
{
  Geometrium g;
  Surface& s = g.new_surface();
  auto pressure = g.add_vertex_attribute<float>("pressure",
    Hooks<float>().on_interpolate([](float a, float b, double t) { return float(a + t * (b - a)); }));
  auto stroke = g.add_vertex_attribute<int>("stroke");
  auto crease = g.add_edge_attribute<std::uint8_t>("crease");
  auto colour = g.add_face_attribute<int>("colour");
  Tet t = tetrahedron(g, s);
  for (int i = 0; i < 4; ++i) { pressure[t.v[i]] = float(i); stroke[t.v[i]] = 10 + i; colour[t.f[i]] = 100 + i; }
  crease[t.e[5]] = 1;
  VertexId m = g.interpolate_vertex(t.v[1], t.v[3], 0.25);
  CHECK(pressure[m] == 1.5f && stroke[m] == 11 && g.vertex(m).p == (Vec3{0.75, 0, 0.25}));

  s.remove_face(t.f[0]); s.remove_face(t.f[1]); s.remove_face(t.f[2]);
  g.collect(); // left: face 3 (1-2-3) with its three edges and vertices
  CHECK(g.vertex_capacity() == 3 && g.edge_capacity() == 3 && g.face_capacity() == 1);
  float sum = 0; int ssum = 0;
  for (Index i = 0; i < 3; ++i) { sum += pressure[VertexId(i)]; ssum += stroke[VertexId(i)]; }
  CHECK(sum == 1 + 2 + 3 && ssum == 11 + 12 + 13);
  CHECK(colour[FaceId(0)] == 103);
  int creases = 0;
  for (Index i = 0; i < 3; ++i) creases += crease[EdgeId(i)];
  CHECK(creases == 1);
}

static void test_import()
{
  Geometrium g1, g2;
  Surface& s = g1.new_surface("src");
  auto p1 = g1.add_vertex_attribute<float>("pressure");
  Tet t = tetrahedron(g1, s);
  for (int i = 0; i < 4; ++i) p1[t.v[i]] = float(i + 1);
  auto p2 = g2.add_vertex_attribute<float>("pressure");
  Surface& d = g2.import(g1, s);
  CHECK(d.face_number() == 4 && d.edge_number() == 6 && d.vertex_number() == 4 && d.name() == "src");
  float sum = 0;
  for (Index i = 0; i < g2.vertex_capacity(); ++i) sum += p2[VertexId(i)];
  CHECK(sum == 10 && std::abs(volume(d) - volume(s)) < 1e-15);
}

static void test_io()
{
  Geometrium g;
  Surface& s = g.new_surface();
  tetrahedron(g, s);
  std::stringstream ss;
  write_gts(s, ss, { .precision = 17 });

  Geometrium h;
  auto r = read_gts(h, ss);
  CHECK(r.has_value());
  if (r)
  {
    Surface& t = **r;
    CHECK(t.vertex_number() == 4 && t.edge_number() == 6 && t.face_number() == 4);
    CHECK(volume(t) == volume(s) && area(t) == area(s));
    std::stringstream again;
    write_gts(t, again, { .precision = 17 });
    std::stringstream first;
    write_gts(s, first, { .precision = 17 });
    CHECK(again.str() == first.str());
  }

  // A comment line, and extra tokens on a line, are skipped.
  std::istringstream in("# a comment\n3 3 1 GtsSurface GtsFace GtsEdge GtsVertex\n"
                        "0 0 0 extra\n1 0 0\n! another comment\n0 1 0\n1 2\n2 3\n3 1\n1 2 3\n");
  auto c = read_gts(h, in);
  CHECK(c.has_value() && (*c)->face_number() == 1);

  // A comment at the end of a line ends the line; GTS read the next line as
  // part of it.
  std::istringstream trailing("3 3 1\n0 0 0\n1 0 0\n0 1 0 ! trailing\n1 2 # edge\n2 3\n3 1\n1 2 3\n");
  auto tr = read_gts(h, trailing);
  CHECK(tr.has_value() && (*tr)->face_number() == 1);

  // An error leaves the Geometrium as it was.
  const std::size_t nv = h.vertex_capacity(), ns = h.surface_number();
  std::istringstream bad("3 3 1\n0 0 0\n1 0 0\n0 1 0\n1 2\n2 7\n3 1\n1 2 3\n");
  auto e = read_gts(h, bad);
  CHECK( ! e.has_value() && h.vertex_capacity() == nv && h.surface_number() == ns);
  if ( ! e) CHECK(e.error().line == 6 && e.error().message.find("out of range") != std::string::npos);

  std::istringstream notri("4 3 1\n0 0 0\n1 0 0\n0 1 0\n5 5 5\n1 2\n2 3\n1 4\n1 2 3\n");
  auto n = read_gts(h, notri);
  CHECK( ! n.has_value() && h.vertex_capacity() == nv);
}

static void test_threads()
{
  std::vector<std::thread> th;
  for (int k = 0; k < 4; ++k)
    th.emplace_back([k] {
      Geometrium g;
      Surface& s = g.new_surface();
      for (int i = 0; i < 200; ++i) tetrahedron(g, s, Vec3{double(i), double(k), 0});
      for (int i = 0; i < 100; ++i) s.remove_face(FaceId(Index(4 * i)));
      g.collect();
      CHECK(g.face_capacity() == 700 && s.face_number() == 700);
    });
  for (auto& x : th) x.join();
  th.clear();

  Geometrium g;
  Surface& s = g.new_surface();
  for (int i = 0; i < 200; ++i) tetrahedron(g, s, Vec3{double(i), 0, 0});
  for (int k = 0; k < 4; ++k)
    th.emplace_back([&] {
      std::size_t n = 0;
      for (Index i = 0; i < g.vertex_capacity(); ++i) n += count(g.edges_of(VertexId(i)));
      CHECK(n == 2 * g.edge_capacity() && s.vertex_number() == 800 && is_closed(s));
    });
  for (auto& x : th) x.join();
}

static void test_set_face_edges()
{
  // Three faces on edge ab (non-manifold), then the middle one moves to
  // another edge and its kept edge changes slot.
  Geometrium g;
  VertexId a = g.add_vertex({0, 0, 0}), b = g.add_vertex({1, 0, 0});
  VertexId c = g.add_vertex({0, 1, 0}), d = g.add_vertex({0, 0, 1}), e = g.add_vertex({0, -1, 0});
  VertexId m = g.add_vertex({0.5, 0, 0});
  EdgeId ab = g.add_edge(a, b);
  EdgeId bc = g.add_edge(b, c), ca = g.add_edge(c, a);
  EdgeId bd = g.add_edge(b, d), da = g.add_edge(d, a);
  EdgeId be = g.add_edge(b, e), ea = g.add_edge(e, a);
  FaceId f1 = g.add_face(ab, bc, ca), f2 = g.add_face(ab, bd, da), f3 = g.add_face(ab, be, ea);
  auto ring = [&](EdgeId x) { std::vector<FaceId> r; for (FaceId f : g.faces_of(x)) r.push_back(f); return r; };
  CHECK(ring(ab) == (std::vector<FaceId>{ f3, f2, f1 }));

  // f2 = (ab, bd, da) becomes (da, dm, ma): it leaves ab and bd, keeps da
  // (from slot 2 to slot 0), joins dm and ma.
  EdgeId am = g.add_edge(a, m), dm = g.add_edge(d, m);
  FaceId other = g.add_face(da, bd, g.add_edge(a, b)); // a second face on da and bd
  CHECK(ring(da) == (std::vector<FaceId>{ other, f2 }));
  g.set_face_edges(f2, da, dm, am);
  CHECK(ring(ab) == (std::vector<FaceId>{ f3, f1 }));
  CHECK(ring(bd) == (std::vector<FaceId>{ other }));
  CHECK(ring(da) == (std::vector<FaceId>{ other, f2 })); // kept its place
  CHECK(ring(dm) == (std::vector<FaceId>{ f2 }) && ring(am) == (std::vector<FaceId>{ f2 }));
  CHECK(g.next_around(f2, da) == other && g.face(f2).e[0] == da);
  // Edges that do not form a triangle are refused, and nothing changes.
  bool threw = false;
  try { g.set_face_edges(f1, ab, bd, ea); } catch (const std::invalid_argument&) { threw = true; }
  CHECK(threw && ring(ab) == (std::vector<FaceId>{ f3, f1 }) && g.face(f1).e[0] == ab);
}

static void test_eheap()
{
  // Smallest key first; equal keys come off in GtsEHeap's order.
  EHeap<int> h([](const int& x) { return double(x / 10); });
  for (int x : { 31, 10, 22, 11, 30, 21, 12 }) h.insert(x);
  std::vector<int> out;
  while (auto t = h.remove_top()) out.push_back(t->first);
  // The order of GtsEHeap for the same insertions, from a C program.
  CHECK(out == (std::vector<int>{ 10, 12, 11, 22, 21, 30, 31 }));

  EHeap<int> f([](const int& x) { return double(x); });
  f.freeze();
  for (int x : { 5, 3, 8, 1, 9, 2 }) f.insert(x);
  f.thaw();
  const auto h9 = f.insert(9, 9.0);
  f.decrease_key(h9, 0.0);
  CHECK(f.top()->first == 9 && f.top()->second == 0.0);
  const auto h7 = f.insert(7);
  CHECK(f.remove(h7) == 7 && f.size() == 7);
  out.clear();
  while (auto t = f.remove_top()) out.push_back(t->first);
  CHECK(out == (std::vector<int>{ 9, 1, 2, 3, 5, 8, 9 }));
}

static void test_default_hooks()
{
  Geometrium g;
  auto d = g.add_vertex_attribute<double>("d");
  auto i = g.add_vertex_attribute<int>("i");
  auto p = g.add_vertex_attribute<Vec3>("p");
  auto s = g.add_vertex_attribute<std::string>("s");
  auto l = g.add_vertex_attribute<Label>("l");
  auto h = g.add_vertex_attribute<double>("h", Hooks<double>().on_interpolate([](double, double, double) { return -1.0; }));
  const VertexId a = g.add_vertex({0, 0, 0}), b = g.add_vertex({1, 0, 0});
  d[a] = 1; d[b] = 3; i[a] = 1; i[b] = 3; p[a] = {0, 0, 0}; p[b] = {2, 4, 6};
  s[a] = "a"; s[b] = "b"; l[a] = {1}; l[b] = {3}; h[a] = 1; h[b] = 3;
  const VertexId m = g.interpolate_vertex(a, b, 0.75);
  CHECK(d[m] == 2.5);              // affine: linear
  CHECK(i[m] == 1);                // integers are labels: the first source
  CHECK(p[m] == (Vec3{1.5, 3, 4.5}));
  CHECK(s[m] == "a");
  CHECK(l[m].id == 3);             // the specialization
  CHECK(h[m] == -1);               // the hook given for this attribute
}

static void test_split_edge()
{
  // Split edge 01 of a tetrahedron that is in two surfaces. The faces on it
  // are cut, the pieces join both surfaces, the hooks set the attributes,
  // and a subscriber gives the cuts the colour of their face.
  Geometrium g;
  auto colour = g.add_face_attribute<int>("colour");
  auto w = g.add_edge_attribute<double>("w", Hooks<double>().on_split([](double x) { return x / 2; }));
  auto ecolour = g.add_edge_attribute<int>("colour");
  Surface& s = g.new_surface();
  Tet t = tetrahedron(g, s);
  Surface& s2 = g.new_surface();
  s2.add_face(t.f[1]);
  for (int k = 0; k < 4; ++k) colour[t.f[k]] = 10 + k;
  w[t.e[0]] = 8;

  std::vector<EdgeSplitEvent> edges;
  std::vector<std::pair<FaceId, std::vector<FaceId>>> faces;
  Subscription sub = g.on_split([&](const SplitEvent& ev) {
    edges.insert(edges.end(), ev.edges.begin(), ev.edges.end());
    for (const auto& f : ev.faces)
    {
      faces.emplace_back(f.parent, std::vector<FaceId>(f.pieces.begin(), f.pieces.end()));
      for (EdgeId c : f.cuts) ecolour[c] = colour[f.parent];
    }
  });

  const VertexId m = g.interpolate_vertex(t.v[0], t.v[1], 0.5);
  std::vector<EdgeId> made;
  const auto halves = g.split_edge(t.e[0], m, [&](EdgeId e) { made.push_back(e); });
  CHECK(made.size() == 4 && made[0] == halves[0] && made[1] == halves[1]);
  CHECK(g.edge(halves[0]).v[0] == t.v[0] && g.edge(halves[1]).v[0] == t.v[1]);
  CHECK(g.dead(t.e[0]) && w[halves[0]] == 4 && w[halves[1]] == 4);
  CHECK(s.face_number() == 6 && is_closed(s) && is_orientable(s) && std::abs(volume(s) - 1.0 / 6) < 1e-15);
  CHECK(s2.face_number() == 2); // f1 and its new piece
  CHECK(edges.size() == 1 && edges[0].parent == t.e[0] && edges[0].vertex == m && edges[0].halves == halves);
  CHECK(faces.size() == 2);
  for (const auto& [parent, pieces] : faces)
  {
    CHECK(pieces.size() == 2 && pieces[0] == parent);
    for (FaceId f : pieces) CHECK(colour[f] == colour[parent]);
  }
  CHECK(ecolour[made[2]] != 0 && ecolour[made[3]] != 0 && ecolour[made[2]] != ecolour[made[3]]);

  // Split a face in four: three edge splits (one edge is a half already
  // split above, but the memo is new), one face split with three cuts.
  edges.clear(); faces.clear();
  g.split_face4(t.f[3], midvertex);
  CHECK(edges.size() == 3 && faces.size() == 1 && faces[0].second.size() == 4);
  CHECK(s.face_number() == 9 && ! is_closed(s)); // the neighbours keep the old edges: open seams
  sub.reset();
  edges.clear();
  g.split_face4(*s.faces().begin(), midvertex);
  CHECK(edges.empty() && s.face_number() == 12);
}

// An octahedron with its faces oriented outwards: vertices +x +y -x -y +z -z.
static std::array<VertexId, 6> octahedron(Geometrium& g, Surface& s)
{
  std::array<VertexId, 6> v = { g.add_vertex({1, 0, 0}), g.add_vertex({0, 1, 0}), g.add_vertex({-1, 0, 0}),
                                g.add_vertex({0, -1, 0}), g.add_vertex({0, 0, 1}), g.add_vertex({0, 0, -1}) };
  auto edge = [&](VertexId a, VertexId b) { EdgeId e = g.find_edge(a, b); return e.valid() ? e : g.add_edge(a, b); };
  for (int k = 0; k < 4; ++k)
  {
    const VertexId a = v[k], b = v[(k + 1) % 4];
    s.add_face(g.add_face(edge(a, b), edge(b, v[4]), edge(v[4], a)));
    s.add_face(g.add_face(edge(b, a), edge(a, v[5]), edge(v[5], b)));
  }
  return v;
}

static void test_collapse()
{
  {
    // Every edge of a tetrahedron belongs to a tetrahedron.
    Geometrium g;
    Surface& s = g.new_surface();
    Tet t = tetrahedron(g, s);
    for (EdgeId e : t.e) CHECK( ! g.collapsible(e));
  }

  Geometrium g;
  auto w = g.add_vertex_attribute<double>("w");
  auto label = g.add_vertex_attribute<int>("label");
  auto n = g.add_edge_attribute<int>("n", Hooks<int>().on_merge([](int a, int b) { return a + b; }));
  Surface& s = g.new_surface();
  const auto v = octahedron(g, s);
  CHECK(is_closed(s) && is_orientable(s) && volume(s) > 0);
  for (int k = 0; k < 6; ++k) { w[v[k]] = k; label[v[k]] = 10 + k; }
  for (Index i = 0; i < g.edge_capacity(); ++i) n[EdgeId(i)] = 1;

  const EdgeId e = g.find_edge(v[4], v[0]);
  CHECK(g.collapsible(e));
  // Faces folding back onto each other; the C library agrees on this and on
  // a grid of 9261 positions around the octahedron.
  CHECK(g.collapse_folds(e, {-3, -3, 0}, std::pow(std::cos(std::numbers::pi / 180), 2)));
  CHECK( ! g.collapse_folds(e, {0.5, 0, 0.5}, std::pow(std::cos(std::numbers::pi / 180), 2)));

  std::vector<CollapseEvent> events;
  std::vector<std::array<EdgeId, 2>> merged;
  Subscription sub = g.on_collapse([&](const CollapseEvent& ev) {
    events.push_back(ev);
    merged.assign(ev.merged_edges.begin(), ev.merged_edges.end());
  });
  std::vector<EdgeId> removed;
  const VertexId m = g.collapse_edge(e, {0.5, 0, 0.5}, [&](EdgeId x) { removed.push_back(x); });
  CHECK(m.valid() && g.dead(v[4]) && g.dead(v[0]) && g.dead(e));
  CHECK(s.vertex_number() == 5 && s.edge_number() == 9 && s.face_number() == 6);
  CHECK(is_closed(s) && is_manifold(s) && is_orientable(s) && volume(s) > 0);
  CHECK(w[m] == 2 && label[m] == 14); // the mean of 4 and 0; the first source
  CHECK(events.size() == 1 && events[0].faces.size() == 2 && events[0].vertex == m);
  CHECK(merged.size() == 2 && removed.size() == 2 && events[0].removed_edges.empty());
  for (auto [kept, gone] : merged)
  {
    CHECK( ! g.dead(kept) && g.dead(gone) && n[kept] == 2);
    CHECK(std::ranges::find(removed, gone) != removed.end());
  }
  for (EdgeId x : g.edges_of(m)) CHECK( ! g.dead(x) && count(g.faces_of(x)) == 2);
}

static void test_coarsen()
{
  Geometrium g;
  Surface& s = generate_sphere(g.new_surface(), 4);
  const double v0 = volume(s);
  coarsen(s, stop_below_edges(300), [](const Geometrium& g, EdgeId e) { return volume_optimized_cost(g, e); },
          [](const Geometrium& g, EdgeId e) { return volume_optimized_vertex(g, e); }, std::numbers::pi / 180);
  // Stops before a collapse would leave fewer than 300 edges.
  CHECK(s.edge_number() >= 300 && s.edge_number() < 303);
  CHECK(is_closed(s) && is_manifold(s) && is_orientable(s));
  CHECK(std::abs(volume(s) - v0) < 1e-3 * v0); // volume-optimized collapses keep the volume
  for (VertexId x : s.vertices()) CHECK( ! g.dead(x));
  g.collapse_folds(*s.edges().begin(), {}, 1.0); // queries work on the coarsened surface
  g.collect();
  CHECK(g.edge_capacity() == s.edge_number() && is_closed(s));
}

static void test_iso()
{
  // sample() evaluates plane k at grid_point(): g.x + i * g.dx and so on,
  // which adding up the steps would not give with a step of 0.1.
  {
    const CartesianGrid h = { 21, 21, 21, -1, 0.1, -1, 0.1, -1, 0.1 };
    PlaneValues plane(h.nx, h.ny);
    sample([](Vec3 p) { return p.x + 10 * p.y + 100 * p.z; })(plane, h, 13);
    const double x = h.x + 20 * h.dx, y = h.y + 7 * h.dy, z = h.z + 13 * h.dz;
    CHECK(plane(20, 7) == x + 10 * y + 100 * z);
    double xs = h.x;
    for (unsigned i = 0; i < 20; ++i) xs += h.dx;
    CHECK(xs != x); // the steps add up to another value
  }
  const CartesianGrid g = { 21, 21, 21, -1.25, 0.125, -1.25, 0.125, -1.25, 0.125 };

  // A potential, larger inside: the faces point outwards, towards f < iso.
  auto sphere = sample([](Vec3 p) { return - (p.x * p.x + p.y * p.y + p.z * p.z); });
  {
    Geometrium geo;
    Surface& s = geo.new_surface();
    isosurface_cartesian(s, g, sphere, -1.0);
    CHECK(is_closed(s) && is_manifold(s) && is_orientable(s));
    CHECK(volume(s) > 0 && std::abs(volume(s) - 4 * std::numbers::pi / 3) < 0.05);
    // Vertices are made only where faces use them.
    CHECK(geo.vertex_capacity() == s.vertex_number());
    for (VertexId v : s.vertices()) CHECK(std::abs(norm(geo.vertex(v).p) - 1) < 0.02);
  }
  {
    Geometrium geo;
    Surface& s = geo.new_surface();
    isosurface_tetra(s, g, sphere, -1.0);
    CHECK(is_closed(s) && is_manifold(s) && volume(s) > 0);
  }
  {
    Geometrium geo;
    Surface& s = geo.new_surface();
    bool threw = false;
    try { isosurface_cartesian(s, { 1, 2, 2, 0, 1, 0, 1, 0, 1 }, sphere, -1.0); }
    catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw && s.face_number() == 0);
  }
}

// An axis-aligned box from lo to hi, its faces oriented outwards.
static Surface& box(Geometrium& g, Vec3 lo, Vec3 hi)
{
  Surface& s = g.new_surface();
  VertexId v[8];
  for (int k = 0; k < 8; ++k)
    v[k] = g.add_vertex({ k & 1 ? hi.x : lo.x, k & 2 ? hi.y : lo.y, k & 4 ? hi.z : lo.z });
  auto edge = [&](VertexId a, VertexId b) { EdgeId e = g.find_edge(a, b); return e.valid() ? e : g.add_edge(a, b); };
  auto tri = [&](int a, int b, int c) { s.add_face(g.add_face(edge(v[a], v[b]), edge(v[b], v[c]), edge(v[c], v[a]))); };
  const int quads[6][4] = { {0, 2, 3, 1}, {4, 5, 7, 6}, {0, 1, 5, 4}, {2, 6, 7, 3}, {0, 4, 6, 2}, {1, 3, 7, 5} };
  for (auto& q : quads) { tri(q[0], q[1], q[2]); tri(q[0], q[2], q[3]); }
  return s;
}

static void test_predicates_and_bbtree()
{
  const Vec3 a{0, 0, 0}, b{1, 0, 0}, c{0, 1, 0};
  CHECK(orient3d(a, b, c, Vec3{0, 0, -1}) > 0 && orient3d(a, b, c, Vec3{0, 0, 1}) < 0);
  CHECK(orient3d(a, b, c, Vec3{0.3, 0.3, 0}) == 0);
  // Simulation of Simplicity never answers zero, and swapping two points
  // turns the answer over.
  const Vec3 d{0.3, 0.3, 0};
  const Vec3* p1[4] = { &a, &b, &c, &d };
  const Vec3* p2[4] = { &b, &a, &c, &d };
  const unsigned r1[4] = { 0, 1, 2, 3 }, r2[4] = { 1, 0, 2, 3 };
  CHECK(orient3d_sos(p1, r1) != 0 && orient3d_sos(p1, r1) == - orient3d_sos(p2, r2));

  Geometrium g;
  Surface& s1 = box(g, {0, 0, 0}, {1, 1, 1});
  Surface& s2 = box(g, {3, 0, 0}, {4, 1, 1});
  Surface& s3 = box(g, {0.5, 0.5, 0.5}, {2, 2, 2});
  std::size_t n12 = 0, n13 = 0, n11 = 0;
  const BBTree t1(s1), t2(s2), t3(s3);
  BBTree::traverse_overlapping(t1, t2, [&](FaceId, FaceId) { ++n12; });
  BBTree::traverse_overlapping(t1, t3, [&](FaceId f1, FaceId f2) { CHECK(s1.contains(f1) && s3.contains(f2)); ++n13; });
  BBTree::traverse_overlapping(t1, t1, [&](FaceId, FaceId) { ++n11; });
  CHECK(n12 == 0 && n13 > 0 && n11 >= 12);
}

static void test_boolean()
{
  // Boxes in general position: no vertex in the plane of a face of the other.
  const double overlap = 0.5 * 0.75 * 0.875;
  struct Case { const char* name; BooleanOp op1, op2; double volume; };
  for (const Case& c : { Case{ "union", BooleanOp::A_out_B, BooleanOp::B_out_A, 2 - overlap },
                         Case{ "inter", BooleanOp::A_in_B,  BooleanOp::B_in_A,  overlap },
                         Case{ "diff",  BooleanOp::A_out_B, BooleanOp::B_in_A,  1 - overlap } })
  {
    Geometrium g;
    Surface& a = box(g, {0, 0, 0}, {1, 1, 1});
    Surface& b = box(g, {0.5, 0.25, 0.125}, {1.5, 1.25, 1.125});
    CHECK(volume(a) > 0 && volume(b) > 0 && ! is_self_intersecting(a));
    Surface& r = g.new_surface();
    {
      SurfaceInter si(a, b);
      const auto check = si.check();
      CHECK(check.ok && check.closed && si.sos_ties() == 0 && ! si.edges().empty());
      si.boolean(r, c.op1);
      si.boolean(r, c.op2);
      // The difference keeps the inside of b, turned over.
      if (std::string_view(c.name) == "diff")
      {
        for (FaceId f : si.pieces_b().faces()) g.revert_face(f);
        for (FaceId f : b.faces())             g.revert_face(f);
      }
    }
    g.remove_surface(a);
    g.remove_surface(b);
    CHECK(is_closed(r) && is_manifold(r) && is_orientable(r));
    CHECK(std::abs(volume(r) - c.volume) < 1e-12);
    CHECK( ! is_self_intersecting(r));
    g.collect();
    CHECK(std::abs(volume(r) - c.volume) < 1e-12);
  }

  // Boxes sharing the plane x = 1: degenerate, decided by Simulation of
  // Simplicity, and still valid.
  {
    Geometrium g;
    Surface& a = box(g, {0, 0, 0}, {1, 1, 1});
    Surface& b = box(g, {1, 0.25, 0.25}, {2, 0.75, 0.75});
    Surface& r = g.new_surface();
    {
      SurfaceInter si(a, b);
      CHECK(si.sos_ties() > 0);
      const auto check = si.check();
      if (check.ok && check.closed)
      {
        si.boolean(r, BooleanOp::A_out_B);
        si.boolean(r, BooleanOp::B_out_A);
      }
    }
    if (r.face_number() > 0)
      CHECK(is_closed(r) && is_manifold(r) && is_orientable(r) && std::abs(volume(r) - 1.25) < 1e-12);
  }

  // Two boxes through each other in one surface intersect themselves.
  {
    Geometrium g;
    Surface& a = box(g, {0, 0, 0}, {1, 1, 1});
    Surface& b = box(g, {0.5, 0.25, 0.125}, {1.5, 1.25, 1.125});
    std::size_t pairs = 0;
    CHECK( ! foreach_intersecting_face(a, [&](FaceId, FaceId) { ++pairs; }) && pairs == 0);
    for (FaceId f : b.faces()) a.add_face(f);
    CHECK(foreach_intersecting_face(a, [&](FaceId, FaceId) { ++pairs; }) && pairs > 0);
  }
}

static void test_small_functions()
{
  Geometrium g;
  Surface& s = g.new_surface();
  Tet t = tetrahedron(g, s);
  // Face 0 is 0-2-1 with edges 02, 12, 01; face 1 is 0-1-3.
  CHECK(g.opposite_edge(t.f[0], t.v[0]) == t.e[3]);
  CHECK(g.opposite_edge(t.f[0], t.v[1]) == t.e[1]);
  CHECK( ! g.opposite_edge(t.f[0], t.v[3]).valid());
  CHECK(g.common_edge(t.f[0], t.f[1]) == t.e[0] && g.common_edge(t.f[1], t.f[0]) == t.e[0]);
  // The bottom face, 0-2-1, runs clockwise seen from +z.
  CHECK(g.orientation(t.f[0]) < 0);

  // gts_surface_copy: new elements, the same shape.
  Surface& c = g.new_surface();
  g.import(g, s, c);
  CHECK(c.face_number() == 4 && is_closed(c) && std::abs(volume(c) - volume(s)) < 1e-15);
  for (FaceId f : c.faces()) CHECK( ! s.contains(f));
  bool threw = false;
  try { g.import(g, s, s); } catch (const std::invalid_argument&) { threw = true; }
  CHECK(threw);

  // gts_surface_merge
  Surface& m = g.new_surface();
  m.merge(s);
  m.merge(c);
  CHECK(m.face_number() == 8 && t.f[0].i < g.face_capacity());

  std::ostringstream out;
  print_stats(s, out);
  CHECK(out.str().starts_with("# vertices: 4 edges: 6 faces: 4\n"));
}

static void test_cluster()
{
  Geometrium g;
  Surface& in = generate_sphere(g.new_surface(), 4);
  BBox b;
  for (VertexId v : in.vertices()) b.add(g.vertex(v).p);
  Surface& out = g.new_surface();
  Range r;
  const std::size_t nv0 = g.vertex_capacity();
  {
    ClusterGrid grid(out, b, 0.3);
    for (FaceId f : in.faces())
    {
      const auto v = g.vertices_of(f);
      grid.add_triangle(g.vertex(v[0]).p, g.vertex(v[1]).p, g.vertex(v[2]).p);
    }
    r = grid.update();
  }
  CHECK(out.face_number() > 0 && out.face_number() < in.face_number());
  CHECK(r.n > 0 && r.min >= 1 && r.n == g.vertex_capacity() - nv0);
  // Each vertex is the mean of points on the unit sphere, so inside it.
  for (VertexId v : out.vertices()) CHECK(norm(g.vertex(v).p) <= 1 + 1e-12 && norm(g.vertex(v).p) > 0.7);
  // Clusters that no face uses are removed with the grid.
  for (Index i = Index(nv0); i < g.vertex_capacity(); ++i)
    CHECK(g.dead(VertexId(i)) || g.vertex(VertexId(i)).edge.valid());

  bool threw = false;
  try { ClusterGrid bad(out, b, 0.); } catch (const std::invalid_argument&) { threw = true; }
  CHECK(threw);
}

static void test_sphere()
{
  for (unsigned order = 1; order <= 4; ++order)
  {
    Geometrium g;
    Surface& s = generate_sphere(g.new_surface(), order);
    const std::size_t nf = 20u << (2 * (order - 1));
    CHECK(s.face_number() == nf && s.edge_number() == nf * 3 / 2 && s.vertex_number() == nf / 2 + 2);
    CHECK(is_closed(s) && is_manifold(s) && is_orientable(s) && volume(s) > 0);
    for (VertexId v : s.vertices()) CHECK(std::abs(norm(g.vertex(v).p) - 1) < 1e-15);
    // No dead element is left in a ring, and every edge has two faces.
    for (EdgeId e : s.edges()) CHECK( ! g.dead(e) && count(g.faces_of(e)) == 2);
  }
  Geometrium g;
  bool threw = false;
  try { generate_sphere(g.new_surface(), 0); } catch (const std::invalid_argument&) { threw = true; }
  CHECK(threw);
}

static void test_tessellate()
{
  // Mid-edge tessellation of a tetrahedron keeps the area and the volume,
  // interpolates a vertex attribute, and leaves the face it splits in place.
  Geometrium g;
  auto w = g.add_vertex_attribute<double>("w", Hooks<double>().on_interpolate([](double a, double b, double t) { return a + t * (b - a); }));
  Surface& s = g.new_surface();
  Tet t = tetrahedron(g, s);
  for (int k = 0; k < 4; ++k) w[t.v[k]] = k;
  const double a0 = area(s), v0 = volume(s);
  tessellate(s, midvertex);
  CHECK(s.face_number() == 16 && s.edge_number() == 24 && s.vertex_number() == 10);
  CHECK(std::abs(area(s) - a0) < 1e-14 && std::abs(volume(s) - v0) < 1e-15);
  CHECK(is_closed(s) && is_orientable(s));
  for (FaceId f : t.f) CHECK(s.contains(f) && ! g.dead(f));
  for (EdgeId e : t.e) CHECK(g.dead(e));
  for (VertexId v : s.vertices())
  {
    if (v.i < 4) continue;
    // A midvertex has the mean weight of the two corners it lies between.
    double sum = 0;
    int k = 0;
    for (int i = 0; i < 4; ++i)
      for (int j = i + 1; j < 4; ++j)
        if (distance(0.5 * (g.vertex(t.v[i]).p + g.vertex(t.v[j]).p), g.vertex(v).p) < 1e-15) { sum = (i + j) / 2.; ++k; }
    CHECK(k == 1 && w[v] == sum);
  }

  // The default split is the midpoint: a tessellated cube stays a cube.
  {
    Geometrium c;
    Surface& cube = box(c, { 0, 0, 0 }, { 1, 1, 1 });
    tessellate(cube);
    CHECK(cube.face_number() == 48 && std::abs(area(cube) - 6) < 1e-14 && std::abs(volume(cube) - 1) < 1e-14);
  }

  // A face of another surface on an edge of s is split in two at the same
  // midvertex, so the two surfaces keep sharing the seam.
  Geometrium h;
  Surface& s1 = h.new_surface();
  Surface& s2 = h.new_surface();
  Tet u = tetrahedron(h, s1);
  const VertexId x = h.add_vertex({1, 1, -1});
  const FaceId outside = h.add_face(u.e[0], h.add_edge(u.v[1], x), h.add_edge(x, u.v[0]));
  s2.add_face(outside);
  const double a2 = area(s2);
  tessellate(s1, midvertex);
  CHECK(h.dead(u.e[0]) && s2.face_number() == 2 && std::abs(area(s2) - a2) < 1e-15);
  CHECK(s1.face_number() == 16 && is_closed(s1));
  unsigned shared = 0;
  for (EdgeId e : s2.edges()) if (face_number(s1, e) == 2) ++shared;
  CHECK(shared == 2); // the halves of the old edge
}

static void test_refine()
{
  Geometrium g;
  Surface& s = g.new_surface();
  tetrahedron(g, s);
  const double a0 = area(s), v0 = volume(s);
  double last = -1e300;
  std::size_t n_split = 0;
  // Costs come off the heap in increasing order: the longest edge first.
  refine(s, [&](double cost, std::size_t) { CHECK(cost >= last); last = cost; return n_split++ == 10; });
  // Each split of an edge between two faces adds a vertex, three edges and two faces.
  CHECK(s.vertex_number() == 14 && s.edge_number() == 36 && s.face_number() == 24);
  CHECK(std::abs(area(s) - a0) < 1e-14 && std::abs(volume(s) - v0) < 1e-15);
  CHECK(is_closed(s) && is_manifold(s) && is_orientable(s));
  // The longest edge left is no longer than the shortest edge split last.
  double longest = 0;
  for (EdgeId e : s.edges()) longest = std::max(longest, distance(g.segment(e).a, g.segment(e).b));
  CHECK(longest * longest <= -last + 1e-15);
}

static void test_swap_edge()
{
  // A square 0-1-2-3 cut along 0-2, both faces counterclockwise from +z.
  Geometrium g;
  Surface& s = g.new_surface();
  auto ids = g.add_face_attribute<int>("id");
  const VertexId v[4] = { g.add_vertex({0, 0, 0}), g.add_vertex({1, 0, 0}), g.add_vertex({1, 1, 0}), g.add_vertex({0, 1, 0}) };
  const EdgeId e01 = g.add_edge(v[0], v[1]), e12 = g.add_edge(v[1], v[2]), e23 = g.add_edge(v[2], v[3]),
               e30 = g.add_edge(v[3], v[0]), e02 = g.add_edge(v[0], v[2]);
  const FaceId f1 = g.add_face(e01, e12, e02), f2 = g.add_face(e02, e23, e30);
  CHECK(g.normal(f1).z > 0 && g.normal(f2).z > 0);
  ids[f1] = 7; ids[f2] = 8;
  s.add_face(f1); s.add_face(f2);

  swap_edge(s, e02);
  CHECK(s.face_number() == 2 && g.dead(e02) && g.dead(f1) && g.dead(f2));
  const EdgeId e13 = g.find_edge(v[1], v[3]);
  CHECK(e13.valid() && count(g.faces_of(e13)) == 2);
  int sum = 0;
  for (FaceId f : s.faces()) { CHECK(g.normal(f).z > 0); sum += ids[f]; }
  CHECK(sum == 15);
  CHECK(is_orientable(s) && std::abs(area(s) - 1) < 1e-15);

  // A third face on the edge, and an edge on the boundary, refuse the swap.
  const VertexId w = g.add_vertex({1, 1, 1});
  s.add_face(g.add_face(e13, g.add_edge(v[3], w), g.add_edge(w, v[1])));
  bool threw = false;
  try { swap_edge(s, e13); } catch (const std::invalid_argument&) { threw = true; }
  CHECK(threw);
  threw = false;
  try { swap_edge(s, e01); } catch (const std::invalid_argument&) { threw = true; }
  CHECK(threw);
}

static void test_distance_inside_curvature_merge()
{
  Geometrium g;
  Surface& s = g.new_surface();
  Tet t = tetrahedron(g, s);

  // Point to edge and face, and the closest points.
  CHECK(std::abs(distance(g, {0.5, -1, 0}, t.e[0]) - 1) < 1e-15);
  CHECK(closest_point(g, {2, 0, 0}, t.e[0]) == (Vec3{1, 0, 0}));
  CHECK(std::abs(distance(g, {0.2, 0.2, -2}, t.f[0]) - 2) < 1e-15);
  const Vec3 c = closest_point(g, {1, 1, 1}, t.f[3]);
  CHECK(std::abs(c.x - 1./3) < 1e-15 && std::abs(c.y - 1./3) < 1e-15 && std::abs(c.z - 1./3) < 1e-15);
  CHECK( ! intersection(g, {0.1, 0.1, 1}, {0.1, 0.1, 2}, t.f[0], false));
  const auto x = intersection(g, {0.1, 0.2, 1}, {0.1, 0.2, -1}, t.f[0], false);
  CHECK(x && std::abs(x->x - 0.1) < 1e-15 && std::abs(x->y - 0.2) < 1e-15 && x->z == 0);

  // Inside, also where the ray meets an edge or a vertex exactly.
  const BBTree tree(s);
  CHECK(is_inside(tree, {0.1, 0.1, 0.1}) && ! is_inside(tree, {1, 1, 1}));
  CHECK(is_inside(tree, {0.1, 0.2, 0.3}, true) == false);
  CHECK(is_inside(tree, {-1, 0.25, 0.25}) == false);  // through the left face's middle
  CHECK(is_inside(tree, {0.2, 0.0, 0.0}) == is_inside(tree, {0.2, 0.0, 0.0})); // on an edge: repeatable
  CHECK(is_inside(tree, {-1, 0, 0}) == false);        // the ray runs along an edge
  CHECK(std::abs(tree.point_distance({0.2, 0.2, -2}) - 2) < 1e-15);
  double d2;
  const Vec3 cp = tree.point_closest({2, 0, 0}, {}, &d2);
  CHECK(cp == (Vec3{1, 0, 0}) && d2 == 1);

  // Two concentric spheres: the distance is the difference of the radii.
  Surface& a = g.new_surface();
  Surface& b = g.new_surface();
  generate_sphere(a, 4);
  generate_sphere(b, 4);
  for (VertexId v : b.vertices()) g.position(v) = 1.5 * g.vertex(v).p;
  const SurfaceDistance sd = distance(b, a, 0.05);
  // The inside of the outer sphere's flat faces is a little closer.
  CHECK(sd.faces.min > 0.49 && std::abs(sd.faces.max - 0.5) < 1e-12 && sd.boundary.n == 0);

  // Curvature of the unit sphere: mean 1, Gaussian 1, Kh along the inward
  // normal.
  for (VertexId v : a.vertices())
  {
    const auto kh = mean_curvature_normal(a, v);
    const auto kg = gaussian_curvature(a, v);
    CHECK(kh && kg);
    if ( ! kh || ! kg) break;
    const Vec3 p = g.vertex(v).p;
    CHECK(std::abs(norm(*kh) / 2 - 1) < 5e-3 && dot(*kh, p) < 0 && std::abs(*kg - 1) < 1e-2);
    const auto [e1, e2] = principal_directions(a, v, *kh, *kg);
    CHECK(std::abs(dot(e1, p)) < 1e-2 && std::abs(dot(e1, e2)) < 1e-12);
  }

  // Centres and boundary.
  const CenterOfMass cm = center_of_mass(s);
  CHECK(std::abs(cm.volume - 1./6) < 1e-15 && std::abs(cm.center.x - 0.25) < 1e-15 && std::abs(cm.center.z - 0.25) < 1e-15);
  CHECK(std::abs(center_of_area(a).center.x) < 1e-15 && boundary(a).empty());

  // Faces on vertices of their own, merged back into a closed tetrahedron.
  Geometrium h;
  Surface& u = h.new_surface();
  std::vector<VertexId> vs;
  std::vector<EdgeId> es;
  for (FaceId f : s.faces())
  {
    const auto o = g.vertices_of(f);
    VertexId v[3];
    for (int k = 0; k < 3; ++k) vs.push_back(v[k] = h.add_vertex(g.vertex(o[k]).p + Vec3{1e-12 * k, 0, 0}));
    EdgeId e[3];
    for (int k = 0; k < 3; ++k) es.push_back(e[k] = h.add_edge(v[k], v[(k + 1) % 3]));
    u.add_face(h.add_face(e[0], e[1], e[2]));
  }
  CHECK(KdTree(h, vs).range({ {-0.1, -0.1, -0.1}, {0.1, 0.1, 0.1} }).size() == 3);
  CHECK( ! is_closed(u));
  vs = h.merge_vertices(vs, 1e-9);
  es = h.merge_edges(es);
  CHECK(vs.size() == 4 && es.size() == 6 && is_closed(u) && is_orientable(u));
  CHECK(u.vertex_number() == 4 && u.edge_number() == 6);

  // A merge across an edge removes the edge and its faces: with vertex 1
  // next to vertex 0, the bottom and front faces go, and the other two
  // become duplicates on duplicate edges.
  Geometrium k;
  Surface& w = k.new_surface();
  const Tet tk = tetrahedron(k, w);
  k.position(tk.v[1]) = { 1e-12, 0, 0 };
  std::vector<VertexId> all(tk.v.begin(), tk.v.end());
  all = k.merge_vertices(all, 1e-9);
  CHECK(all.size() == 3 && k.dead(tk.v[1]) && k.dead(tk.e[0]) && k.dead(tk.f[0]) && k.dead(tk.f[1]));
  CHECK(w.face_number() == 2 && k.merge_edges({ tk.e[1], tk.e[2], tk.e[3], tk.e[4], tk.e[5] }).size() == 3);
  CHECK(remove_duplicate_faces(w) == 1 && w.face_number() == 1 && remove_duplicate_faces(w) == 0);
}

// Errors of GTS that gts-cxx fixes; see "Fixed" in QUIRKS.md.
static void test_fixes()
{
  // The bounded and bcl tetrahedra close the surface, also on a grid that
  // is not square. The sphere crosses the first and the last layer of cubes
  // in z.
  const CartesianGrid grid = { 13, 11, 9, -0.75, 0.125, -0.625, 0.125, -0.5, 0.125 };
  auto ball = sample([](Vec3 p) { return p.x * p.x + p.y * p.y + p.z * p.z; });
  for (auto method : { isosurface_tetra_bounded, isosurface_tetra_bcl })
  {
    Geometrium geo;
    Surface& s = geo.new_surface();
    method(s, grid, ball, 0.16);
    CHECK(s.face_number() > 0 && is_closed(s) && is_manifold(s) && is_orientable(s));
    CHECK(std::isfinite(area(s)));
  }

  // refine() costs a new edge once it has its faces.
  {
    Geometrium g;
    Surface& s = g.new_surface();
    tetrahedron(g, s);
    unsigned costed = 0, without_faces = 0, splits = 0;
    refine(s, [&](double, std::size_t) { return splits++ == 10; },
           [&](const Geometrium& g, EdgeId e) {
             ++costed;
             if (face_number(s, e) != 2) ++without_faces;
             const Segment x = g.segment(e);
             return - distance2(x.a, x.b);
           });
    CHECK(costed > 6 && without_faces == 0);
  }

  // Edge angles and boundaries are those within the surface.
  {
    Geometrium g;
    Surface& s = g.new_surface();
    Surface& t = g.new_surface();
    const Tet a = tetrahedron(g, s);
    t.add_face(a.f[0]);
    const VertexId far = g.add_vertex({ 0.3, 0.3, -1 });
    const EdgeId e1 = g.add_edge(a.v[0], far), e2 = g.add_edge(a.v[1], far);
    s.add_face(g.add_face(a.e[0], e2, e1)); // a fin on edge 0-1, so it has three faces in s
    CHECK(quality_stats(t).edge_angle.n == 0);

    Surface& u = g.new_surface();
    const VertexId w0 = g.add_vertex({ 5, 0, 0 }), w1 = g.add_vertex({ 6, 0, 0 }), w2 = g.add_vertex({ 5, 1, 0 });
    u.add_face(g.add_face(g.add_edge(w0, w1), g.add_edge(w1, w2), g.add_edge(w2, w0)));
    // t's edges have faces of s too, but t's single face makes them t's boundary.
    CHECK(distance(t, u, 0.1).boundary.n > 0);
  }

  // The ray of is_inside() leaves the box also when the box ends at x = 0.
  {
    Geometrium g;
    const BBTree tree(box(g, { -2, -1, -1 }, { 0, 1, 1 }));
    CHECK(is_inside(tree, { -1, 0, 0 }) && ! is_inside(tree, { -1, 2, 0 }));
  }

  // A collapse that leaves an edge without faces, through a duplicate face,
  // keeps the new vertex while it has other edges. The edges are made so
  // that the faceless edge comes last in the new vertex's ring.
  {
    Geometrium g;
    Surface& s = g.new_surface();
    const VertexId a = g.add_vertex({ 0, 0, 0 }), b = g.add_vertex({ 1, 0, 0 }), c = g.add_vertex({ 0.5, 1, 0 });
    const VertexId d = g.add_vertex({ 2, -1, 0 }), f = g.add_vertex({ 2, 1, 0 });
    const EdgeId bd = g.add_edge(b, d), df = g.add_edge(d, f), fb = g.add_edge(f, b);
    const EdgeId ab = g.add_edge(a, b), ca = g.add_edge(c, a), bc = g.add_edge(b, c);
    s.add_face(g.add_face(bd, df, fb));
    s.add_face(g.add_face(ab, bc, ca));
    s.add_face(g.add_face(ab, bc, ca));
    const VertexId m = g.collapse_edge(ab, { 0.5, 0, 0 });
    CHECK(m.valid() && count(g.edges_of(m)) == 2 && s.face_number() == 1);
  }

  // Edges of equal cost are taken in an order that does not depend on the
  // order of the faces: a symmetric sphere, with its faces given in reverse,
  // coarsens and refines to the same vertices.
  {
    auto reversed_copy = [](const Surface& s, Geometrium& h) -> Surface& {
      const Geometrium& g = s.geometrium();
      Surface& r = h.new_surface();
      std::vector<FaceId> faces(s.faces().begin(), s.faces().end());
      std::vector<VertexId> vmap(g.vertex_capacity());
      std::vector<EdgeId> emap(g.edge_capacity());
      for (auto it = faces.rbegin(); it != faces.rend(); ++it)
      {
        std::array<EdgeId, 3> fe;
        for (int k = 0; k < 3; ++k)
        {
          const EdgeId e = g.face(*it).e[2 - k];
          if ( ! emap[e.i].valid())
          {
            for (VertexId v : g.edge(e).v)
              if ( ! vmap[v.i].valid()) vmap[v.i] = h.add_vertex(g.vertex(v).p);
            emap[e.i] = h.add_edge(vmap[g.edge(e).v[1].i], vmap[g.edge(e).v[0].i]);
          }
          fe[k] = emap[e.i];
        }
        r.add_face(h.add_face(fe[0], fe[1], fe[2]));
      }
      return r;
    };
    auto points = [](const Surface& s) {
      std::vector<std::array<double, 3>> p;
      for (VertexId v : s.vertices()) { const Vec3 x = s.geometrium().vertex(v).p; p.push_back({ x.x, x.y, x.z }); }
      std::sort(p.begin(), p.end());
      return p;
    };
    Geometrium g1, g2, g3, g4;
    Surface& a = generate_sphere(g1.new_surface(), 3);
    Surface& b = reversed_copy(a, g2);
    coarsen(a, stop_below_edges(200));
    coarsen(b, stop_below_edges(200));
    CHECK(a.edge_number() == b.edge_number() && points(a) == points(b));
    Surface& c = generate_sphere(g3.new_surface(), 2);
    Surface& d = reversed_copy(c, g4);
    refine(c, [](double, std::size_t n) { return n > 600; });
    refine(d, [](double, std::size_t n) { return n > 600; });
    CHECK(c.edge_number() == d.edge_number() && points(c) == points(d));
  }

  // The mean curvature has a sign: positive on a sphere whose faces face
  // outwards, negative once they are reverted, and then both principal
  // curvatures are negative.
  {
    Geometrium g;
    Surface& s = generate_sphere(g.new_surface(), 4);
    const VertexId v = *s.vertices().begin();
    const auto h = mean_curvature(s, v);
    const auto k = gaussian_curvature(s, v);
    CHECK(h && k && std::abs(*h - 1) < 0.02);
    for (FaceId f : s.faces()) g.revert_face(f);
    const auto hr = mean_curvature(s, v);
    CHECK(hr && std::abs(*hr + *h) < 1e-12);
    const auto [k1, k2] = principal_curvatures(*hr, *k);
    CHECK(k1 < 0 && k2 < 0 && k1 >= k2);
  }

  // On a tetrahedron the normals around an edge's ends sum to rounding
  // noise. That row is dropped as nearly zero for its terms' size, at any
  // scale; GTS took its direction as a constraint and put the vertex some
  // 1e15 away.
  for (double scale : { 1.0, 1e-3, 1e3 })
  {
    Geometrium g;
    Surface& s = g.new_surface();
    const VertexId v0 = g.add_vertex(scale * Vec3{ 0.1, 0.2, 0.3 }), v1 = g.add_vertex(scale * Vec3{ 1.3, 0.1, 0.2 });
    const VertexId v2 = g.add_vertex(scale * Vec3{ 0.2, 1.1, 0.4 }), v3 = g.add_vertex(scale * Vec3{ 0.3, 0.2, 1.7 });
    const EdgeId e01 = g.add_edge(v0, v1), e02 = g.add_edge(v0, v2), e03 = g.add_edge(v0, v3);
    const EdgeId e12 = g.add_edge(v1, v2), e13 = g.add_edge(v1, v3), e23 = g.add_edge(v2, v3);
    s.add_face(g.add_face(e02, e12, e01));
    s.add_face(g.add_face(e01, e13, e03));
    s.add_face(g.add_face(e03, e23, e02));
    s.add_face(g.add_face(e12, e23, e13));
    const Vec3 p = volume_optimized_vertex(g, e01);
    CHECK(distance(p, scale * Vec3{ 0.475, 0.4, 0.65 }) < 1e-9 * scale);
  }

  // A comment at the end of a line ends the line.
  {
    Geometrium g;
    std::istringstream in("# a closed tetrahedron\n4 6 4\n0 0 0 # origin\n1 0 0\n0 1 0\n0 0 1\n"
                          "1 2\n1 3\n1 4 ! last edge from 1\n2 3\n2 4\n3 4\n2 4 1\n1 5 3\n3 6 2\n4 6 5\n");
    auto r = read_gts(g, in);
    CHECK(r.has_value() && (*r)->face_number() == 4 && is_closed(**r));
  }
}

//==============================================================================
// Ported on 2026-10-11: neighbourhoods, Delaunay, matrices, boxes, writers,
// strips, parked elements, progressive surfaces, graphs, the slice API
//==============================================================================

static void test_topology()
{
  Geometrium g;
  Surface& s = g.new_surface();
  Tet t = tetrahedron(g, s);
  CHECK(vertex_faces(s, t.v[0]).size() == 3 && vertex_neighbors(g, t.v[0], &s).size() == 3);
  CHECK(fan_oriented(s, t.v[0]).size() == 3);
  CHECK(face_neighbors(g, t.f[0], &s).size() == 3 && is_compatible(s, t.f[0]));
  CHECK(manifold_faces(s, t.e[0]).has_value() && belongs_to_tetrahedron(g, t.e[0]));
  CHECK(edge_contacts(g, t.e[0]) == 1);
  CHECK(edges_from_vertices(g, std::array{ t.v[0] }).size() == 3);
  CHECK(faces_from_edges(g, std::array{ t.e[0], t.e[5] }).size() == 4);
  CHECK(vertices_from_edges(g, std::array{ t.e[0], t.e[5] }).size() == 4);

  // Breadth first from a face: itself, then the three others.
  SurfaceTraverse tr(s, t.f[0]);
  std::vector<unsigned> levels;
  while (auto n = tr.next()) levels.push_back(n->second);
  CHECK(levels == (std::vector<unsigned>{ 0, 1, 1, 1 }));

  // Two tetrahedra sharing a vertex: one surface, two components, a
  // contact vertex that sever separates.
  Tet u = tetrahedron(g, s, { 0, 0, -1 });
  const VertexId shared = u.v[3];
  g.replace_vertex(t.v[0], shared);
  CHECK(g.vertex_contacts(shared) == 2);
  const std::size_t nv = s.vertex_number();
  CHECK(g.vertex_contacts(shared, true) == 2 && s.vertex_number() == nv + 1);
  CHECK(g.vertex_contacts(shared) == 1);
  const auto parts = split_components(s);
  CHECK(parts.size() == 2 && parts[0]->face_number() == 4 && parts[1]->face_number() == 4);
  for (Surface* p : parts) CHECK(is_closed(*p) && is_orientable(*p));
  for (Surface* p : parts) g.remove_surface(*p);

  CHECK(remove_faces_if(s, [&](FaceId f) { return f == u.f[0] || f == u.f[1]; }) == 2 && s.face_number() == 6);
  // gts_object_destroy on a vertex takes its faces with it.
  g.destroy(t.v[1]);
  CHECK(s.face_number() == 3 && g.dead(t.v[1]) && g.dead(t.e[0]) && g.dead(t.f[0]));
}

static void test_delaunay()
{
  // The planar tests.
  const Triangle tr{ { 0, 0, 0 }, { 2, 0, 0 }, { 0, 2, 0 } };
  CHECK(point_in_triangle({ 0.5, 0.5, 9 }, tr) == Intersect::In);
  CHECK(point_in_triangle({ 1, 0, 0 }, tr) == Intersect::On && point_in_triangle({ 3, 3, 0 }, tr) == Intersect::Out);
  CHECK(in_triangle_circle({ 1, 1, 0 }, tr) > 0 && in_triangle_circle({ 3, 3, 0 }, tr) < 0);
  CHECK(in_triangle_circle({ 2, 2, 0 }, tr) == 0); // on the circle
  const auto c = circumcircle_center(tr);
  CHECK(c && c->x == 1 && c->y == 1 && c->z == 0);
  CHECK( ! circumcircle_center({ { 0, 0, 0 }, { 1, 1, 0 }, { 2, 2, 0 } }));
  CHECK(segments_intersect({ { 0, 0, 0 }, { 2, 2, 0 } }, { { 0, 2, 0 }, { 2, 0, 0 } }) == Intersect::In);
  CHECK(segments_intersect({ { 0, 0, 0 }, { 1, 1, 0 } }, { { 1, 1, 0 }, { 2, 0, 0 } }) == Intersect::On);
  CHECK(segments_intersect({ { 0, 0, 0 }, { 1, 0, 0 } }, { { 0, 1, 0 }, { 1, 1, 0 } }) == Intersect::Out);
  {
    const double a[3] = { 0, 0, 0 }, b[3] = { 1, 0, 0 }, cc[3] = { 0, 1, 0 }, d[3] = { 0, 0, 1 }, e[3] = { 0.1, 0.1, 0.1 };
    CHECK(orient3d(a, b, cc, d) < 0 && insphere(a, cc, b, d, e) > 0); // a, c, b, d has a positive orientation
  }

  // A triangulation of points on a parabola-free scatter, from an enclosing
  // triangle.
  Geometrium g;
  std::vector<Vec3> pts;
  unsigned x = 12345;
  auto rnd = [&] { x = x * 1103515245u + 12345u; return double((x >> 8) & 0xffff) / 65536.; };
  for (int i = 0; i < 200; ++i) { const double px = rnd(), py = rnd(); pts.push_back({ px, py, px + 2 * py }); }
  Surface& s = g.new_surface();
  const FaceId t = enclosing_triangle(g, pts, 100.);
  const auto hull = g.vertices_of(t);
  s.add_face(t);
  std::vector<VertexId> vs;
  for (Vec3 p : pts) { vs.push_back(g.add_vertex(p)); CHECK( ! delaunay_add_vertex(s, vs.back()).valid()); }
  CHECK( ! delaunay_check(s).valid());
  // The same x and y again: refused, the vertex there returned.
  CHECK(delaunay_add_vertex(s, g.add_vertex({ pts[7].x, pts[7].y, 5 })) == vs[7]);
  // A constraint across the square.
  const EdgeId ce = g.add_edge(vs[0], vs[1]);
  CHECK(delaunay_add_constraint(s, ce).empty());
  CHECK(is_constraint(g, ce) && count(g.faces_of(ce)) == 2);
  {
    KeepOrphans keep(g, { .vertices = true });
    for (VertexId h : hull) g.destroy(h);
  }
  CHECK(s.vertex_number() == 200 && is_orientable(s) && is_manifold(s));
  // The plane z = x + 2y is interpolated exactly.
  const FaceId f0 = *s.faces().begin();
  const Triangle q = g.triangle(f0);
  const Vec3 m{ (q.a.x + q.b.x + q.c.x) / 3, (q.a.y + q.b.y + q.c.y) / 3, 0 };
  CHECK(std::abs(interpolate_height(g, f0, m) - (m.x + 2 * m.y)) < 1e-12);
  // Removing a vertex keeps it Delaunay.
  delaunay_remove_vertex(s, vs[50]);
  CHECK(s.vertex_number() == 199 && is_orientable(s) && is_manifold(s));
  // Refinement of a square with constrained sides and points inside, a
  // domain where Ruppert's refinement ends. (On the hull of the random
  // points it does not: GTS and gts-cxx insert circumcentres ever closer to
  // a vertex until the circumcircle cannot be computed.)
  {
    Geometrium q;
    std::vector<Vec3> sq = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 } };
    for (int i = 0; i < 30; ++i) sq.push_back({ 0.1 + 0.8 * rnd(), 0.1 + 0.8 * rnd(), 0 });
    Surface& d = q.new_surface();
    const FaceId et = enclosing_triangle(q, sq, 100.);
    const auto eh = q.vertices_of(et);
    d.add_face(et);
    std::vector<VertexId> v;
    for (Vec3 p : sq) { v.push_back(q.add_vertex(p)); delaunay_add_vertex(d, v.back()); }
    for (int k = 0; k < 4; ++k) delaunay_add_constraint(d, q.add_edge(v[k], v[(k + 1) % 4]));
    {
      KeepOrphans keep(q, { .vertices = true });
      for (VertexId h : eh) q.destroy(h);
    }
    delaunay_remove_hull(d);
    CHECK(std::abs(area(d) - 1) < 1e-12);
    CHECK(delaunay_conform(d, -1) == 0);
    const double qmin = 0.5, amax = 0.002;
    const FaceCost cost = [&](const Geometrium& g, FaceId f) {
      return g.quality(f) < qmin || g.area(f) > amax ? g.quality(f) : 0.;
    };
    CHECK(delaunay_refine(d, -1, cost) == 0);
    for (FaceId f : d.faces()) CHECK(q.area(f) <= amax && q.quality(f) >= qmin);
    CHECK(std::abs(area(d) - 1) < 1e-12 && is_orientable(d) && is_manifold(d));
  }
}

static void test_matrix()
{
  const Matrix4 r = Matrix4::rotate({ 0, 0, 2 }, std::numbers::pi / 2);
  const Vec3 p = transform(r, { 1, 0, 0 });
  CHECK(std::abs(p.x) < 1e-15 && std::abs(p.y - 1) < 1e-15 && p.z == 0);
  const Matrix4 m = Matrix4::translate({ 1, 2, 3 }) * r * Matrix4::scale({ 2, 3, 4 });
  CHECK(std::abs(determinant(m) - 24) < 1e-12);
  const auto mi = inverse(m);
  CHECK(mi.has_value());
  const Matrix4 id = *mi * m;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) CHECK(std::abs(id[i][j] - (i == j)) < 1e-14);
  CHECK( ! inverse(Matrix4{}).has_value() && ! inverse3(Matrix4{}).has_value());
  CHECK(transpose(transpose(m)) == m);
  Geometrium g;
  Surface& s = g.new_surface();
  Tet t = tetrahedron(g, s);
  const Matrix4 pr = Matrix4::projection(g, t.f[3]);
  const Matrix4 ptp = transpose(pr) * pr;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) CHECK(std::abs(ptp[i][j] - (i == j)) < 1e-14);
  const double v0 = volume(s);
  transform(s, Matrix4::scale({ 2, 2, 2 }));
  CHECK(std::abs(volume(s) - 8 * v0) < 1e-14);
  CHECK(to_string(Matrix4::identity()).starts_with("[[              1"));
}

static void test_boxes_and_writers()
{
  BBox b;
  b.add({ 0, 0, 0 });
  b.add({ 1, 1, 1 });
  CHECK(overlaps(b, Triangle{ { 0.5, 0.5, -1 }, { 0.5, 0.5, 2 }, { 3, 3, 0.5 } }));
  CHECK( ! overlaps(b, Triangle{ { 2, 2, 2 }, { 3, 2, 2 }, { 2, 3, 2 } }));
  CHECK(overlaps(b, Segment{ { -1, 0.5, 0.5 }, { 2, 0.5, 0.5 } }) && ! overlaps(b, Segment{ { 2, 0, 0 }, { 3, 0, 0 } }));
  CHECK(is_stabbed(b, { -5, 0.5, 0.5 }) && ! is_stabbed(b, { 5, 0.5, 0.5 }));

  Geometrium g;
  Surface& s = g.new_surface();
  tetrahedron(g, s);
  tetrahedron(g, s, { 5, 0, 0 });
  const BBTree tree(s);
  BBox q;
  q.add({ 4.9, -0.1, -0.1 });
  q.add({ 5.2, 0.2, 0.2 });
  CHECK(tree.overlap(q).size() == 4 && tree.is_overlapping(q));
  std::ostringstream boxes;
  tree.write_oogl(1, boxes);
  CHECK(boxes.str().starts_with("{ LIST\n{ OFF 8 6 12\n0 0 0\n6 0 0\n") && boxes.str().ends_with("}\n}\n"));

  std::ostringstream off, vtk, bnd;
  write_oogl(s, off);
  CHECK(off.str().starts_with("OFF 8 8 12\n0 0 0\n1 0 0\n"));
  write_vtk(s, vtk);
  CHECK(vtk.str().find("POINTS 8 float\n") != std::string::npos && vtk.str().find("POLYGONS 8 32\n") != std::string::npos);
  write_oogl_boundary(s, bnd);
  CHECK(bnd.str() == "LIST {\n}\n");
  std::ostringstream coff;
  write_oogl(s, coff, { .vertex = [](VertexId) { return Rgb{ 1, 0, 0 }; } });
  CHECK(coff.str().starts_with("COFF 8 8 12\n0 0 0 1 0 0 1.0\n"));
}

static void test_strips()
{
  Geometrium g;
  Surface& s = g.new_surface();
  generate_sphere(s, 3);
  const auto list = strips(s);
  std::vector<int> seen(g.face_capacity());
  std::size_t n = 0;
  for (const auto& strip : list)
  {
    for (std::size_t i = 0; i < strip.size(); ++i)
    {
      ++seen[strip[i].i];
      ++n;
      if (i + 1 < strip.size()) CHECK(g.common_edge(strip[i], strip[i + 1]).valid());
    }
  }
  CHECK(n == s.face_number() && std::ranges::all_of(s.faces(), [&](FaceId f) { return seen[f.i] == 1; }));
  CHECK(list.size() < s.face_number() / 4);
}

static void test_park()
{
  Geometrium g;
  Surface& s = g.new_surface();
  Tet t = tetrahedron(g, s);
  auto colour = g.add_face_attribute<int>("colour");
  colour[t.f[3]] = 7;
  const FaceId f = t.f[3];
  const auto e = g.face(f).e;
  g.park(f);
  CHECK(g.dead(f) && g.parked(f) && ! s.contains(f) && count(g.faces_of(e[0])) == 1);
  // Collection keeps the parked face, renumbered, with its attribute.
  g.add_vertex({ 9, 9, 9 }); // an orphan, which goes
  FaceId pf = f;
  auto sub = g.on_renumber([&](const Renumbering& r) { pf = r(pf); });
  g.collect();
  CHECK(g.parked(pf) && colour[pf] == 7 && g.vertex_capacity() == 4);
  // Back on the edges of the hole, those with one face.
  std::array<EdgeId, 3> ne;
  int k = 0;
  for (Index i = 0; i < g.edge_capacity(); ++i)
    if (count(g.faces_of(EdgeId(i))) == 1) ne[k++] = EdgeId(i);
  CHECK(k == 3);
  g.unpark(pf, ne[0], ne[1], ne[2]);
  s.add_face(pf);
  CHECK( ! g.dead(pf) && is_closed(s) && colour[pf] == 7);
  // A vertex without edges.
  const VertexId v = g.add_vertex({ 1, 2, 3 });
  g.park(v);
  g.collect();
  CHECK(g.vertex_capacity() == 5);
}

static void test_progressive()
{
  Geometrium g;
  Surface& s = g.new_surface();
  generate_sphere(s, 4);
  const double a0 = area(s), v0 = volume(s);
  const std::size_t n0 = s.vertex_number();
  PSurface ps(s, stop_below_edges(60), {}, {}, std::numbers::pi / 180);
  CHECK(ps.max_vertex_number() == n0 && ps.vertex_number() == ps.min_vertex_number() && s.vertex_number() == ps.min_vertex_number());
  CHECK(is_closed(s) && is_manifold(s) && is_orientable(s));
  ps.set_vertex_number(300);
  CHECK(s.vertex_number() == 300 && is_closed(s) && is_orientable(s));
  g.collect(); // the splits follow
  ps.set_vertex_number(ps.min_vertex_number() + 10);
  ps.set_vertex_number(ps.max_vertex_number());
  CHECK(s.vertex_number() == n0 && area(s) == a0 && volume(s) == v0);

  // Written and read back.
  std::stringstream io;
  ps.write(io);
  Geometrium h;
  auto r = PSurface::read(h, io);
  CHECK(r.has_value());
  if (r)
  {
    PSurface& p2 = **r;
    CHECK(p2.max_vertex_number() == n0 && p2.min_vertex_number() == ps.min_vertex_number());
    CHECK(std::abs(area(p2.surface()) - a0) < 1e-8);
    p2.set_vertex_number(200);
    ps.set_vertex_number(200);
    CHECK(std::abs(area(p2.surface()) - area(s)) < 1e-8 && p2.surface().face_number() == s.face_number());
    ps.set_vertex_number(ps.max_vertex_number());
  }

  // A hierarchical surface, its splits taken in a scrambled order.
  auto key = [&](const VertexSplit& vs) { return g.vertex(vs.v).p.y; };
  HSurface hs(ps, key, key);
  CHECK(hs.vertex_number() == n0 && ps.splits().empty() && hs.height() > 1);
  unsigned x = 7;
  for (int step = 0; step < 3000; ++step)
  {
    x = x * 1103515245u + 12345u;
    const int i = int((x >> 8) % hs.splits().size());
    if (hs.can_expand(i)) hs.expand(i);
    else if (hs.can_collapse(i)) hs.collapse(i);
  }
  CHECK(is_closed(s) && is_manifold(s) && is_orientable(s) && s.vertex_number() == hs.vertex_number());
  while (auto top = hs.top_expandable()) hs.expand(top->first);
  CHECK(s.vertex_number() == n0 && std::abs(area(s) - a0) < 1e-12 && std::abs(volume(s) - v0) < 1e-12);
  int visited = 0;
  hs.traverse(true, -1, [&](int) { ++visited; return false; });
  CHECK(visited == int(hs.splits().size()));
  hs.foreach(true, [](int) { return false; });
  CHECK(s.vertex_number() == n0 - hs.splits().size());
}

static void test_graph()
{
  Geometrium geo;
  Surface& s = geo.new_surface();
  generate_sphere(s, 3);
  Graph gr;
  NodeSet g = surface_graph(gr, s);
  CHECK(g.size() == s.face_number() && edges_of(g).size() == s.edge_number() && g.weight() == float(g.size()));
  for (GNodeId n : g.nodes()) CHECK(degree(gr, n, &g) == 3);
  CHECK(edges_cut(g) == 0 && distance_sum(g, g.nodes()[0]) > 0);

  auto covers = [&](std::span<const NodeSet> parts) {
    std::size_t n = 0;
    for (const NodeSet& p : parts) n += p.size();
    for (GNodeId x : g.nodes())
    {
      int k = 0;
      for (const NodeSet& p : parts) k += p.contains(x);
      if (k != 1) return false;
    }
    return n == g.size();
  };
  const auto bubbles = bubble_partition(g, 4, 10);
  CHECK(bubbles.size() == 4 && covers(bubbles));
  GraphBisection bg = ggg_bisection(g, 10);
  CHECK(check(bg) && bg.g1.size() + bg.g2.size() == g.size());
  while (bkl_refine(bg, 50, 0.1f) != 0.) {}
  CHECK(check(bg) && std::abs(bg.g1.weight() - bg.g2.weight()) <= 0.1f * g.weight());
  NodeSet whole = g;
  const auto parts = recursive_bisection(gr, whole, 2, 10, 50, 100, 0.1f);
  CHECK(parts.size() == 4 && covers(parts));
  for (const NodeSet& p : parts) CHECK(p.size() > g.size() / 8);
  Surface& part = graph_surface(parts[0], geo);
  CHECK(part.face_number() == parts[0].size());

  {
    PGraph pg(gr, g, 50);
    CHECK(g.size() <= 50 && pg.min_node_number() == g.size() && pg.max_node_number() == s.face_number());
    CHECK(g.weight() == float(s.face_number()));
    std::size_t expanded = 0;
    while (pg.down([&](const GNodeSplit&) { ++expanded; })) {}
    CHECK(g.size() == s.face_number() && expanded == pg.splits().size() && edges_of(g).size() == s.edge_number());
  }

  std::stringstream io;
  write_graph(g, io);
  Graph gr2;
  auto g2 = read_graph(gr2, io);
  CHECK(g2 && g2->size() == g.size() && edges_of(*g2).size() == edges_of(g).size());

  std::istringstream jostle("4 2\n2\n1 3\n2\n\n");
  Graph gr3;
  auto g3 = read_jostle(gr3, jostle);
  CHECK(g3 && g3->size() == 4 && edges_of(*g3).size() == 2);
  std::istringstream bad("3 1\n2\n1\n");
  Graph gr4;
  CHECK( ! read_jostle(gr4, bad).has_value());
}

static void test_isoslice()
{
  // The slice API on a cartesian grid makes the surface isosurface_cartesian
  // makes.
  const unsigned n = 20;
  const CartesianGrid cg{ n + 1, n + 1, n + 1, -1.25, 2.5 / n, -1.25, 2.5 / n, -1.25, 2.5 / n };
  auto f = [](Vec3 p) { return p.x * p.x + p.y * p.y + p.z * p.z; };
  Geometrium g;
  Surface& s1 = g.new_surface();
  isosurface_cartesian(s1, cg, sample(f), 1.);
  Surface& s2 = g.new_surface();
  const unsigned m = n + 1;
  GridPlane planes[2] = { GridPlane(m, m), GridPlane(m, m) };
  PlaneValues values[2] = { PlaneValues(m, m), PlaneValues(m, m) };
  auto fill = [&](int w, unsigned k) {
    for (unsigned i = 0; i < m; ++i)
      for (unsigned j = 0; j < m; ++j)
      {
        planes[w](i, j) = grid_point(cg, i, j, k);
        values[w](i, j) = f(planes[w](i, j));
      }
  };
  IsoSlice a(m, m), b(m, m);
  IsoSlice* sl[2] = { &a, &b };
  int x = 0, y = 1;
  fill(x, 0);
  fill(y, 1);
  sl[0]->fill(planes[x], &planes[y], values[x], &values[y], 1.);
  for (unsigned k = 2; k < m; ++k)
  {
    std::swap(x, y);
    fill(y, k);
    sl[1]->fill(planes[x], &planes[y], values[x], &values[y], 1.);
    isosurface_slice(s2, *sl[0], *sl[1]);
    std::swap(sl[0], sl[1]);
  }
  sl[1]->fill(planes[y], nullptr, values[y], nullptr, 1.);
  isosurface_slice(s2, *sl[0], *sl[1]);
  CHECK(s1.face_number() == s2.face_number() && std::abs(area(s1) - area(s2)) < 1e-12 && is_closed(s2));
}

// Whether two surfaces have the same faces, in the same order, through the
// same vertices.
static bool same_faces(const Surface& a, const Surface& b)
{
  if (a.face_number() != b.face_number()) return false;
  std::vector<std::array<Vec3, 3>> fa, fb;
  for (auto [s, out] : { std::pair(&a, &fa), std::pair(&b, &fb) })
    for (FaceId f : s->faces())
    {
      const auto ve = s->geometrium().vertices_edges(f);
      out->push_back({ s->geometrium().vertex(ve.v[0]).p, s->geometrium().vertex(ve.v[1]).p,
                       s->geometrium().vertex(ve.v[2]).p });
    }
  return std::ranges::equal(fa, fb, [](const auto& x, const auto& y) {
    return std::ranges::equal(x, y, [](Vec3 p, Vec3 q) { return p.x == q.x && p.y == q.y && p.z == q.z; });
  });
}

static void test_iso_seeded()
{
  const unsigned n = 40;
  const CartesianGrid cg{ n + 1, n + 1, n + 1, -1.25, 2.5 / n, -1.25, 2.5 / n, -1.25, 2.5 / n };
  unsigned evaluations = 0;
  auto sphere = [&](Vec3 p) { ++evaluations; return p.x * p.x + p.y * p.y + p.z * p.z; };

  // One part: the same faces as the full method, from a seed found by
  // bisection, with far fewer evaluations.
  Geometrium g1, g2;
  Surface& full = g1.new_surface();
  isosurface_cartesian(full, cg, sample(sphere), 1.);
  const unsigned full_evaluations = evaluations;
  evaluations = 0;
  Surface& seeded = g2.new_surface();
  const auto seeds = cells_crossing(cg, sphere, 1., { 0, 0, 0 }, { 1.2, 0.1, 0.05 });
  CHECK( ! seeds.empty());
  isosurface_cartesian(seeded, cg, sphere, 1., seeds);
  CHECK(same_faces(full, seeded) && full.vertex_number() == seeded.vertex_number() && full.edge_number() == seeded.edge_number());
  CHECK(is_closed(seeded) && is_orientable(seeded));
  CHECK(evaluations < full_evaluations / 4);

  // The surface of a previous step as the seed, with a band of one cell.
  Geometrium g3, g4;
  Surface& bigger = g3.new_surface();
  isosurface_cartesian(bigger, cg, sample(sphere), 1.1);
  Surface& tracked = g4.new_surface();
  isosurface_cartesian(tracked, cg, sphere, 1.1, cells_near(cg, seeded, 1));
  CHECK(same_faces(bigger, tracked));

  // Seeds the surface does not cross, or none.
  Geometrium g5;
  Surface& empty = g5.new_surface();
  CHECK(cells_crossing(cg, sphere, 1., { 0, 0, 0 }, { 0.5, 0, 0 }).empty());
  const std::array<Vec3, 1> centre{};
  isosurface_cartesian(empty, cg, sphere, 1., cells_near(cg, centre));
  isosurface_cartesian(empty, cg, sphere, 1., {});
  CHECK(empty.face_number() == 0 && empty.vertex_number() == 0);
  CHECK(cells_near(cg, std::array{ Vec3{ 5, 0, 0 } }).empty());
  CHECK(cells_near(cg, std::array{ Vec3{ 1.3, 1.25, -1.25 } }, 1).size() == 4);
  // A point on the far planes of the grid is in the last cell.
  const auto far = cells_near(cg, std::array{ Vec3{ 1.25, 1.25, 1.25 } });
  CHECK(far.size() == 1 && far[0] == (GridCell{ n - 1, n - 1, n - 1 }));
  {
    Geometrium gc;
    const CartesianGrid huge{ 1u << 31, 1u << 31, 8, 0, 1, 0, 1, 0, 1 };
    bool threw = false;
    try { isosurface_cartesian(gc.new_surface(), huge, sphere, 1., std::array{ GridCell{} }); }
    catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);
  }

  // Two parts: a seed on one makes that part only; seeds on both make what
  // the full method makes.
  auto blobs = [](Vec3 p) {
    auto b = [](Vec3 d) { return 1. / (d.x * d.x + d.y * d.y + d.z * d.z + 1e-12); };
    return b(p - Vec3{ 0.6, 0, 0 }) + b(p - Vec3{ -0.6, 0.1, 0 });
  };
  Geometrium g6, g7, g8;
  Surface& both_full = g6.new_surface();
  isosurface_cartesian(both_full, cg, sample(blobs), 16.);
  Surface& one = g7.new_surface();
  const auto right = cells_crossing(cg, blobs, 16., { 0.6, 0, 0 }, { 1.2, 0, 0 });
  isosurface_cartesian(one, cg, blobs, 16., right);
  CHECK(is_closed(one) && one.face_number() > 0 && one.face_number() < both_full.face_number());
  CHECK(std::ranges::all_of(one.vertices(), [&](VertexId v) { return g7.vertex(v).p.x > 0.; }));
  auto seeds2 = cells_crossing(cg, blobs, 16., { -0.6, 0.1, 0 }, { -1.2, 0.1, 0 });
  seeds2.insert(seeds2.end(), right.begin(), right.end());
  Surface& both = g8.new_surface();
  isosurface_cartesian(both, cg, blobs, 16., seeds2);
  CHECK(same_faces(both_full, both));
}

static void test_face_events()
{
  Geometrium g;
  Surface& s = g.new_surface();
  Surface& other = g.new_surface();
  std::vector<std::pair<FaceId, bool>> ev;
  Subscription sub = s.on_face([&](const FaceEvent& e) {
    // A face that leaves is still linked.
    if ( ! e.added) CHECK( ! g.dead(e.face) && ! s.contains(e.face));
    ev.emplace_back(e.face, e.added);
  });
  auto count = [&](bool added) { return std::ranges::count(ev, added, &std::pair<FaceId, bool>::second); };

  const Tet t = tetrahedron(g, s);
  CHECK(ev.size() == 4 && count(true) == 4);
  s.add_face(t.f[0]);
  other.add_face(t.f[0]);
  CHECK(ev.size() == 4);
  s.remove_face(t.f[0]);
  s.remove_face(t.f[0]);
  CHECK(ev.size() == 5 && ev.back() == std::pair(t.f[0], false) && ! g.dead(t.f[0]));
  s.add_face(t.f[0]);
  ev.clear();

  // Splits add the pieces; the parent keeps its index and stays.
  g.split_face4(t.f[3], [](Geometrium& g, EdgeId e) { return g.add_vertex(midpoint(g, e)); });
  CHECK(count(true) > 0 && count(false) == 0);
  ev.clear();
  // A swap adds two faces and removes two.
  Geometrium g2;
  Surface& s2 = g2.new_surface();
  std::vector<std::pair<FaceId, bool>> ev2;
  const std::array<VertexId, 4> q = { g2.add_vertex({ 0, 0, 0 }), g2.add_vertex({ 1, 0, 0 }),
                                      g2.add_vertex({ 1, 1, 0 }), g2.add_vertex({ 0, 1, 0 }) };
  const EdgeId diagonal = g2.add_edge(q[0], q[2]);
  s2.add_face(g2.add_face(g2.add_edge(q[0], q[1]), g2.add_edge(q[1], q[2]), diagonal));
  s2.add_face(g2.add_face(diagonal, g2.add_edge(q[2], q[3]), g2.add_edge(q[3], q[0])));
  Subscription sub2 = s2.on_face([&](const FaceEvent& e) { ev2.emplace_back(e.face, e.added); });
  swap_edge(s2, diagonal);
  CHECK(ev2.size() == 4 && ev2[0].second && ev2[1].second && ! ev2[2].second && ! ev2[3].second);
  // A reset subscription hears nothing.
  sub2.reset();
  ev2.clear();
  s2.remove_face(*s2.faces().begin());
  CHECK(ev2.empty());
  // A collapse and a destroy remove the faces on the edge.
  const std::size_t before = s.face_number();
  const EdgeId e = *g.edges_of(*s.vertices().begin()).begin();
  std::size_t on_e = 0;
  for (FaceId f : g.faces_of(e)) on_e += s.contains(f);
  g.collapse_edge(e, midpoint(g, e));
  CHECK(count(false) == std::ptrdiff_t(on_e) && s.face_number() == before - on_e);
  ev.clear();
  const EdgeId e2 = *g.edges_of(*s.vertices().begin()).begin();
  g.destroy(e2);
  CHECK(count(false) > 0 && count(true) == 0);
  // Removing a surface removes each of its faces; the subscription may
  // outlive the surface.
  ev.clear();
  const std::size_t n = s.face_number();
  g.remove_surface(other);
  CHECK(ev.empty());
  g.remove_surface(s);
  CHECK(n > 0 && count(false) == std::ptrdiff_t(n));
  sub.reset();
}

// Whether every point not removed is in the list of exactly one face of s,
// whose projection holds it, and there are n of them.
static bool points_consistent(FacePoints& fp, const Surface& s, std::size_t n)
{
  const Geometrium& g = s.geometrium();
  std::vector<int> seen(fp.size());
  std::size_t total = 0;
  bool ok = true;
  for (FaceId f : s.faces())
    for (FacePoints::PointId i : fp.points(f))
    {
      ++total;
      if (++seen[i] > 1 || point_in_triangle(fp.point(i), g.triangle(f)) == Intersect::Out) ok = false;
    }
  return ok && total == n;
}

static void test_face_points()
{
  Geometrium g;
  Surface& s = g.new_surface();
  unsigned x = 4321;
  auto rnd = [&] { x = x * 1103515245u + 12345u; return double((x >> 8) & 0xffff) / 65536.; };
  std::vector<Vec3> pts;
  for (int i = 0; i < 2000; ++i) pts.push_back({ rnd(), rnd(), 0. });
  const FaceId t = enclosing_triangle(g, pts, 10.);
  s.add_face(t);
  FacePoints fp(s);
  for (const Vec3& p : pts) fp.add(p, t);
  CHECK(fp.size() == 2000 && points_consistent(fp, s, 2000));

  // Delaunay insertions of some of the points, as vertices.
  std::size_t n = 2000;
  std::vector<VertexId> vs;
  auto insert = [&](FacePoints::PointId i) {
    const FaceId f = locate(s, fp.point(i));
    fp.remove(i);
    --n;
    vs.push_back(g.add_vertex(fp.point(i)));
    delaunay_add_vertex_to_face(s, vs.back(), f);
  };
  for (FacePoints::PointId i = 0; i < 300; ++i) insert(i);
  CHECK(points_consistent(fp, s, n));

  // Removals of vertices, which add faces before removing the old ones, in
  // turn with insertions and then alone, with no reads in between.
  for (int k = 0; k < 100; ++k)
  {
    delaunay_remove_vertex(s, vs[k]);
    insert(FacePoints::PointId(300 + k));
  }
  CHECK(points_consistent(fp, s, n));
  for (int k = 100; k < 200; ++k) delaunay_remove_vertex(s, vs[k]);
  CHECK(points_consistent(fp, s, n));

  // Changes in place, by the core's primitives.
  auto inner_edge = [&] {
    for (FaceId f : s.faces())
      for (EdgeId e : g.face(f).e)
      {
        int k = 0;
        for (FaceId h : g.faces_of(e)) k += s.contains(h);
        if (k == 2) return e;
      }
    return EdgeId();
  };
  tessellate(s, [](Geometrium& g, EdgeId e) { return g.add_vertex(midpoint(g, e)); });
  CHECK(points_consistent(fp, s, n));
  swap_edge(s, inner_edge());
  CHECK(points_consistent(fp, s, n));
  const EdgeId c = inner_edge();
  g.collapse_edge(c, midpoint(g, c));
  CHECK(points_consistent(fp, s, n));

  // Renumbering, and faces that leave with no face to take their points.
  g.collect();
  CHECK(points_consistent(fp, s, n));
  const std::size_t faces = s.face_number();
  delaunay_remove_hull(s);
  CHECK(s.face_number() < faces);
  std::size_t left = 0;
  for (FaceId f : s.faces()) left += count(fp.points(f));
  CHECK(left < n && points_consistent(fp, s, left));
}

// Random sequences of changes to a Delaunay triangulation, with the points
// read at random moments; the store must stay consistent after each.
static void test_face_points_random()
{
  for (unsigned seed : { 1u, 2u, 3u, 4u, 5u })
  {
    unsigned x = seed * 2654435761u;
    auto rnd = [&] { x = x * 1103515245u + 12345u; return double((x >> 8) & 0xffff) / 65536.; };
    auto pick = [&](std::size_t n) { return std::size_t(rnd() * double(n)); };
    Geometrium g;
    Surface& s = g.new_surface();
    std::vector<Vec3> pts;
    for (int i = 0; i < 1500; ++i) pts.push_back({ rnd(), rnd(), 0. });
    const FaceId t = enclosing_triangle(g, pts, 10.);
    const auto hull = g.vertices_of(t);
    s.add_face(t);
    FacePoints fp(s);
    for (const Vec3& p : pts) fp.add(p, t);

    std::size_t n = pts.size();
    FacePoints::PointId next = 0;
    std::vector<VertexId> inner;
    auto interior = [&](VertexId v) { return v != hull[0] && v != hull[1] && v != hull[2]; };
    // An edge of s between two faces whose quadrilateral is convex, so that
    // a swap keeps the faces counterclockwise.
    auto swappable = [&]() -> EdgeId {
      for (int tries = 0; tries < 50; ++tries)
      {
        const FaceId f = *std::next(s.faces().begin(), std::ptrdiff_t(pick(s.face_number())));
        const EdgeId e = g.face(f).e[pick(3)];
        FaceId other;
        for (FaceId h : g.faces_of(e)) if (h != f && s.contains(h)) other = h;
        if ( ! other.valid()) continue;
        const VertexId a = g.edge(e).v[0], b = g.edge(e).v[1];
        const VertexId c = g.opposite_vertex(f, e), d = g.opposite_vertex(other, e);
        const Vec3 pa = g.vertex(a).p, pb = g.vertex(b).p, pc = g.vertex(c).p, pd = g.vertex(d).p;
        if (orient2d(pc, pd, pa) * orient2d(pc, pd, pb) < 0. && orient2d(pa, pb, pc) * orient2d(pa, pb, pd) < 0.)
          return e;
      }
      return {};
    };
    // An inner edge whose collapse to its middle keeps every face around it
    // counterclockwise.
    auto collapsible_edge = [&]() -> EdgeId {
      for (int tries = 0; tries < 50; ++tries)
      {
        const FaceId f = *std::next(s.faces().begin(), std::ptrdiff_t(pick(s.face_number())));
        const EdgeId e = g.face(f).e[pick(3)];
        const VertexId a = g.edge(e).v[0], b = g.edge(e).v[1];
        if ( ! interior(a) || ! interior(b) || ! g.collapsible(e)) continue;
        const Vec3 m = midpoint(g, e);
        bool ok = true;
        for (VertexId w : { a, b })
          for (FaceId h : vertex_faces(s, w))
          {
            const auto vs = g.vertices_of(h);
            if (std::ranges::count(vs, a) + std::ranges::count(vs, b) == 2) continue;
            std::array<Vec3, 3> q;
            for (int k = 0; k < 3; ++k) q[k] = vs[k] == w ? m : g.vertex(vs[k]).p;
            if ( ! (orient2d(q[0], q[1], q[2]) > 0.)) ok = false;
          }
        if (ok) return e;
      }
      return {};
    };

    for (int step = 0; step < 400; ++step)
    {
      const double r = rnd();
      if (r < 0.45 && next < pts.size())
      {
        const FacePoints::PointId i = next++;
        fp.remove(i);
        --n;
        const VertexId v = g.add_vertex(pts[i]);
        if (delaunay_add_vertex(s, v).valid()) CHECK(false);
        inner.push_back(v);
      }
      else if (r < 0.65 && ! inner.empty())
      {
        const std::size_t k = pick(inner.size());
        delaunay_remove_vertex(s, inner[k]);
        inner.erase(inner.begin() + std::ptrdiff_t(k));
      }
      else if (r < 0.8)
      {
        if (const EdgeId e = swappable(); e.valid()) swap_edge(s, e);
      }
      else if (r < 0.9)
      {
        if (const EdgeId e = collapsible_edge(); e.valid())
        {
          const VertexId a = g.edge(e).v[0], b = g.edge(e).v[1];
          const VertexId v = g.collapse_edge(e, midpoint(g, e));
          std::erase_if(inner, [&](VertexId w) { return w == a || w == b; });
          if (v.valid()) inner.push_back(v);
        }
      }
      else if (r < 0.93)
      {
        // Vertex ids change; keep track of the interior ones through the map.
        std::vector<VertexId> renumbered;
        Subscription sub = g.on_renumber([&](const Renumbering& rn) {
          for (VertexId v : inner) if (const VertexId w = rn(v); w.valid()) renumbered.push_back(w);
        });
        g.collect();
        inner = renumbered;
      }
      if (rnd() < 0.2 && ! points_consistent(fp, s, n))
      {
        CHECK(false);
        std::println("  seed {} step {}", seed, step);
        break;
      }
    }
    CHECK(points_consistent(fp, s, n));
  }
}

// Cases from the reviews of 2026-10-13, each of which failed before.
static void test_robustness()
{
  // An edge at a contact vertex, where two sets of faces meet only at that
  // vertex, is not collapsible; a progressive surface keeps the vertex and
  // expands back to the input.
  {
    Geometrium g;
    Surface& s = g.new_surface();
    std::array<VertexId, 6> v;
    for (int i = 0; i < 6; ++i) v[i] = g.add_vertex({ double(i % 3), double(i / 3), double(i % 2) });
    auto edge = [&](int i, int j) { const EdgeId e = g.find_edge(v[i], v[j]); return e.valid() ? e : g.add_edge(v[i], v[j]); };
    for (const char* t : { "514", "124", "032" })
    {
      const int i = t[0] - '0', j = t[1] - '0', k = t[2] - '0';
      s.add_face(g.add_face(edge(i, j), edge(j, k), edge(k, i)));
    }
    CHECK( ! g.collapsible(g.find_edge(v[2], v[4])));
    const std::size_t nf = s.face_number(), nv = s.vertex_number();
    PSurface ps(s, stop_below_edges(0));
    bool ok = true;
    try { ps.set_vertex_number(ps.max_vertex_number()); }
    catch (const std::exception&) { ok = false; }
    CHECK(ok && s.face_number() == nf && s.vertex_number() == nv);
  }
  // A split whose flags do not fit the surface is a ReadError, and g is
  // left as it was.
  {
    Geometrium g0;
    Surface& s0 = g0.new_surface();
    generate_sphere(s0, 2);
    PSurface ps(s0, stop_below_edges(0));
    std::ostringstream out;
    ps.write(out);
    std::vector<std::string> lines;
    {
      std::istringstream in(out.str());
      for (std::string l; std::getline(in, l); ) lines.push_back(l);
    }
    unsigned nv, ne, nf;
    std::istringstream(lines[0]) >> nv >> ne >> nf;
    const std::size_t h = 1 + nv + ne + nf;
    const std::string t = lines[h + 4].substr(0, lines[h + 4].find(' '));
    Geometrium g;
    Surface& other = g.new_surface();
    tetrahedron(g, other);
    int failed = 0;
    bool clean = true;
    for (unsigned flags = 0; flags < 64; ++flags)
    {
      std::vector<std::string> bad = lines;
      bad[h + 4] = std::format("{} {}", t, flags);
      std::string text;
      for (const auto& l : bad) text += l + '\n';
      std::istringstream in(text);
      const std::size_t cv = g.vertex_capacity(), ce = g.edge_capacity(), cf = g.face_capacity();
      const std::size_t ns = g.surface_number();
      try
      {
        if (PSurface::read(g, in)) continue;
        ++failed;
        if (g.vertex_capacity() != cv || g.edge_capacity() != ce || g.face_capacity() != cf || g.surface_number() != ns)
          clean = false;
      }
      catch (const std::exception&) { clean = false; }
    }
    CHECK(failed > 0 && clean && other.face_number() == 4);
  }
  // A callout that cancels another subscription during a call: the
  // cancelled one is not called, and its state may be gone.
  {
    Geometrium g;
    Surface& s = g.new_surface();
    int late = 0;
    std::unique_ptr<Subscription> b;
    Subscription a = s.on_face([&](const FaceEvent&) { b.reset(); });
    b = std::make_unique<Subscription>(s.on_face([&](const FaceEvent&) { if ( ! b) ++late; }));
    tetrahedron(g, s);
    CHECK( ! b && late == 0);
  }
  // import() into a surface of another Geometrium.
  {
    Geometrium a, b;
    Surface& sa = a.new_surface();
    tetrahedron(a, sa);
    Surface& sb = b.new_surface();
    bool threw = false;
    try { a.import(a, sa, sb); }
    catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw && sb.face_number() == 0);
  }
  // A SplitMemo whose half was split again.
  {
    Geometrium g;
    Surface& s = g.new_surface();
    const std::array<VertexId, 4> q = { g.add_vertex({ 0, 0, 0 }), g.add_vertex({ 1, 0, 0 }),
                                        g.add_vertex({ 1, 1, 0 }), g.add_vertex({ 0, 1, 0 }) };
    const EdgeId diagonal = g.add_edge(q[0], q[2]);
    const FaceId f1 = g.add_face(g.add_edge(q[0], q[1]), g.add_edge(q[1], q[2]), diagonal);
    s.add_face(f1);
    s.add_face(g.add_face(diagonal, g.add_edge(q[2], q[3]), g.add_edge(q[3], q[0])));
    SplitMemo memo;
    const Splitter mid = [](Geometrium& g, EdgeId e) { return g.add_vertex(midpoint(g, e)); };
    g.split_face4(f1, mid, memo);
    EdgeId half;
    for (EdgeId e : g.edges_of(q[0]))
    {
      const VertexId o = g.edge(e).v[0] == q[0] ? g.edge(e).v[1] : g.edge(e).v[0];
      if (g.vertex(o).p == Vec3{ 0.5, 0.5, 0 }) half = e;
    }
    CHECK(half.valid());
    g.split_edge(half, mid(g, half));
    bool threw = false;
    try { g.split_remaining(memo); }
    catch (const std::logic_error&) { threw = true; }
    CHECK(threw);
  }
  // collapse_edge() of a lone edge leaves no vertex; of an edge with a
  // duplicate, it throws before changing anything.
  {
    Geometrium g;
    const VertexId a = g.add_vertex({ 0, 0, 0 }), b = g.add_vertex({ 1, 0, 0 });
    const EdgeId e = g.add_edge(a, b);
    CHECK( ! g.collapse_edge(e, { 0.5, 0, 0 }).valid());
    std::size_t live = 0;
    for (Index i = 0; i < g.vertex_capacity(); ++i) live += ! g.dead(VertexId(i));
    CHECK(live == 0);
    const VertexId c = g.add_vertex({ 0, 0, 0 }), d = g.add_vertex({ 1, 0, 0 });
    const EdgeId e1 = g.add_edge(c, d);
    g.add_edge(c, d);
    bool threw = false;
    try { g.collapse_edge(e1, { 0.5, 0, 0 }); }
    catch (const std::logic_error&) { threw = true; }
    CHECK(threw && ! g.dead(e1) && ! g.dead(c) && ! g.dead(d));
  }
  // split_components() takes every face once, also a non-manifold face
  // that leaves a component after the walk passed it.
  {
    Geometrium g;
    Surface& s = g.new_surface();
    std::array<VertexId, 7> v;
    for (int i = 0; i < 7; ++i) v[i] = g.add_vertex({ double(i), double(i * i % 5), double(i % 3) });
    auto edge = [&](int i, int j) { const EdgeId e = g.find_edge(v[i], v[j]); return e.valid() ? e : g.add_edge(v[i], v[j]); };
    for (const char* t : { "041", "254", "061", "063", "025", "206", "450", "413", "430", "536", "256" })
    {
      const int i = t[0] - '0', j = t[1] - '0', k = t[2] - '0';
      s.add_face(g.add_face(edge(i, j), edge(j, k), edge(k, i)));
    }
    std::map<FaceId, int> seen;
    for (Surface* c : split_components(s))
      for (FaceId f : c->faces()) ++seen[f];
    CHECK(seen.size() == 11 && std::ranges::all_of(seen, [](const auto& p) { return p.second == 1; }));
  }
  // Readers: a Jostle file whose first node has no neighbours, and counts
  // that are negative, huge, or more than the file holds.
  {
    std::istringstream jostle("3 1\n\n3\n2\n");
    Graph gr;
    const auto r = read_jostle(gr, jostle);
    CHECK(r && r->size() == 3 && edges_of(*r).size() == 1);
    for (const char* text : { "2000000000 0\n", "-1 0\n", "3 0\n\n\n" })
    {
      std::istringstream in(text);
      Graph gr2;
      CHECK( ! read_graph(gr2, in).has_value());
    }
    for (const char* text : { "-5 6 4\n", "4 -6 4\n", "999999999999 6 4\n0 0 0\n" })
    {
      Geometrium g;
      std::istringstream in(text);
      CHECK( ! read_gts(g, in).has_value());
    }
  }
  // Progressive surfaces read back: a face of the split being read, a file
  // cut anywhere after the surface, and a split without faces.
  {
    Geometrium g;
    Surface& s = g.new_surface();
    generate_sphere(s, 2);
    PSurface ps(s, stop_below_edges(0));
    std::ostringstream out;
    ps.write(out);
    std::vector<std::string> lines;
    {
      std::istringstream in(out.str());
      for (std::string l; std::getline(in, l); ) lines.push_back(l);
    }
    unsigned nv, ne, nf;
    std::istringstream(lines[0]) >> nv >> ne >> nf;
    const std::size_t h = 1 + nv + ne + nf; // the number of splits
    auto read = [&](const std::vector<std::string>& ls) {
      std::string text;
      for (const auto& l : ls) text += l + '\n';
      Geometrium g2;
      std::istringstream in(text);
      return PSurface::read(g2, in).has_value();
    };
    CHECK(read(lines));
    for (std::size_t cut = h + 1; cut < lines.size(); ++cut)
      if (read(std::vector(lines.begin(), lines.begin() + std::ptrdiff_t(cut)))) { CHECK(false); break; }
    std::vector<std::string> bad = lines;
    std::istringstream(lines[h + 4]) >> ne; // the first face, and its flags
    unsigned flags;
    std::istringstream(lines[h + 4]).ignore(64, ' ') >> flags;
    bad[h + 7] = std::format("{} {}", nf + 1, flags); // the split's own first face
    CHECK( ! read(bad));
    bad = lines;
    bad[h + 1] = std::format("{} 0", lines[h + 1].substr(0, lines[h + 1].find(' ')));
    CHECK( ! read(bad));
  }
  // An invalid index stays invalid through renumbering; a dead face cannot
  // join a surface; expand() refuses a split that is not collapsed.
  {
    Geometrium g;
    Surface& s = g.new_surface();
    const Tet t = tetrahedron(g, s);
    Renumbering r;
    r.faces.assign(g.face_capacity(), 0);
    CHECK( ! r(FaceId()).valid() && ! r(VertexId()).valid() && ! r(EdgeId()).valid());
    s.remove_face(t.f[0]);
    bool threw = false;
    try { s.add_face(t.f[0]); }
    catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);
  }
  {
    Geometrium g;
    Surface& s = g.new_surface();
    generate_sphere(s, 2);
    PSurface ps(s, stop_below_edges(0));
    ps.set_vertex_number(ps.max_vertex_number());
    VertexSplit vs = ps.splits().front();
    bool threw = false;
    try { expand(g, vs, s); }
    catch (const std::logic_error&) { threw = true; }
    CHECK(threw);
  }
}

static void test_height_field()
{
  // A plane: the corners make it exactly, so nothing else is inserted.
  {
    Geometrium g;
    Surface& s = g.new_surface();
    const std::size_t w = 30, h = 20;
    std::vector<double> z(w * h);
    for (std::size_t j = 0; j < h; ++j)
      for (std::size_t i = 0; i < w; ++i) z[j * w + i] = 2. * i + 3. * j + 1.;
    approximate_grid(s, w, h, z, stop_below_cost(1e-9));
    CHECK(s.vertex_number() == 4 && s.face_number() == 2);
  }
  // A terrain: stopping below a cost leaves every grid point within it.
  {
    Geometrium g;
    Surface& s = g.new_surface();
    const std::size_t w = 60, h = 40;
    std::vector<double> z(w * h);
    for (std::size_t j = 0; j < h; ++j)
      for (std::size_t i = 0; i < w; ++i) z[j * w + i] = 10. * std::sin(0.2 * i) * std::cos(0.15 * j) + 0.1 * i;
    approximate_grid(s, w, h, z, stop_below_cost(0.5));
    CHECK(s.vertex_number() > 4 && s.vertex_number() < w * h / 4);
    double worst = 0.;
    for (std::size_t j = 0; j < h; ++j)
      for (std::size_t i = 0; i < w; ++i)
      {
        const Vec3 p{ double(i), double(j), z[j * w + i] };
        const FaceId f = locate(s, p);
        if (f.valid()) worst = std::max(worst, std::abs(interpolate_height(g, f, p) - p.z));
      }
    CHECK(worst < 0.5 + 1e-9);
    CHECK(is_orientable(s) && boundary(s).size() > 0);

    Geometrium g2;
    Surface& s2 = g2.new_surface();
    approximate_grid(s2, w, h, z, stop_at_vertices(100));
    CHECK(s2.vertex_number() == 100);
    bool threw = false;
    try { approximate_grid(s2, w, h, std::span(z).first(10), stop_at_vertices(100)); }
    catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);
  }
  // Scattered points, closed by a face on the enclosing triangle.
  {
    Geometrium g;
    Surface& s = g.new_surface();
    unsigned x = 777;
    auto rnd = [&] { x = x * 1103515245u + 12345u; return double((x >> 8) & 0xffff) / 65536.; };
    std::vector<Vec3> pts;
    for (int i = 0; i < 500; ++i) { const double px = rnd(), py = rnd(); pts.push_back({ px, py, std::sin(5 * px) * py }); }
    approximate_points(s, pts, stop_at_vertices(60), height_cost, Enclosing::Close);
    CHECK(is_closed(s) && is_orientable(s) && s.vertex_number() == 59);
    bool threw = false;
    try { approximate_points(s, {}, stop_at_vertices(60)); }
    catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);
  }
}

int main()
{
  std::println("sizeof Vertex {}, Edge {}, Face {}", sizeof(Vertex), sizeof(Edge), sizeof(Face));
  test_build_and_rings();
  test_orientation();
  test_release_cascade();
  test_keep_orphans();
  test_shared_faces();
  test_non_manifold_and_duplicates();
  test_ring_order();
  test_attributes();
  test_import();
  test_io();
  test_threads();
  test_set_face_edges();
  test_eheap();
  test_default_hooks();
  test_split_edge();
  test_sphere();
  test_tessellate();
  test_refine();
  test_collapse();
  test_swap_edge();
  test_coarsen();
  test_iso();
  test_predicates_and_bbtree();
  test_boolean();
  test_small_functions();
  test_cluster();
  test_distance_inside_curvature_merge();
  test_fixes();
  test_topology();
  test_delaunay();
  test_matrix();
  test_boxes_and_writers();
  test_strips();
  test_park();
  test_progressive();
  test_graph();
  test_isoslice();
  test_iso_seeded();
  test_face_events();
  test_face_points();
  test_face_points_random();
  test_height_field();
  test_robustness();
  std::println("{}", g_failed ? "FAILED" : "all passed");
  return g_failed != 0;
}
