// gts-ref: runs one GTS operation through the C library (Fedora's libgts)
// and prints a summary of the resulting surface that does not depend on the
// order of vertices, edges or faces. Optionally writes the surface.
//
// The operations and their parameters mirror how GledGTS calls GTS.
//
// Copyright (C) 2026 Matevz Tadel
// SPDX-License-Identifier: LGPL-2.1-or-later

#include <gts.h>

#include <map>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <format>
#include <numbers>
#include <print>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace
{

//==============================================================================
// Errors and glib messages
//==============================================================================

[[noreturn]] void fail(std::string_view msg)
{
  std::println("status error {}", msg);
  std::exit(1);
}

// GTS reports failed argument checks as glib CRITICAL messages and goes on.
// They are counted and printed, because a C++ version that asserts instead
// would otherwise look like a different result.
int g_n_critical = 0;
int g_n_warning  = 0;

void log_counter(const gchar* domain, GLogLevelFlags level, const gchar* msg, gpointer)
{
  if (level & G_LOG_LEVEL_CRITICAL) ++g_n_critical;
  if (level & G_LOG_LEVEL_WARNING)  ++g_n_warning;
  std::fprintf(stderr, "glib %s: %s\n", domain ? domain : "-", msg);
}

//==============================================================================
// Surface I/O
//==============================================================================

GtsSurface* new_surface()
{
  return gts_surface_new(gts_surface_class(), gts_face_class(),
                         gts_edge_class(), gts_vertex_class());
}

GtsSurface* read_surface(const char* path)
{
  FILE* fp = std::fopen(path, "r");
  if ( ! fp) fail(std::string("cannot open ") + path);
  GtsSurface* s  = new_surface();
  GtsFile*    gf = gts_file_new(fp);
  if (gts_surface_read(s, gf) != 0)
    fail(std::string("cannot read ") + path + ": " + (gf->error ? gf->error : "?"));
  gts_file_destroy(gf);
  std::fclose(fp);
  return s;
}

void write_surface(GtsSurface* s, const char* path)
{
  FILE* fp = std::fopen(path, "w");
  if ( ! fp) fail(std::string("cannot write ") + path);
  gts_surface_write(s, fp);
  std::fclose(fp);
}

//==============================================================================
// Summary
//==============================================================================

void print_range(const char* name, const GtsRange& r)
{
  std::println("{} n={} min={:.12e} mean={:.12e} max={:.12e}", name, r.n, r.min, r.mean, r.max);
}

void bbox_vertex(GtsVertex* v, double* b)
{
  const GtsPoint& p = v->p;
  b[0] = std::fmin(b[0], p.x); b[3] = std::fmax(b[3], p.x);
  b[1] = std::fmin(b[1], p.y); b[4] = std::fmax(b[4], p.y);
  b[2] = std::fmin(b[2], p.z); b[5] = std::fmax(b[5], p.z);
}

void print_summary(GtsSurface* s)
{
  std::println("vertices {}", gts_surface_vertex_number(s));
  std::println("edges {}",    gts_surface_edge_number(s));
  std::println("faces {}",    gts_surface_face_number(s));
  if (gts_surface_face_number(s) == 0)
    return;

  std::println("area {:.12e}",   gts_surface_area(s));
  std::println("volume {:.12e}", gts_surface_volume(s));
  std::println("manifold {}",    gts_surface_is_manifold(s)   ? 1 : 0);
  std::println("closed {}",      gts_surface_is_closed(s)     ? 1 : 0);
  std::println("orientable {}",  gts_surface_is_orientable(s) ? 1 : 0);

  GtsSurface* si = gts_surface_is_self_intersecting(s);
  std::println("self_intersecting {}", si ? 1 : 0);
  if (si) gts_object_destroy(GTS_OBJECT(si));

  double b[6] = { HUGE_VAL, HUGE_VAL, HUGE_VAL, -HUGE_VAL, -HUGE_VAL, -HUGE_VAL };
  gts_surface_foreach_vertex(s, (GtsFunc) bbox_vertex, b);
  std::println("bbox {:.12e} {:.12e} {:.12e} {:.12e} {:.12e} {:.12e}", b[0], b[1], b[2], b[3], b[4], b[5]);

  GtsSurfaceStats st;
  gts_surface_stats(s, &st);
  std::println("incompatible_faces {}", st.n_incompatible_faces);
  std::println("duplicate_faces {}",    st.n_duplicate_faces);
  std::println("duplicate_edges {}",    st.n_duplicate_edges);
  std::println("boundary_edges {}",     st.n_boundary_edges);
  std::println("non_manifold_edges {}", st.n_non_manifold_edges);
  print_range("edges_per_vertex", st.edges_per_vertex);
  print_range("faces_per_edge",   st.faces_per_edge);

  GtsSurfaceQualityStats qs;
  gts_surface_quality_stats(s, &qs);
  print_range("face_quality", qs.face_quality);
  print_range("face_area",    qs.face_area);
  print_range("edge_length",  qs.edge_length);
  print_range("edge_angle",   qs.edge_angle);
}

//==============================================================================
// Coarsening and refinement, as in GledGTS's GTSRetriangulator
//==============================================================================

// Cost from the angle opposite to the edge, from GTSRetriangulator.cxx:
// 1 for 60 degrees, falling linearly to 0 at 0 and 180 degrees. The cost of
// an edge is the smallest over its triangles.
gdouble edge_angle(GtsTriangle* t, GtsEdge* e)
{
  static const gdouble r60 = std::numbers::pi / 3;
  GtsVertex* v = gts_triangle_vertex_opposite(t, e);
  if ( ! v) return 1.0;
  GtsVector a, b, c;
  gts_vector_init(a, &v->p, &e->segment.v1->p);
  gts_vector_init(b, &v->p, &e->segment.v2->p);
  gts_vector_cross(c, a, b);
  gdouble phi = std::atan2(gts_vector_norm(c), gts_vector_scalar(a, b));
  return phi <= r60 ? phi / r60 : 1.5 - 0.5 * phi / r60;
}

gdouble cost_angle(GtsEdge* e)
{
  gdouble cost = 1.0;
  for (GSList* i = e->triangles; i; i = i->next)
    cost = std::fmin(edge_angle((GtsTriangle*) i->data, e), cost);
  return cost;
}

gboolean refine_stop_number(gdouble, guint number, guint* max)  { return number > *max; }
gboolean refine_stop_cost  (gdouble cost, guint, gdouble* min)  { return cost < *min; }

struct RetriangulateArgs
{
  enum Cost_e { C_Length, C_Volume, C_Angle } cost = C_Length;
  enum Stop_e { S_Number, S_Cost } stop = S_Number;
  guint   stop_number = 1000;
  gdouble stop_cost   = 0.5;
  bool    midvert_volume = true;
  gdouble min_angle_deg  = 1;
  GtsVolumeOptimizedParams vo = { 0.5, 0.5, 0 };
};

RetriangulateArgs parse_retriangulate(std::span<char*> a, bool for_coarsen)
{
  // COST STOP VALUE [MIDVERT [MINANGLE_DEG]]
  if (a.size() < 3) fail("expected COST STOP VALUE");
  RetriangulateArgs r;
  std::string_view c = a[0], s = a[1];
  if      (c == "length") r.cost = RetriangulateArgs::C_Length;
  else if (c == "volume") r.cost = RetriangulateArgs::C_Volume;
  else if (c == "angle")  r.cost = RetriangulateArgs::C_Angle;
  else fail("COST must be length, volume or angle");
  if      (s == "number") { r.stop = RetriangulateArgs::S_Number; r.stop_number = std::strtoul(a[2], nullptr, 10); }
  else if (s == "cost")   { r.stop = RetriangulateArgs::S_Cost;   r.stop_cost   = std::strtod(a[2], nullptr); }
  else fail("STOP must be number or cost");
  if (for_coarsen && a.size() > 3)
  {
    std::string_view m = a[3];
    if      (m == "mid")    r.midvert_volume = false;
    else if (m == "volume") r.midvert_volume = true;
    else fail("MIDVERT must be mid or volume");
  }
  if (for_coarsen && a.size() > 4) r.min_angle_deg = std::strtod(a[4], nullptr);
  // The length cost works on the square of the edge length.
  if (r.cost == RetriangulateArgs::C_Length) r.stop_cost *= r.stop_cost;
  return r;
}

GtsKeyFunc cost_func(RetriangulateArgs& r, gpointer& data)
{
  data = nullptr;
  switch (r.cost)
  {
    case RetriangulateArgs::C_Length: return nullptr;
    case RetriangulateArgs::C_Volume: data = &r.vo; return (GtsKeyFunc) gts_volume_optimized_cost;
    case RetriangulateArgs::C_Angle:  return (GtsKeyFunc) cost_angle;
  }
  return nullptr;
}

void coarsen(GtsSurface* s, RetriangulateArgs r)
{
  GtsStopFunc stop_func; gpointer stop_data;
  if (r.stop == RetriangulateArgs::S_Number) { stop_func = (GtsStopFunc) gts_coarsen_stop_number; stop_data = &r.stop_number; }
  else                                       { stop_func = (GtsStopFunc) gts_coarsen_stop_cost;   stop_data = &r.stop_cost; }
  gpointer   c_data;
  GtsKeyFunc c_func = cost_func(r, c_data);
  GtsCoarsenFunc m_func = r.midvert_volume ? (GtsCoarsenFunc) gts_volume_optimized_vertex : nullptr;
  gpointer       m_data = r.midvert_volume ? &r.vo : nullptr;
  gts_surface_coarsen(s, c_func, c_data, m_func, m_data, stop_func, stop_data,
                      r.min_angle_deg * std::numbers::pi / 180);
}

void refine(GtsSurface* s, RetriangulateArgs r)
{
  GtsStopFunc stop_func; gpointer stop_data;
  if (r.stop == RetriangulateArgs::S_Number) { stop_func = (GtsStopFunc) refine_stop_number; stop_data = &r.stop_number; }
  else                                       { stop_func = (GtsStopFunc) refine_stop_cost;   stop_data = &r.stop_cost; }
  gpointer   c_data;
  GtsKeyFunc c_func = cost_func(r, c_data);
  gts_surface_refine(s, c_func, c_data, nullptr, nullptr, stop_func, stop_data);
}

//==============================================================================
// Tessellation, as in GTSurf::Tessellate
//==============================================================================

GtsVertex* mid_edge_splitter(GtsEdge* e, GtsVertexClass* k, gpointer)
{
  GtsPoint &a = e->segment.v1->p, &b = e->segment.v2->p;
  return gts_vertex_new(k, a.x + 0.5 * (b.x - a.x), a.y + 0.5 * (b.y - a.y), a.z + 0.5 * (b.z - a.z));
}

//==============================================================================
// Iso-surfaces, as in GTSIsoMaker::MakeSurface
//==============================================================================

using field_t = double (*)(double x, double y, double z);

double field_sphere(double x, double y, double z) { return x*x + y*y + z*z; }

double field_torus(double x, double y, double z)
{
  // Major radius 0.7, minor radius 0.3; the surface is at 0.09.
  const double q = std::sqrt(x*x + y*y) - 0.7;
  return q*q + z*z;
}

double field_blobs(double x, double y, double z)
{
  // Two blobs, 1/r^2 potentials summed; the surface is at 4.
  auto p = [](double dx, double dy, double dz) { return 1.0 / (dx*dx + dy*dy + dz*dz + 1e-12); };
  return p(x - 0.4, y, z) + p(x + 0.4, y - 0.1, z + 0.05);
}

void field_to_plane(gdouble** a, GtsCartesianGrid g, guint, gpointer data)
{
  field_t f = (field_t) data;
  double x = g.x;
  for (guint i = 0; i < g.nx; ++i, x += g.dx)
  {
    double y = g.y;
    for (guint j = 0; j < g.ny; ++j, y += g.dy)
      a[i][j] = f(x, y, g.z);
  }
}

GtsSurface* make_iso(std::string_view algo, std::string_view func, guint n)
{
  field_t f; double iso;
  if      (func == "sphere") { f = field_sphere; iso = 1;    }
  else if (func == "torus")  { f = field_torus;  iso = 0.09; }
  else if (func == "blobs")  { f = field_blobs;  iso = 4;    }
  else fail("FUNC must be sphere, torus or blobs");

  const double lo = -1.25, d = 2.5 / n;
  GtsCartesianGrid g = { n + 1, n + 1, n + 1, lo, d, lo, d, lo, d };

  GtsSurface* s = new_surface();
  if      (algo == "cartesian")     gts_isosurface_cartesian    (s, g, field_to_plane, (gpointer) f, iso);
  else if (algo == "tetra")         gts_isosurface_tetra        (s, g, field_to_plane, (gpointer) f, iso);
  else if (algo == "tetra_bounded") gts_isosurface_tetra_bounded(s, g, field_to_plane, (gpointer) f, iso);
  else if (algo == "tetra_bcl")     gts_isosurface_tetra_bcl    (s, g, field_to_plane, (gpointer) f, iso);
  else fail("ALGO must be cartesian, tetra, tetra_bounded or tetra_bcl");
  return s;
}

//==============================================================================
// Booleans, as in GTSBoolOpHelper with post-processing off
//==============================================================================

GtsSurface* make_boolean(std::string_view op, GtsSurface* a, GtsSurface* b)
{
  if (op != "union" && op != "inter" && op != "diff")
    fail("OP must be union, inter or diff");

  GNode* ta = gts_bb_tree_surface(a);
  GNode* tb = gts_bb_tree_surface(b);
  const gboolean a_open = gts_surface_volume(a) < 0;
  const gboolean b_open = gts_surface_volume(b) < 0;
  GtsSurfaceInter* si = gts_surface_inter_new(gts_surface_inter_class(), a, b, ta, tb, a_open, b_open);
  gts_bb_tree_destroy(ta, TRUE);
  gts_bb_tree_destroy(tb, TRUE);

  gboolean closed = FALSE;
  const gboolean orientable = gts_surface_inter_check(si, &closed);
  std::println("inter_edges {}",      g_slist_length(si->edges));
  std::println("inter_orientable {}", orientable ? 1 : 0);
  std::println("inter_closed {}",     closed     ? 1 : 0);
  if ( ! orientable || ! closed)
    fail("intersection curve not orientable or not closed");

  GtsSurface* r = new_surface();
  if (op == "union")
  {
    gts_surface_inter_boolean(si, r, GTS_1_OUT_2);
    gts_surface_inter_boolean(si, r, GTS_2_OUT_1);
  }
  else if (op == "inter")
  {
    gts_surface_inter_boolean(si, r, GTS_1_IN_2);
    gts_surface_inter_boolean(si, r, GTS_2_IN_1);
  }
  else
  {
    gts_surface_inter_boolean(si, r, GTS_1_OUT_2);
    gts_surface_inter_boolean(si, r, GTS_2_IN_1);
    gts_surface_foreach_face(si->s2, (GtsFunc) gts_triangle_revert, nullptr);
    gts_surface_foreach_face(b,      (GtsFunc) gts_triangle_revert, nullptr);
  }
  gts_object_destroy(GTS_OBJECT(si));
  return r;
}

//==============================================================================
// Cluster grid (out-of-core simplification), as in
// GTSRetriangulator::OutOfCoreSimplification
//==============================================================================

void cluster_add_face(GtsFace* f, GtsClusterGrid* g)
{
  GtsVertex *v1, *v2, *v3;
  gts_triangle_vertices(&f->triangle, &v1, &v2, &v3);
  gts_cluster_grid_add_triangle(g, &v1->p, &v2->p, &v3->p, nullptr);
}

GtsSurface* make_cluster(GtsSurface* s, double delta)
{
  GtsSurface* d = new_surface();
  double b[6] = { HUGE_VAL, HUGE_VAL, HUGE_VAL, -HUGE_VAL, -HUGE_VAL, -HUGE_VAL };
  gts_surface_foreach_vertex(s, (GtsFunc) bbox_vertex, b);
  GtsBBox* bbox = gts_bbox_new(gts_bbox_class(), d, b[0], b[1], b[2], b[3], b[4], b[5]);
  GtsClusterGrid* g = gts_cluster_grid_new(gts_cluster_grid_class(), gts_cluster_class(), d, bbox, delta);
  gts_surface_foreach_face(s, (GtsFunc) cluster_add_face, g);
  GtsRange r = gts_cluster_grid_update(g);
  print_range("clusters", r);
  gts_object_destroy(GTS_OBJECT(g));
  gts_object_destroy(GTS_OBJECT(bbox));
  return d;
}

//==============================================================================

//==============================================================================
// Distance, curvature, inside, merging and centres
//==============================================================================

// A key for a point that is the same in gts-cxx, so that results can be
// compared point by point whatever the order of the vertices.
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

struct CurvatureArg { GtsSurface* s; Lines* lines; };

void curvature_vertex(GtsVertex* v, CurvatureArg* arg)
{
  GtsSurface* s = arg->s;
  GtsVector kh;
  gdouble kg;
  const GtsPoint& p = v->p;
  const std::string k = key("curv", p.x, p.y, p.z);
  const bool hm = gts_vertex_mean_curvature_normal(v, s, kh);
  const bool hg = gts_vertex_gaussian_curvature(v, s, &kg);
  if ( ! hm || ! hg)
  {
    (*arg->lines)[k].push_back(std::format("mean={} gauss={}", int(hm), int(hg)));
    return;
  }
  const double nkh = std::sqrt(kh[0]*kh[0] + kh[1]*kh[1] + kh[2]*kh[2]);
  gdouble k1, k2;
  gts_vertex_principal_curvatures(0.5*nkh, kg, &k1, &k2);
  GtsVector e1, e2;
  gts_vertex_principal_directions(v, s, kh, kg, e1, e2);
  (*arg->lines)[k].push_back(std::format(
    "khx={:.12e} khy={:.12e} khz={:.12e} kg={:.12e} k1={:.12e} k2={:.12e} "
    "e1x={:.12e} e1y={:.12e} e1z={:.12e} e2x={:.12e} e2y={:.12e} e2z={:.12e}",
    kh[0], kh[1], kh[2], kg, k1, k2, e1[0], e1[1], e1[2], e2[0], e2[1], e2[2]));
}

GtsPoint* triangle_closest(GtsPoint* p, gpointer bounded)
{
  GtsPoint* c = gts_point_new(gts_point_class(), 0., 0., 0.);
  gts_point_triangle_closest(p, GTS_TRIANGLE(bounded), c);
  return c;
}

// Points on an n^3 grid over the box of s, enlarged by a tenth.
void inside_cmd(GtsSurface* s, unsigned n)
{
  double b[6] = { G_MAXDOUBLE, G_MAXDOUBLE, G_MAXDOUBLE, -G_MAXDOUBLE, -G_MAXDOUBLE, -G_MAXDOUBLE };
  gts_surface_foreach_vertex(s, (GtsFunc) bbox_vertex, b);
  GNode* tree = gts_bb_tree_surface(s);
  GtsPoint* p = gts_point_new(gts_point_class(), 0., 0., 0.);
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
        gts_point_set(p, q[0], q[1], q[2]);
        const bool in = gts_point_is_inside_surface(p, tree, FALSE);
        n_in += in;
        const double d = gts_bb_tree_point_distance(tree, p, (GtsBBoxDistFunc) gts_point_triangle_distance, nullptr);
        double d2;
        GtsPoint* c = gts_bb_tree_point_closest(tree, p, triangle_closest, &d2);
        // Closest points can tie, so only their distance is printed.
        std::println("{} in={} d={:.12e} d2={:.12e}", key("pt", q[0], q[1], q[2]), int(in), d, d2);
        gts_object_destroy(GTS_OBJECT(c));
      }
  std::println("inside n={} in={}", n*n*n, n_in);
  gts_object_destroy(GTS_OBJECT(p));
  gts_bb_tree_destroy(tree, TRUE);
}

void prepend(gpointer data, GSList** list)
{
  *list = g_slist_prepend(*list, data);
}

struct Unweld
{
  GtsSurface* s;
  GList*      vertices = nullptr;
  GList*      edges    = nullptr;
};

void unweld_face(GtsTriangle* t, Unweld* u)
{
  GtsVertex *a, *b, *c;
  gts_triangle_vertices(t, &a, &b, &c);
  GtsVertex* v[3];
  GtsVertex* o[3] = { a, b, c };
  for (int k = 0; k < 3; ++k)
  {
    v[k] = gts_vertex_new(u->s->vertex_class, o[k]->p.x, o[k]->p.y, o[k]->p.z);
    u->vertices = g_list_prepend(u->vertices, v[k]);
  }
  GtsEdge* e[3];
  for (int k = 0; k < 3; ++k)
  {
    e[k] = gts_edge_new(u->s->edge_class, v[k], v[(k + 1) % 3]);
    u->edges = g_list_prepend(u->edges, e[k]);
  }
  gts_surface_add_face(u->s, gts_face_new(u->s->face_class, e[0], e[1], e[2]));
}

// Every face of in on vertices of its own, then merged back together.
GtsSurface* merge_cmd(GtsSurface* in, double eps)
{
  Unweld u{ new_surface() };
  gts_surface_foreach_face(in, (GtsFunc) unweld_face, &u);
  u.vertices = g_list_reverse(u.vertices);
  u.edges    = g_list_reverse(u.edges);
  std::println("unwelded vertices={} edges={}", g_list_length(u.vertices), g_list_length(u.edges));
  u.vertices = gts_vertices_merge(u.vertices, eps, nullptr);
  u.edges    = gts_edges_merge(u.edges);
  std::println("merged vertices={} edges={}", g_list_length(u.vertices), g_list_length(u.edges));
  // triangle_cleanup of GTS's example cleanup.c
  GSList* faces = nullptr;
  gts_surface_foreach_face(u.s, (GtsFunc) prepend, &faces);
  unsigned n = 0;
  for (GSList* i = faces; i; i = i->next)
    if (gts_triangle_is_duplicate(GTS_TRIANGLE(i->data)))
    {
      gts_object_destroy(GTS_OBJECT(i->data));
      ++n;
    }
  g_slist_free(faces);
  std::println("duplicate_faces_removed {}", n);
  g_list_free(u.vertices);
  g_list_free(u.edges);
  return u.s;
}

void centers_cmd(GtsSurface* s)
{
  GtsVector cm, ca;
  const double vol  = gts_surface_center_of_mass(s, cm);
  const double area = gts_surface_center_of_area(s, ca);
  std::println("center_of_mass volume={:.12e} x={:.12e} y={:.12e} z={:.12e}", vol, cm[0], cm[1], cm[2]);
  std::println("center_of_area area={:.12e} x={:.12e} y={:.12e} z={:.12e}", area, ca[0], ca[1], ca[2]);
  GSList* bnd = gts_surface_boundary(s);
  std::println("boundary edges={}", g_slist_length(bnd));
  g_slist_free(bnd);
}

void usage()
{
  std::print(stderr,
    "usage: gts-ref [-o OUT.gts] COMMAND ARGS...\n"
    "  stats      IN\n"
    "  coarsen    IN COST STOP VALUE [MIDVERT [MINANGLE_DEG]]\n"
    "  refine     IN COST STOP VALUE\n"
    "               COST: length | volume | angle   STOP: number | cost\n"
    "               MIDVERT: mid | volume (default volume)\n"
    "  tessellate IN ORDER [mid]\n"
    "  sphere     ORDER\n"
    "  iso        ALGO FUNC N\n"
    "               ALGO: cartesian | tetra | tetra_bounded | tetra_bcl\n"
    "               FUNC: sphere | torus | blobs, on [-1.25, 1.25]^3 with N cells per axis\n"
    "  bool       OP A B           OP: union | inter | diff\n"
    "  cluster    IN DELTA\n"
    "  distance   A B [DELTA]      distance from A to B and from B to A\n"
    "  curvature  IN               curvature at each vertex\n"
    "  inside     IN N             inside, distance and closest point on an N^3 grid\n"
    "  merge      IN EPS           faces on vertices of their own, merged back\n"
    "  centers    IN               centres of mass and area, boundary\n");
  std::exit(2);
}

} // namespace

int main(int argc, char* argv[])
{
  g_log_set_handler(nullptr, GLogLevelFlags(G_LOG_LEVEL_MASK | G_LOG_FLAG_FATAL | G_LOG_FLAG_RECURSION),
                    log_counter, nullptr);
  g_log_set_handler("Gts", GLogLevelFlags(G_LOG_LEVEL_MASK | G_LOG_FLAG_FATAL | G_LOG_FLAG_RECURSION),
                    log_counter, nullptr);

  std::vector<char*> args(argv + 1, argv + argc);
  const char* out = nullptr;
  if (args.size() >= 2 && std::string_view(args[0]) == "-o")
  {
    out = args[1];
    args.erase(args.begin(), args.begin() + 2);
  }
  if (args.empty()) usage();

  const std::string_view cmd = args[0];
  std::span<char*> a(args.begin() + 1, args.end());

  std::print("command");
  for (char* x : args) std::print(" {}", x);
  std::println();

  GtsSurface* s = nullptr;

  if (cmd == "stats" && a.size() == 1)
  {
    s = read_surface(a[0]);
  }
  else if (cmd == "coarsen" && a.size() >= 4)
  {
    s = read_surface(a[0]);
    coarsen(s, parse_retriangulate(a.subspan(1), true));
  }
  else if (cmd == "refine" && a.size() == 4)
  {
    s = read_surface(a[0]);
    refine(s, parse_retriangulate(a.subspan(1), false));
  }
  else if (cmd == "tessellate" && (a.size() == 2 || a.size() == 3))
  {
    s = read_surface(a[0]);
    const bool mid = a.size() == 3 && std::string_view(a[2]) == "mid";
    for (unsigned long n = std::strtoul(a[1], nullptr, 10); n > 0; --n)
      gts_surface_tessellate(s, mid ? mid_edge_splitter : nullptr, nullptr);
  }
  else if (cmd == "sphere" && a.size() == 1)
  {
    s = new_surface();
    gts_surface_generate_sphere(s, std::strtoul(a[0], nullptr, 10));
  }
  else if (cmd == "iso" && a.size() == 3)
  {
    s = make_iso(a[0], a[1], std::strtoul(a[2], nullptr, 10));
  }
  else if (cmd == "bool" && a.size() == 3)
  {
    GtsSurface* sa = read_surface(a[1]);
    GtsSurface* sb = read_surface(a[2]);
    s = make_boolean(a[0], sa, sb);
    gts_object_destroy(GTS_OBJECT(sa));
    gts_object_destroy(GTS_OBJECT(sb));
  }
  else if (cmd == "cluster" && a.size() == 2)
  {
    GtsSurface* in = read_surface(a[0]);
    s = make_cluster(in, std::strtod(a[1], nullptr));
    gts_object_destroy(GTS_OBJECT(in));
  }
  else if (cmd == "distance" && (a.size() == 2 || a.size() == 3))
  {
    GtsSurface* sa = read_surface(a[0]);
    GtsSurface* sb = read_surface(a[1]);
    const double delta = a.size() == 3 ? std::strtod(a[2], nullptr) : 0.1;
    // gts_range_reset leaves mean and stddev as they were.
    GtsRange face{}, boundary{};
    gts_surface_distance(sa, sb, delta, &face, &boundary);
    print_range("distance_a_to_b", face);
    print_range("boundary_distance_a_to_b", boundary);
    face = boundary = GtsRange{};
    gts_surface_distance(sb, sa, delta, &face, &boundary);
    print_range("distance_b_to_a", face);
    print_range("boundary_distance_b_to_a", boundary);
    gts_object_destroy(GTS_OBJECT(sa));
    gts_object_destroy(GTS_OBJECT(sb));
  }
  else if (cmd == "curvature" && a.size() == 1)
  {
    GtsSurface* in = read_surface(a[0]);
    Lines lines;
    CurvatureArg arg{ in, &lines };
    gts_surface_foreach_vertex(in, (GtsFunc) curvature_vertex, &arg);
    print_lines(lines);
    gts_object_destroy(GTS_OBJECT(in));
  }
  else if (cmd == "inside" && a.size() == 2)
  {
    GtsSurface* in = read_surface(a[0]);
    inside_cmd(in, std::strtoul(a[1], nullptr, 10));
    gts_object_destroy(GTS_OBJECT(in));
  }
  else if (cmd == "merge" && a.size() == 2)
  {
    GtsSurface* in = read_surface(a[0]);
    s = merge_cmd(in, std::strtod(a[1], nullptr));
    gts_object_destroy(GTS_OBJECT(in));
  }
  else if (cmd == "centers" && a.size() == 1)
  {
    GtsSurface* in = read_surface(a[0]);
    centers_cmd(in);
    gts_object_destroy(GTS_OBJECT(in));
  }
  else
  {
    usage();
  }

  if (s)
  {
    print_summary(s);
    if (out) write_surface(s, out);
    gts_object_destroy(GTS_OBJECT(s));
  }

  std::println("glib_critical {}", g_n_critical);
  std::println("glib_warning {}",  g_n_warning);
  std::println("status ok");
  return 0;
}
