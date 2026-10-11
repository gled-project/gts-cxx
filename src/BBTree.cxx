// Copyright (C) 1999 Stéphane Popinet (GTS: bbtree.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/BBTree.h"
#include "gts/Distance.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <format>
#include <ostream>

namespace gts {

BBox face_bbox(const Geometrium& g, FaceId f)
{
  const Edge& e1 = g.edge(g.face(f).e[0]);
  BBox b;
  b.min = b.max = g.vertex(e1.v[0]).p;
  for (VertexId v : { e1.v[1], g.opposite_vertex(f) })
  {
    const Vec3 p = g.vertex(v).p;
    if (p.x > b.max.x) b.max.x = p.x;
    if (p.x < b.min.x) b.min.x = p.x;
    if (p.y > b.max.y) b.max.y = p.y;
    if (p.y < b.min.y) b.min.y = p.y;
    if (p.z > b.max.z) b.max.z = p.z;
    if (p.z < b.min.z) b.min.z = p.z;
  }
  return b;
}

BBox edge_bbox(const Geometrium& g, EdgeId e)
{
  const Segment s = g.segment(e);
  BBox b;
  b.min = { std::fmin(s.a.x, s.b.x), std::fmin(s.a.y, s.b.y), std::fmin(s.a.z, s.b.z) };
  b.max = { std::fmax(s.a.x, s.b.x), std::fmax(s.a.y, s.b.y), std::fmax(s.a.z, s.b.z) };
  return b;
}

DistanceBounds point_distance2(const BBox& b, Vec3 p)
{
  const double xd1 = (b.min.x - p.x)*(b.min.x - p.x), xd2 = (p.x - b.max.x)*(p.x - b.max.x);
  const double yd1 = (b.min.y - p.y)*(b.min.y - p.y), yd2 = (p.y - b.max.y)*(p.y - b.max.y);
  const double zd1 = (b.min.z - p.z)*(b.min.z - p.z), zd2 = (p.z - b.max.z)*(p.z - b.max.z);

  double dmin = p.x < b.min.x ? xd1 : p.x > b.max.x ? xd2 : 0.0;
  dmin += p.y < b.min.y ? yd1 : p.y > b.max.y ? yd2 : 0.0;
  dmin += p.z < b.min.z ? zd1 : p.z > b.max.z ? zd2 : 0.0;

  auto minmax = [](double a, double c, double& lo, double& hi) { if (a < c) { lo = a; hi = c; } else { lo = c; hi = a; } };
  double mx, Mx, my, My, mz, Mz;
  minmax(xd1, xd2, mx, Mx);
  minmax(yd1, yd2, my, My);
  minmax(zd1, zd2, mz, Mz);

  double dmax = mx + My + Mz;
  dmax = std::min(dmax, Mx + my + Mz);
  dmax = std::min(dmax, Mx + My + mz);
  return { dmin, dmax };
}

namespace {

double volume(const BBox& b) { return (b.max.x - b.min.x) * (b.max.y - b.min.y) * (b.max.z - b.min.z); }

double lo(const BBox& b, int dir) { return dir == 0 ? b.min.x : dir == 1 ? b.min.y : b.min.z; }
double hi(const BBox& b, int dir) { return dir == 0 ? b.max.x : dir == 1 ? b.max.y : b.max.z; }

} // namespace

// GTS prepended the boxes of the faces in the order its foreach visited
// them, so its list runs in the reverse order of the faces.
BBTree::BBTree(const Surface& s) : m_geo(&s.geometrium())
{
  const Geometrium& g = s.geometrium();
  std::vector<int> leaves;
  for (FaceId f : s.faces())
  {
    m_nodes.push_back({ face_bbox(g, f), { f, EdgeId() }, { -1, -1 } });
    leaves.push_back(int(m_nodes.size()) - 1);
  }
  std::reverse(leaves.begin(), leaves.end());
  if ( ! leaves.empty())
    m_root = build(std::move(leaves));
}

BBTree::BBTree(const Geometrium& g, const std::vector<EdgeId>& edges) : m_geo(&g)
{
  std::vector<int> leaves;
  for (EdgeId e : edges)
  {
    m_nodes.push_back({ edge_bbox(g, e), { FaceId(), e }, { -1, -1 } });
    leaves.push_back(int(m_nodes.size()) - 1);
  }
  std::reverse(leaves.begin(), leaves.end());
  if ( ! leaves.empty())
    m_root = build(std::move(leaves));
}

// gts_bb_tree_new. The lists are GSLists in GTS, built by prepending; the
// vectors here hold them in the same order.
int BBTree::build(std::vector<int> list)
{
  if (list.size() == 1)
    return list[0];

  // gts_bbox_bboxes
  BBox bbox = m_nodes[list[0]].box;
  for (std::size_t i = 1; i < list.size(); ++i)
  {
    const BBox& bb = m_nodes[list[i]].box;
    if (bb.min.x < bbox.min.x) bbox.min.x = bb.min.x;
    if (bb.min.y < bbox.min.y) bbox.min.y = bb.min.y;
    if (bb.min.z < bbox.min.z) bbox.min.z = bb.min.z;
    if (bb.max.x > bbox.max.x) bbox.max.x = bb.max.x;
    if (bb.max.y > bbox.max.y) bbox.max.y = bb.max.y;
    if (bb.max.z > bbox.max.z) bbox.max.z = bb.max.z;
  }

  int dir;
  if (bbox.max.x - bbox.min.x > bbox.max.y - bbox.min.y)
    dir = bbox.max.z - bbox.min.z > bbox.max.x - bbox.min.x ? 2 : 0;
  else if (bbox.max.z - bbox.min.z > bbox.max.y - bbox.min.y)
    dir = 2;
  else
    dir = 1;

  const double cut = (lo(bbox, dir) + hi(bbox, dir)) / 2.;
  std::vector<int> positive, negative;
  for (int n : list)
  {
    const BBox& bb = m_nodes[n].box;
    if ((lo(bb, dir) + hi(bb, dir)) / 2. > cut) positive.insert(positive.begin(), n);
    else                                        negative.insert(negative.begin(), n);
  }
  // All on one side: the first half of the list stays, the rest goes over.
  if (positive.empty())
  {
    const std::size_t keep = (negative.size() - 1) / 2 + 1;
    positive.assign(negative.begin() + keep, negative.end());
    negative.resize(keep);
  }
  else if (negative.empty())
  {
    const std::size_t keep = (positive.size() - 1) / 2 + 1;
    negative.assign(positive.begin() + keep, positive.end());
    positive.resize(keep);
  }

  // g_node_prepend of the positive side, then of the negative one.
  const int p = build(std::move(positive));
  const int q = build(std::move(negative));
  m_nodes.push_back({ bbox, {}, { q, p } });
  return int(m_nodes.size()) - 1;
}

void BBTree::traverse_overlapping(const BBTree& a, const BBTree& b, const std::function<void(FaceId, FaceId)>& f)
{
  if (a.m_root >= 0 && b.m_root >= 0)
    traverse(a, a.m_root, b, b.m_root, f);
}

void BBTree::traverse(const BBTree& a, int na, const BBTree& b, int nb, const std::function<void(FaceId, FaceId)>& f)
{
  const Node& n1 = a.m_nodes[na];
  const Node& n2 = b.m_nodes[nb];
  // GTS compared the boxes as pointers first: a box always overlaps itself.
  if ( ! (&n1 == &n2 || overlapping(n1.box, n2.box)))
    return;

  const bool leaf1 = n1.child[0] < 0, leaf2 = n2.child[0] < 0;
  if (leaf1 && leaf2)
    f(n1.leaf.face, n2.leaf.face);
  else if (leaf2 || ( ! leaf1 && volume(n1.box) > volume(n2.box)))
  {
    for (int c : n1.child) traverse(a, c, b, nb, f);
  }
  else
  {
    for (int c : n2.child) traverse(a, na, b, c, f);
  }
}

//==============================================================================
// Queries by point
//==============================================================================

bool is_stabbed(const BBox& b, Vec3 p)
{
  return ! (p.x > b.max.x || p.y < b.min.y || p.y > b.max.y || p.z < b.min.z || p.z > b.max.z);
}

std::vector<BBLeaf> BBTree::stabbed(Vec3 p) const
{
  std::vector<BBLeaf> list;
  if (m_root >= 0) stabbed(m_root, p, list);
  return list;
}

void BBTree::stabbed(int n, Vec3 p, std::vector<BBLeaf>& list) const
{
  const Node& node = m_nodes[n];
  if ( ! is_stabbed(node.box, p)) return;
  if (node.child[0] < 0) { list.push_back(node.leaf); return; }
  for (int c : node.child) stabbed(c, p, list);
}

// bb_tree_min_max: list gets the leaves in the order GTS prepended them,
// so reversed.
double BBTree::min_max(int n, Vec3 p, double min_max, std::vector<int>& list) const
{
  const Node& node = m_nodes[n];
  if (node.child[0] < 0)
  {
    list.push_back(n);
    return min_max;
  }
  const int t1 = node.child[0], t2 = node.child[1];
  const DistanceBounds d1 = point_distance2(m_nodes[t1].box, p);
  if (d1.max < min_max) min_max = d1.max;
  const DistanceBounds d2 = point_distance2(m_nodes[t2].box, p);
  if (d2.max < min_max) min_max = d2.max;

  if (d1.min < d2.min)
  {
    if (d1.min <= min_max)
    {
      min_max = this->min_max(t1, p, min_max, list);
      if (d2.min <= min_max) min_max = this->min_max(t2, p, min_max, list);
    }
  }
  else
  {
    if (d2.min <= min_max)
    {
      min_max = this->min_max(t2, p, min_max, list);
      if (d1.min <= min_max) min_max = this->min_max(t1, p, min_max, list);
    }
  }
  return min_max;
}

std::vector<BBLeaf> BBTree::closest_leaves(Vec3 p) const
{
  std::vector<BBLeaf> leaves;
  if (m_root < 0) return leaves;
  std::vector<int> list;
  double mm = point_distance2(m_nodes[m_root].box, p).max;
  mm = min_max(m_root, p, mm, list);
  for (auto i = list.rbegin(); i != list.rend(); ++i)
    if ( ! (point_distance2(m_nodes[*i].box, p).min > mm))
      leaves.push_back(m_nodes[*i].leaf);
  return leaves;
}

double BBTree::point_distance(Vec3 p, const DistanceFunc& d, BBLeaf* closest) const
{
  const Geometrium& g = *m_geo;
  double dmin = DBL_MAX;
  for (const BBLeaf& l : closest_leaves(p))
  {
    const double x = d ? d(p, l) : l.face.valid() ? distance(g, p, l.face) : distance(g, p, l.edge);
    if (std::fabs(x) < std::fabs(dmin))
    {
      dmin = x;
      if (closest) *closest = l;
    }
  }
  return dmin;
}

Vec3 BBTree::point_closest(Vec3 p, const ClosestFunc& c, double* distance2) const
{
  const Geometrium& g = *m_geo;
  double dmin = DBL_MAX;
  Vec3 np;
  for (const BBLeaf& l : closest_leaves(p))
  {
    const Vec3 tp = c ? c(p, l) : l.face.valid() ? closest_point(g, p, l.face) : closest_point(g, p, l.edge);
    const double x = gts::distance2(tp, p);
    if (x < dmin)
    {
      np = tp;
      dmin = x;
    }
  }
  if (distance2) *distance2 = dmin;
  return np;
}

//==============================================================================
// Queries by box
//==============================================================================

std::vector<BBLeaf> BBTree::overlap(const BBox& b) const
{
  std::vector<BBLeaf> list;
  if (m_root >= 0) overlap(m_root, b, list);
  return list;
}

void BBTree::overlap(int n, const BBox& b, std::vector<BBLeaf>& list) const
{
  const Node& node = m_nodes[n];
  if ( ! overlapping(b, node.box)) return;
  if (node.child[0] < 0) { list.push_back(node.leaf); return; }
  for (int c : node.child) overlap(c, b, list);
}

bool BBTree::is_overlapping(const BBox& b) const
{
  return m_root >= 0 && is_overlapping(m_root, b);
}

bool BBTree::is_overlapping(int n, const BBox& b) const
{
  const Node& node = m_nodes[n];
  if ( ! overlapping(b, node.box)) return false;
  if (node.child[0] < 0) return true;
  for (int c : node.child)
    if (is_overlapping(c, b)) return true;
  return false;
}

//==============================================================================
// Box and triangle, from GTS's tribox3.c
//==============================================================================
//
// AABB-triangle overlap test by Tomas Akenine-Möller (2001), changed from
// float to double by Stéphane Popinet (2004). GTS kept fabsf(), which
// rounds the edge components to float; gts-cxx uses fabs(); see QUIRKS.md.

namespace {

bool plane_box_overlap(const double normal[3], const double vert[3], const double maxbox[3])
{
  double vmin[3], vmax[3];
  for (int q = 0; q < 3; q++)
  {
    const double v = vert[q];
    if (normal[q] > 0.0) { vmin[q] = -maxbox[q] - v; vmax[q] =  maxbox[q] - v; }
    else                 { vmin[q] =  maxbox[q] - v; vmax[q] = -maxbox[q] - v; }
  }
  auto dot = [](const double* a, const double* b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; };
  if (dot(normal, vmin) > 0.0) return false;
  return dot(normal, vmax) >= 0.0;
}

bool tri_box_overlap(const double boxcenter[3], const double boxhalfsize[3], const double triverts[3][3])
{
  enum { X, Y, Z };
  double v0[3], v1[3], v2[3], e0[3], e1[3], e2[3];
  for (int k = 0; k < 3; ++k)
  {
    v0[k] = triverts[0][k] - boxcenter[k];
    v1[k] = triverts[1][k] - boxcenter[k];
    v2[k] = triverts[2][k] - boxcenter[k];
  }
  for (int k = 0; k < 3; ++k)
  {
    e0[k] = v1[k] - v0[k];
    e1[k] = v2[k] - v1[k];
    e2[k] = v0[k] - v2[k];
  }
  const double* h = boxhalfsize;
  double min, max, p0, p1, p2, rad;
  auto out = [&] { return min > rad || max < -rad; };
  auto minmax = [&](double a, double b) { if (a < b) { min = a; max = b; } else { min = b; max = a; } };

  // AXISTEST_X01, X2, Y02, Y1, Z12 and Z0 of the original.
  auto x01 = [&](double a, double b, double fa, double fb) {
    p0 = a*v0[Y] - b*v0[Z]; p2 = a*v2[Y] - b*v2[Z]; minmax(p0, p2); rad = fa*h[Y] + fb*h[Z]; return out(); };
  auto x2 = [&](double a, double b, double fa, double fb) {
    p0 = a*v0[Y] - b*v0[Z]; p1 = a*v1[Y] - b*v1[Z]; minmax(p0, p1); rad = fa*h[Y] + fb*h[Z]; return out(); };
  auto y02 = [&](double a, double b, double fa, double fb) {
    p0 = -a*v0[X] + b*v0[Z]; p2 = -a*v2[X] + b*v2[Z]; minmax(p0, p2); rad = fa*h[X] + fb*h[Z]; return out(); };
  auto y1 = [&](double a, double b, double fa, double fb) {
    p0 = -a*v0[X] + b*v0[Z]; p1 = -a*v1[X] + b*v1[Z]; minmax(p0, p1); rad = fa*h[X] + fb*h[Z]; return out(); };
  auto z12 = [&](double a, double b, double fa, double fb) {
    p1 = a*v1[X] - b*v1[Y]; p2 = a*v2[X] - b*v2[Y];
    if (p2 < p1) { min = p2; max = p1; } else { min = p1; max = p2; }
    rad = fa*h[X] + fb*h[Y]; return out(); };
  auto z0 = [&](double a, double b, double fa, double fb) {
    p0 = a*v0[X] - b*v0[Y]; p1 = a*v1[X] - b*v1[Y]; minmax(p0, p1); rad = fa*h[X] + fb*h[Y]; return out(); };

  // The nine cross products of the edges with the axes.
  double fex = std::fabs(e0[X]), fey = std::fabs(e0[Y]), fez = std::fabs(e0[Z]);
  if (x01(e0[Z], e0[Y], fez, fey) || y02(e0[Z], e0[X], fez, fex) || z12(e0[Y], e0[X], fey, fex)) return false;
  fex = std::fabs(e1[X]); fey = std::fabs(e1[Y]); fez = std::fabs(e1[Z]);
  if (x01(e1[Z], e1[Y], fez, fey) || y02(e1[Z], e1[X], fez, fex) || z0(e1[Y], e1[X], fey, fex)) return false;
  fex = std::fabs(e2[X]); fey = std::fabs(e2[Y]); fez = std::fabs(e2[Z]);
  if (x2(e2[Z], e2[Y], fez, fey) || y1(e2[Z], e2[X], fez, fex) || z12(e2[Y], e2[X], fey, fex)) return false;

  // The axes: the box of the triangle against the box.
  for (int k = 0; k < 3; ++k)
  {
    min = max = v0[k];
    if (v1[k] < min) min = v1[k];
    if (v1[k] > max) max = v1[k];
    if (v2[k] < min) min = v2[k];
    if (v2[k] > max) max = v2[k];
    if (min > h[k] || max < -h[k]) return false;
  }

  // The plane of the triangle.
  const double normal[3] = { e0[1]*e1[2] - e0[2]*e1[1], e0[2]*e1[0] - e0[0]*e1[2], e0[0]*e1[1] - e0[1]*e1[0] };
  return plane_box_overlap(normal, v0, h);
}

bool overlaps(const BBox& b, Vec3 p1, Vec3 p2, Vec3 p3)
{
  const double bc[3] = { (b.max.x + b.min.x)/2., (b.max.y + b.min.y)/2., (b.max.z + b.min.z)/2. };
  const double bh[3] = { (b.max.x - b.min.x)/2., (b.max.y - b.min.y)/2., (b.max.z - b.min.z)/2. };
  const double tv[3][3] = { { p1.x, p1.y, p1.z }, { p2.x, p2.y, p2.z }, { p3.x, p3.y, p3.z } };
  return tri_box_overlap(bc, bh, tv);
}

} // namespace

bool overlaps(const BBox& b, const Triangle& t) { return overlaps(b, t.a, t.b, t.c); }
bool overlaps(const BBox& b, const Segment& s)  { return overlaps(b, s.a, s.b, s.a); }

//==============================================================================
// Drawing
//==============================================================================

void write_oogl(const BBox& b, std::ostream& out)
{
  out << "OFF 8 6 12\n";
  const double x[2] = { b.min.x, b.max.x }, y[2] = { b.min.y, b.max.y }, z[2] = { b.min.z, b.max.z };
  const int c[8][3] = { {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1} };
  for (const auto& k : c) out << std::format("{:g} {:g} {:g}\n", x[k[0]], y[k[1]], z[k[2]]);
  out << "4 3 2 1 0\n"
         "4 4 5 6 7\n"
         "4 2 3 7 6\n"
         "4 0 1 5 4\n"
         "4 0 4 7 3\n"
         "4 1 2 6 5\n";
}

void BBTree::write_oogl(unsigned depth, std::ostream& out) const
{
  out << "{ LIST";
  if (m_root >= 0) write_oogl(m_root, 1, depth, out);
  out << "}\n";
}

void BBTree::write_oogl(int n, unsigned d, unsigned depth, std::ostream& out) const
{
  const Node& node = m_nodes[n];
  if (d == depth) gts::write_oogl(node.box, out);
  else if (d < depth && node.child[0] >= 0)
    for (int c : node.child) write_oogl(c, d + 1, depth, out);
}

} // namespace gts
