// Copyright (C) 1999 Stéphane Popinet (GTS: kdtree.c, vertex.c, edge.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/KdTree.h"

#include <algorithm>

namespace gts {

//==============================================================================
// KdTree
//==============================================================================

namespace {

double coord(Vec3 p, int axis) { return axis == 0 ? p.x : axis == 1 ? p.y : p.z; }

} // namespace

KdTree::KdTree(const Geometrium& g, std::vector<VertexId> vertices) : m_geo(&g)
{
  if ( ! vertices.empty())
    m_root = build(vertices, 0);
}

int KdTree::build(std::span<VertexId> points, int axis)
{
  const Geometrium& g = *m_geo;
  std::stable_sort(points.begin(), points.end(), [&](VertexId a, VertexId b) {
    return coord(g.vertex(a).p, axis) < coord(g.vertex(b).p, axis);
  });

  const std::size_t middle = (points.size() - 1)/2;
  const int n = int(m_nodes.size());
  m_nodes.push_back({ points[middle] });
  if (points.size() > 1)
  {
    // g_node_prepend of the lower half, then of the upper one; an empty
    // half is a node without a point.
    const int next = (axis + 1) % 3;
    const int lower = middle > 0 ? build(points.first(middle), next) : -1;
    const int upper = points.size() - middle - 1 > 0 ? build(points.subspan(middle + 1), next) : -1;
    m_nodes[n].child[0] = upper;
    m_nodes[n].child[1] = lower;
  }
  return n;
}

std::vector<VertexId> KdTree::range(const BBox& b) const
{
  std::vector<VertexId> list;
  if (m_root >= 0) range(m_root, b, 0, list);
  return list;
}

void KdTree::range(int n, const BBox& b, int axis, std::vector<VertexId>& list) const
{
  const Node& node = m_nodes[n];
  const Vec3 p = m_geo->vertex(node.v).p;
  // gts_bbox_point_is_inside
  if (p.x >= b.min.x && p.y >= b.min.y && p.z >= b.min.z &&
      p.x <= b.max.x && p.y <= b.max.y && p.z <= b.max.z)
    list.push_back(node.v);

  const double left = coord(b.min, axis), right = coord(b.max, axis), v = coord(p, axis);
  const int next = (axis + 1) % 3;
  if (node.child[0] >= 0 && right >= v) range(node.child[0], b, next, list);
  if (node.child[1] >= 0 && left <= v)  range(node.child[1], b, next, list);
}

//==============================================================================
// Merging
//==============================================================================

std::vector<VertexId> Geometrium::merge_vertices(const std::vector<VertexId>& vertices, double epsilon,
                                                 const std::function<bool(VertexId, VertexId)>& check)
{
  const KdTree tree(*this, vertices);
  std::vector<std::uint8_t> inactive(m_vertices.size());

  for (VertexId v : vertices)
  {
    if (inactive[v.i]) continue;
    const Vec3 p = m_vertices[v.i].p;
    BBox box;
    box.min = { p.x - epsilon, p.y - epsilon, p.z - epsilon };
    box.max = { p.x + epsilon, p.y + epsilon, p.z + epsilon };
    for (VertexId sv : tree.range(box))
    {
      if (sv == v || inactive[sv.i] || (check && ! check(sv, v))) continue;
      // An edge joining sv and v would join v to itself: it goes, with its
      // faces.
      for (EdgeId e; (e = find_edge(sv, v)).valid(); )
      {
        std::vector<FaceId> faces;
        for (FaceId f : faces_of(e)) faces.push_back(f);
        for (FaceId f : faces) destroy_face(f);
        kill_edge(e);
      }
      replace_vertex(sv, v);
      inactive[sv.i] = 1;
    }
  }

  std::vector<VertexId> kept;
  for (VertexId v : vertices)
  {
    if (inactive[v.i])
      m_vertices[v.i].flags |= kDead;
    else
      kept.push_back(v);
  }
  return kept;
}

std::vector<EdgeId> Geometrium::merge_edges(const std::vector<EdgeId>& edges)
{
  std::vector<EdgeId> kept;
  for (EdgeId e : edges)
  {
    if (dead(e)) continue;
    const EdgeId de = duplicate_of(e);
    if (de.valid())
    {
      replace_edge(e, de);
      kill_edge(e);
    }
    else
      kept.push_back(e);
  }
  return kept;
}

} // namespace gts
