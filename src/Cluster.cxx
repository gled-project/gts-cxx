// Copyright (C) 1999 Stéphane Popinet (GTS: oocs.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/Cluster.h"

#include <climits>
#include <cmath>
#include <stdexcept>

namespace gts {

ClusterGrid::ClusterGrid(Surface& s, const BBox& bbox, double delta) :
  m_surface(s), m_bbox(bbox)
{
  if ( ! (delta > 0.))
    throw std::invalid_argument("gts::ClusterGrid: the size of the cubes must be positive");
  m_size = { std::ceil((bbox.max.x - bbox.min.x) / delta),
             std::ceil((bbox.max.y - bbox.min.y) / delta),
             std::ceil((bbox.max.z - bbox.min.z) / delta) };
  for (double n : m_size)
    if ( ! (n <= 2. * INT_MAX + 2.))
      throw std::invalid_argument("gts::ClusterGrid: too many cubes");
}

ClusterGrid::~ClusterGrid()
{
  // cluster_destroy: the vertex goes when nothing uses it.
  Geometrium& g = m_surface.geometrium();
  for (auto& [id, c] : m_clusters)
    if ( ! g.dead(c.v) && ! g.vertex(c.v).edge.valid())
      g.kill_vertex(c.v);
}

// cluster_index and cluster_grid_add_point. A point outside the box goes to
// cluster (0, 0, 0); see QUIRKS.md.
ClusterGrid::Cluster& ClusterGrid::add_point(Vec3 p)
{
  const BBox& bb = m_bbox;
  Id id = { 0, 0, 0 };
  if (p.x >= bb.min.x && p.x <= bb.max.x && p.y >= bb.min.y && p.y <= bb.max.y && p.z >= bb.min.z && p.z <= bb.max.z)
    id = { unsigned(p.x == bb.max.x ? m_size[0] - 1. : m_size[0] * (p.x - bb.min.x) / (bb.max.x - bb.min.x)),
           unsigned(p.y == bb.max.y ? m_size[1] - 1. : m_size[1] * (p.y - bb.min.y) / (bb.max.y - bb.min.y)),
           unsigned(p.z == bb.max.z ? m_size[2] - 1. : m_size[2] * (p.z - bb.min.z) / (bb.max.z - bb.min.z)) };

  auto [it, made] = m_clusters.try_emplace(id);
  Cluster& c = it->second;
  if (made)
    c.v = m_surface.geometrium().add_vertex({ 0., 0., 0. });
  // cluster_add, summed apart from the vertex.
  c.sum.x += p.x;
  c.sum.y += p.y;
  c.sum.z += p.z;
  c.n++;
  return c;
}

void ClusterGrid::add_triangle(Vec3 p1, Vec3 p2, Vec3 p3)
{
  Geometrium& g = m_surface.geometrium();
  const VertexId v1 = add_point(p1).v;
  const VertexId v2 = add_point(p2).v;
  const VertexId v3 = add_point(p3).v;
  if (v1 == v2 || v2 == v3 || v3 == v1)
    return;

  bool new_edge = false;
  auto edge = [&](VertexId a, VertexId b) {
    const EdgeId e = g.find_edge(a, b);
    if (e.valid()) return e;
    new_edge = true;
    return g.add_edge(a, b);
  };
  const EdgeId e1 = edge(v1, v2), e2 = edge(v2, v3), e3 = edge(v3, v1);
  // gts_triangle_use_edges: a face of these three edges, in any surface.
  bool used = false;
  for (FaceId t : g.faces_of(e1))
  {
    const auto& E = g.face(t).e;
    const bool has2 = E[0] == e2 || E[1] == e2 || E[2] == e2;
    const bool has3 = E[0] == e3 || E[1] == e3 || E[2] == e3;
    if (has2 && has3) { used = true; break; }
  }
  if (new_edge || ! used)
    m_surface.add_face(g.add_face(e1, e2, e3));
}

Range ClusterGrid::update()
{
  Geometrium& g = m_surface.geometrium();
  Range stats;
  for (auto& [id, c] : m_clusters)
  {
    // cluster_update
    Vec3 p = c.sum;
    if (c.n > 1)
      p = { p.x / c.n, p.y / c.n, p.z / c.n };
    g.position(c.v) = p;
    stats.add(c.n);
  }
  stats.update();
  return stats;
}

} // namespace gts
