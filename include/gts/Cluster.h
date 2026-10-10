// ClusterGrid: out-of-core simplification by vertex clustering, as GTS's
// GtsClusterGrid (Lindstrom, "Out-of-core simplification of large polygonal
// models", 2000). Space is cut into cubes of a given size. Each triangle
// added maps its corners to the clusters of their cubes, and becomes a face
// of the surface when the three clusters differ. update() moves each
// cluster's vertex to the mean of the points added to it.
//
// Copyright (C) 1999 Stéphane Popinet (GTS: oocs.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "gts/Geometrium.h"
#include "gts/Stats.h"

#include <array>
#include <map>

namespace gts {

class ClusterGrid
{
public:
  // gts_cluster_grid_new: cubes of size delta over bbox, which must contain
  // every point added. The faces go to s. Throws std::invalid_argument for
  // delta <= 0 or a grid too large to index.
  ClusterGrid(Surface& s, const BBox& bbox, double delta);
  // Removes the vertices of clusters that no face uses.
  ~ClusterGrid();
  ClusterGrid(const ClusterGrid&) = delete;
  ClusterGrid& operator=(const ClusterGrid&) = delete;

  // gts_cluster_grid_add_triangle.
  void add_triangle(Vec3 p1, Vec3 p2, Vec3 p3);

  // gts_cluster_grid_update: moves each cluster's vertex to the mean of its
  // points. Returns the statistics of the number of points per cluster.
  Range update();

private:
  using Id = std::array<unsigned, 3>;
  struct Cluster { VertexId v; Vec3 sum; unsigned n = 0; };

  Cluster& add_point(Vec3 p);

  Surface&                m_surface;
  BBox                    m_bbox;
  std::array<double, 3>   m_size;
  std::map<Id, Cluster>   m_clusters;
};

} // namespace gts
