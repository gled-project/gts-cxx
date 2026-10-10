// KdTree: GTS's kd-tree of vertices, for finding the vertices in a box.
//
// Copyright (C) 1999 Stéphane Popinet (GTS: kdtree.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "gts/Geometrium.h"

#include <span>
#include <vector>

namespace gts {

class KdTree
{
public:
  // gts_kdtree_new: split at the median by x, y and z in turn. Vertices
  // with equal coordinates keep the order they are given in.
  KdTree(const Geometrium& g, std::vector<VertexId> vertices);

  // gts_kdtree_range: the vertices inside b, its sides included, in GTS's
  // order.
  std::vector<VertexId> range(const BBox& b) const;

private:
  struct Node
  {
    VertexId v;
    int      child[2] = { -1, -1 }; // GTS's order: the upper half first
  };

  int  build(std::span<VertexId> points, int axis);
  void range(int n, const BBox& b, int axis, std::vector<VertexId>& list) const;

  const Geometrium* m_geo;
  std::vector<Node> m_nodes;
  int               m_root = -1;
};

} // namespace gts
