// BBTree: a tree of bounding boxes of the faces of a surface, or of a set
// of edges, as GTS's bounding-box trees: for finding the pairs of faces
// whose boxes overlap, the faces a ray crosses, and the face or edge
// closest to a point.
//
// Copyright (C) 1999 Stéphane Popinet (GTS: bbtree.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "gts/Geometrium.h"

#include <functional>
#include <iosfwd>
#include <vector>

namespace gts {

// gts_bbox_triangle: the bounding box of a face.
BBox face_bbox(const Geometrium& g, FaceId f);

// gts_bbox_segment: the bounding box of an edge.
BBox edge_bbox(const Geometrium& g, EdgeId e);

// gts_bbox_diagonal2: the square of the length of the box's diagonal.
inline double diagonal2(const BBox& b)
{
  const double x = b.max.x - b.min.x, y = b.max.y - b.min.y, z = b.max.z - b.min.z;
  return x * x + y * y + z * z;
}

// gts_bbox_point_distance2: lower and upper bounds of the square of the
// distance from p to the object in b, which must touch every side of b.
struct DistanceBounds { double min, max; };
DistanceBounds point_distance2(const BBox& b, Vec3 p);

// gts_bboxes_are_overlapping: boxes that touch overlap.
inline bool overlapping(const BBox& a, const BBox& b)
{
  return ! (a.min.x > b.max.x || b.min.x > a.max.x ||
            a.min.y > b.max.y || b.min.y > a.max.y ||
            a.min.z > b.max.z || b.min.z > a.max.z);
}

// gts_bbox_is_stabbed: whether the ray from p towards +x meets b.
bool is_stabbed(const BBox& b, Vec3 p);

// gts_bbox_overlaps_triangle and gts_bbox_overlaps_segment: whether b and
// t, or b and s, overlap, by the separating axis test of Tomas
// Akenine-Möller. GTS took a face's points from its first edge, then the
// vertex opposite it.
bool overlaps(const BBox& b, const Triangle& t);
bool overlaps(const BBox& b, const Segment& s);

// gts_bbox_draw: b as a Geomview OFF object of 8 vertices and 6 faces.
void write_oogl(const BBox& b, std::ostream& out);

// A leaf of a tree: a face, or an edge in a tree of edges.
struct BBLeaf
{
  FaceId face;
  EdgeId edge;
};

class BBTree
{
public:
  // A distance from a point to a leaf, and the point of a leaf closest to
  // a point. By default the Euclidean distance to the face or edge.
  using DistanceFunc = std::function<double(Vec3, BBLeaf)>;
  using ClosestFunc  = std::function<Vec3(Vec3, BBLeaf)>;

  // gts_bb_tree_surface: a tree of the faces of s.
  explicit BBTree(const Surface& s);
  // gts_bb_tree_new of the boxes of the edges, given in the order GTS
  // visited them.
  BBTree(const Geometrium& g, const std::vector<EdgeId>& edges);

  const Geometrium& geometrium() const { return *m_geo; }
  // The box of the root; the tree must not be empty.
  const BBox& box() const { return m_nodes[m_root].box; }

  // gts_bb_tree_stabbed: the leaves whose boxes the ray from p towards +x
  // meets.
  std::vector<BBLeaf> stabbed(Vec3 p) const;

  // gts_bb_tree_point_closest_bboxes: leaves, one of which holds the object
  // closest to p.
  std::vector<BBLeaf> closest_leaves(Vec3 p) const;

  // gts_bb_tree_point_distance: the distance, by d, from p to the closest
  // leaf; that leaf is stored in closest when given. The leaf with the
  // smallest absolute value wins, the first of equal ones.
  double point_distance(Vec3 p, const DistanceFunc& d = {}, BBLeaf* closest = nullptr) const;

  // gts_bb_tree_point_closest: the point of a leaf closest to p, found with
  // c. distance2, when given, gets the square of its distance from p.
  Vec3 point_closest(Vec3 p, const ClosestFunc& c = {}, double* distance2 = nullptr) const;

  // gts_bb_tree_overlap: the leaves whose boxes overlap b, in the tree's
  // order.
  std::vector<BBLeaf> overlap(const BBox& b) const;
  // gts_bb_tree_is_overlapping: whether a leaf's box overlaps b.
  bool is_overlapping(const BBox& b) const;

  // gts_bb_tree_draw: the boxes of the nodes at depth, 1 for the root, as a
  // Geomview LIST.
  void write_oogl(unsigned depth, std::ostream& out) const;

  // gts_bb_tree_traverse_overlapping: calls f(face1, face2) for each pair of
  // a face of a and a face of b whose boxes overlap, in GTS's order. The two
  // trees may be the same, which gives each pair twice and each face with
  // itself.
  static void traverse_overlapping(const BBTree& a, const BBTree& b,
                                   const std::function<void(FaceId, FaceId)>& f);

  bool empty() const { return m_nodes.empty(); }

private:
  struct Node
  {
    BBox   box;
    BBLeaf leaf;                  // a leaf's face or edge
    int    child[2] = { -1, -1 }; // GTS's order: the negative side first
  };

  int build(std::vector<int> leaves);
  static void traverse(const BBTree& a, int na, const BBTree& b, int nb,
                       const std::function<void(FaceId, FaceId)>& f);
  void stabbed(int n, Vec3 p, std::vector<BBLeaf>& list) const;
  void overlap(int n, const BBox& b, std::vector<BBLeaf>& list) const;
  bool is_overlapping(int n, const BBox& b) const;
  void write_oogl(int n, unsigned d, unsigned depth, std::ostream& out) const;
  double min_max(int n, Vec3 p, double min_max, std::vector<int>& list) const;

  const Geometrium* m_geo;
  std::vector<Node> m_nodes;
  int               m_root = -1;
};

} // namespace gts
