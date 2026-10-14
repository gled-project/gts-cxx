// FacePoints: points kept in lists by the faces of a surface whose
// projection on the xy plane holds them, as GTS's GtsListFace, for a
// triangulation in that plane with counterclockwise faces, as Delaunay.h
// makes. When faces leave the surface, their points go to the faces that
// take their place.
//
// The store follows the surface through its face events, and the splits and
// collapses of its Geometrium, which change faces in place. A face split
// at a new vertex into three, and a swap of the edge between two faces, as
// the Delaunay functions make them, divide the points as GTS's cdt.c did,
// keeping their order. After any other change, a point goes to the first
// new or changed face that holds it. A point that none holds waits until
// the points are next read, when the surface is whole again; it then goes
// to the face locate() finds, and is dropped when locate() finds none.
// Faces changed in place other than by the splits and collapses, as by the
// vertex splits of Progressive.h, are not followed.
//
// Copyright (C) 1999 Stéphane Popinet (GTS: cdt.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-13.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "gts/Geometrium.h"

#include <array>
#include <cstdint>
#include <generator>
#include <span>
#include <vector>

namespace gts {

class FacePoints
{
public:
  using PointId = std::uint32_t;

  // Follows s until destroyed; s must outlive it.
  explicit FacePoints(Surface& s);
  FacePoints(const FacePoints&) = delete;
  FacePoints& operator=(const FacePoints&) = delete;

  // Adds p at the end of the list of f, which should hold it.
  PointId add(Vec3 p, FaceId f);
  // Takes point i out of its list, as when it becomes a vertex.
  void remove(PointId i);

  Vec3 point(PointId i) const { return m_p[i]; }
  std::size_t size() const    { return m_p.size(); }

  // The points of f, in order. The surface must not change during the walk.
  std::generator<PointId> points(FaceId f);

private:
  static constexpr PointId kNoPoint = ~PointId(0);

  void added(FaceId f);
  void removed(FaceId f);
  // Faces of s changed in place, whose points are divided again.
  void changed(std::span<const FaceId> faces);
  // Divides the points of the faces that left among the faces that came,
  // from within a callout: it does not walk the surface.
  void divide();
  // Locates the points that divide() could not place.
  void place();
  void append(FaceId f, PointId i);
  void take_list(FaceId f, std::vector<PointId>& to);
  void renumber(const Renumbering& r);

  Surface*             m_s;
  std::vector<Vec3>    m_p;
  std::vector<PointId> m_next;
  std::vector<bool>    m_removed;
  std::vector<PointId> m_head, m_tail; // by face index

  // Since the last divide(): the points of the faces that left, in order,
  // the vertices of those faces, in the order of vertices_edges(), and the
  // faces that came.
  std::vector<PointId>                 m_pending;
  std::vector<std::array<VertexId, 3>> m_left;
  std::vector<FaceId>                  m_came;
  std::vector<PointId>                 m_unplaced;

  Subscription m_on_face, m_on_split, m_on_collapse, m_on_renumber;
};

} // namespace gts
