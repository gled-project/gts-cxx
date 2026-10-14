// Copyright (C) 1999 Stéphane Popinet (GTS: cdt.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-13.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/FacePoints.h"

#include "gts/Delaunay.h"
#include "gts/Predicates.h"
#include "gts/Topology.h"

#include <algorithm>
#include <utility>

namespace gts {

FacePoints::FacePoints(Surface& s) : m_s(&s)
{
  Geometrium& g = s.geometrium();
  m_on_face = s.on_face([this](const FaceEvent& e) { e.added ? added(e.face) : removed(e.face); });
  m_on_split = g.on_split([this](const SplitEvent& e) {
    for (const FaceSplitEvent& f : e.faces) changed(f.pieces);
  });
  m_on_collapse = g.on_collapse([this](const CollapseEvent& e) {
    if (e.vertex.valid()) changed(vertex_faces(*m_s, e.vertex));
    else                  changed({});
  });
  m_on_renumber = g.on_renumber([this](const Renumbering& r) { renumber(r); });
}

FacePoints::PointId FacePoints::add(Vec3 p, FaceId f)
{
  const PointId i = PointId(m_p.size());
  m_p.push_back(p);
  m_next.push_back(kNoPoint);
  m_removed.push_back(false);
  append(f, i);
  return i;
}

void FacePoints::remove(PointId i)
{
  m_removed[i] = true;
}

std::generator<FacePoints::PointId> FacePoints::points(FaceId f)
{
  divide();
  place();
  if (f.i >= m_head.size()) co_return;
  for (PointId i = m_head[f.i]; i != kNoPoint; i = m_next[i])
    if ( ! m_removed[i]) co_yield i;
}

void FacePoints::append(FaceId f, PointId i)
{
  if (m_head.size() <= f.i)
  {
    const std::size_t n = std::max<std::size_t>(f.i + 1, 2 * m_head.size());
    m_head.resize(n, kNoPoint);
    m_tail.resize(n, kNoPoint);
  }
  m_next[i] = kNoPoint;
  if (m_tail[f.i] == kNoPoint) m_head[f.i] = i;
  else                         m_next[m_tail[f.i]] = i;
  m_tail[f.i] = i;
}

void FacePoints::take_list(FaceId f, std::vector<PointId>& to)
{
  if (f.i >= m_head.size()) return;
  for (PointId i = m_head[f.i]; i != kNoPoint; i = m_next[i])
    if ( ! m_removed[i]) to.push_back(i);
  m_head[f.i] = m_tail[f.i] = kNoPoint;
}

void FacePoints::added(FaceId f)
{
  if ( ! m_left.empty()) m_came.push_back(f);
}

// A face that leaves after faces came ends a change: the one before is
// divided first.
void FacePoints::removed(FaceId f)
{
  if ( ! m_came.empty()) divide();
  m_left.push_back(m_s->geometrium().vertices_edges(f).v);
  take_list(f, m_pending);
}

void FacePoints::changed(std::span<const FaceId> faces)
{
  divide();
  for (FaceId f : faces)
    if (m_s->contains(f))
    {
      take_list(f, m_pending);
      m_came.push_back(f);
    }
  divide();
}

void FacePoints::divide()
{
  std::vector<PointId> pending = std::exchange(m_pending, {});
  std::vector<std::array<VertexId, 3>> left = std::exchange(m_left, {});
  std::vector<FaceId> came = std::exchange(m_came, {});
  if (pending.empty()) return;

  const Geometrium& g = m_s->geometrium();
  // A face that came may be leaving now: out of s, but still linked.
  std::erase_if(came, [&](FaceId f) { return g.dead(f); });
  auto pos = [&](VertexId v) { return g.vertex(v).p; };
  auto has = [](const std::array<VertexId, 3>& a, VertexId v) { return a[0] == v || a[1] == v || a[2] == v; };

  // A face split at vertex v into three, as gts_delaunay_add_vertex_to_face:
  // nf[j] joins v to the edge from vertex j to vertex j + 1 of the old face.
  if (left.size() == 1 && came.size() == 3)
  {
    const auto& L = left[0];
    VertexId v;
    for (VertexId w : g.vertices_of(came[0])) if ( ! has(L, w)) v = w;
    std::array<FaceId, 3> nf{};
    for (FaceId f : came)
    {
      const auto W = g.vertices_of(f);
      if (has(W, v))
        for (int j = 0; j < 3; ++j)
          if (has(W, L[j]) && has(W, L[(j + 1) % 3])) nf[j] = f;
    }
    if (v.valid() && nf[0].valid() && nf[1].valid() && nf[2].valid())
    {
      const Vec3 pv = pos(v), p1 = pos(L[0]), p2 = pos(L[1]), p3 = pos(L[2]);
      for (PointId i : pending)
      {
        const Vec3 p = m_p[i];
        int j;
        if (orient2d(pv, p1, p) >= 0.)
        {
          const double o = orient2d(pv, p2, p);
          j = o != 0. ? (o > 0. ? 1 : 0) : (orient2d(pv, p3, p) > 0. ? 0 : 1);
        }
        else
          j = orient2d(pv, p3, p) > 0. ? 2 : 1;
        append(nf[j], i);
      }
      return;
    }
  }

  // A swap, as swap_if_in_circle: v3 and v4 are the vertices of the old
  // faces off the old edge; the points on the left of v3 -> v4, and on that
  // line, go to the first new face.
  if (left.size() == 2 && came.size() == 2)
  {
    VertexId v3, v4;
    for (VertexId w : left[0]) if ( ! has(left[1], w)) v3 = w;
    for (VertexId w : left[1]) if ( ! has(left[0], w)) v4 = w;
    const auto W0 = g.vertices_of(came[0]), W1 = g.vertices_of(came[1]);
    if (v3.valid() && v4.valid() && has(W0, v3) && has(W0, v4) && has(W1, v3) && has(W1, v4))
    {
      const Vec3 p3 = pos(v3), p4 = pos(v4);
      for (PointId i : pending) append(orient2d(p3, p4, m_p[i]) >= 0. ? came[0] : came[1], i);
      return;
    }
  }

  // Any other change: the first face that came and holds the point, or
  // place() later.
  for (PointId i : pending)
  {
    const Vec3 p = m_p[i];
    FaceId to;
    for (FaceId f : came)
      if (point_in_triangle(p, g.triangle(f)) != Intersect::Out) { to = f; break; }
    if (to.valid()) append(to, i);
    else            m_unplaced.push_back(i);
  }
}

void FacePoints::place()
{
  FaceId guess;
  for (PointId i : std::exchange(m_unplaced, {}))
  {
    if (m_removed[i]) continue;
    FaceId f = locate(*m_s, m_p[i], guess);
    if ( ! f.valid() && guess.valid()) f = locate(*m_s, m_p[i]);
    if (f.valid()) { append(f, i); guess = f; }
    else           m_removed[i] = true;
  }
}

void FacePoints::renumber(const Renumbering& r)
{
  const std::size_t n = m_s->geometrium().face_capacity();
  std::vector<PointId> head(n, kNoPoint), tail(n, kNoPoint);
  for (std::size_t i = 0; i < m_head.size() && i < r.faces.size(); ++i)
    if (m_head[i] != kNoPoint)
      if (const FaceId f = r(FaceId(Index(i))); f.valid())
      {
        head[f.i] = m_head[i];
        tail[f.i] = m_tail[i];
      }
  m_head = std::move(head);
  m_tail = std::move(tail);
  for (FaceId& f : m_came) f = r(f);
  std::erase_if(m_came, [](FaceId f) { return ! f.valid(); });
  for (auto& L : m_left)
    for (VertexId& v : L)
      if (v.valid()) v = r(v);
}

} // namespace gts
