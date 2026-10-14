// The split primitives of Geometrium: an edge split in two, with the faces on
// it, and a face split in four.
//
// Copyright (C) 1999 Stéphane Popinet (GTS: surface.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/Geometrium.h"

#include <utility>

namespace gts {

void Geometrium::split_attributes(VertexId piece, VertexId parent)
{
  for (auto& c : m_vertex_columns) c.split(piece.i, parent.i);
}

void Geometrium::split_attributes(EdgeId piece, EdgeId parent)
{
  for (auto& c : m_edge_columns) c.split(piece.i, parent.i);
}

void Geometrium::split_attributes(FaceId piece, FaceId parent)
{
  for (auto& c : m_face_columns) c.split(piece.i, parent.i);
}

void Geometrium::join_surfaces(FaceId piece, FaceId parent)
{
  for (auto& s : m_surfaces)
    if (s->contains(parent)) s->add_face(piece);
}

// As midvertex_insertion in GTS's surface.c, in the same order of operations
// on every ring.
std::array<EdgeId, 2> Geometrium::split_edge(EdgeId e, VertexId m, const std::function<void(EdgeId)>& on_new_edge)
{
  auto made = [&](EdgeId x) { split_attributes(x, e); if (on_new_edge) on_new_edge(x); };
  const EdgeId e1 = add_edge(m_edges[e.i].v[0], m);
  made(e1);
  const EdgeId e2 = add_edge(m_edges[e.i].v[1], m);
  made(e2);
  const std::array<EdgeId, 2> halves = { e1, e2 };
  const EdgeSplitEvent es = { e, m, halves };
  split_faces(e, m, halves, on_new_edge, &es);
  return halves;
}

void Geometrium::split_faces(EdgeId e, VertexId m, std::array<EdgeId, 2> halves,
                             const std::function<void(EdgeId)>& on_new_edge, const EdgeSplitEvent* es)
{
  // The faces leave the ring of e as they are split.
  std::vector<FaceId> faces;
  for (FaceId t : faces_of(e)) faces.push_back(t);

  std::vector<FaceId> pieces;
  std::vector<EdgeId> cuts;
  pieces.reserve(2 * faces.size());
  cuts.reserve(faces.size());
  for (FaceId t : faces)
  {
    EdgeId e1 = halves[0], e2 = halves[1];
    const auto [v, te] = vertices_edges(t, e);
    const EdgeId ne = add_edge(m, v[2]);
    if (on_new_edge) on_new_edge(ne);
    if (m_edges[e1.i].v[0] == v[1]) std::swap(e1, e2);
    set_face_edges(t, e1, ne, te[2]);
    const FaceId nf = add_face(e2, te[1], ne);
    split_attributes(nf, t);
    split_attributes(t, t);
    join_surfaces(nf, t);
    pieces.insert(pieces.end(), { t, nf });
    cuts.push_back(ne);
  }
  kill_edge(e);

  if ( ! m_split_callouts->callouts.empty())
  {
    std::vector<FaceSplitEvent> fs;
    for (std::size_t i = 0; i < faces.size(); ++i)
      fs.push_back({ faces[i], std::span(pieces).subspan(2 * i, 2), std::span(cuts).subspan(i, 1) });
    call(m_split_callouts, SplitEvent{ es ? std::span(es, 1) : std::span<const EdgeSplitEvent>(), fs });
  }
}

void Geometrium::split_remaining(SplitMemo& memo)
{
  for (Index i = 0; i < Index(memo.m_halves.size()); ++i)
  {
    const auto h = memo.m_halves[i];
    const EdgeId e(i);
    if ( ! h[0].valid() || dead(e)) continue;
    if (dead(h[0]) || dead(h[1]))
      throw std::logic_error("gts::Geometrium::split_remaining: a half of a split edge was removed");
    // h[0] runs from one end of e to the vertex that split it.
    const VertexId m = m_edges[h[0].i].v[1];
    const bool first = m_edges[h[0].i].v[0] == m_edges[e.i].v[0];
    split_faces(e, m, first ? h : std::array{ h[1], h[0] }, {}, nullptr);
  }
}

// GTS kept the halves in the edge's reserved pointer, in the order they were
// made, and put the one from a first when it found them again.
Geometrium::EdgeHalves Geometrium::split_or_find(EdgeId e, VertexId a, VertexId b, const Splitter& split,
                                                 SplitMemo& memo, std::vector<EdgeSplitEvent>& events)
{
  if (memo.m_halves.size() <= e.i) memo.m_halves.resize(m_edges.size());
  auto& h = memo.m_halves[e.i];
  if (h[0].valid())
  {
    if (dead(h[0]) || dead(h[1]))
      throw std::logic_error("gts::Geometrium::split_face4: a half of a split edge was removed");
    EdgeHalves r = { m_edges[h[0].i].v[1], h[0], h[1] };
    if (m_edges[r.from_a.i].v[0] == b) std::swap(r.from_a, r.from_b);
    return r;
  }
  const VertexId v = split(*this, e);
  const EdgeId ea = add_edge(a, v);
  split_attributes(ea, e);
  const EdgeId eb = add_edge(b, v);
  split_attributes(eb, e);
  memo.m_halves[e.i] = { ea, eb };
  events.push_back({ e, v, m_edges[e.i].v[0] == a ? std::array{ ea, eb } : std::array{ eb, ea } });
  return { v, ea, eb };
}

// As tessellate_face in GTS's surface.c.
void Geometrium::split_face4(FaceId f, const Splitter& split, SplitMemo& memo)
{
  const auto [e1, e2, e3] = m_faces[f.i].e;
  const Edge E1 = m_edges[e1.i];
  const Edge E2 = m_edges[e2.i];

  VertexId v1, v2, v3;
  if      (E1.v[1] == E2.v[0]) { v1 = E2.v[1]; v2 = E1.v[0]; v3 = E1.v[1]; }
  else if (E1.v[1] == E2.v[1]) { v1 = E2.v[0]; v2 = E1.v[0]; v3 = E1.v[1]; }
  else if (E1.v[0] == E2.v[0]) { v1 = E2.v[1]; v2 = E1.v[1]; v3 = E1.v[0]; }
  else                         { v1 = E2.v[0]; v2 = E1.v[1]; v3 = E1.v[0]; }

  // GTS removes the face from its edges first and destroys the edges left
  // without faces at the end; an edge is left without faces when f is its
  // only face.
  const bool alone[3] = { next_around(f, e1) == f, next_around(f, e2) == f, next_around(f, e3) == f };

  std::vector<EdgeSplitEvent> edge_events;
  const auto [v4, e24, e34] = split_or_find(e1, v2, v3, split, memo, edge_events);
  const auto [v5, e35, e15] = split_or_find(e2, v3, v1, split, memo, edge_events);
  const auto [v6, e16, e26] = split_or_find(e3, v1, v2, split, memo, edge_events);

  const std::array<EdgeId, 3> cuts = { add_edge(v5, v6), add_edge(v6, v4), add_edge(v4, v5) };
  const auto [e56, e64, e45] = cuts;
  set_face_edges(f, e56, e64, e45);

  const EdgeId old[3] = { e1, e2, e3 };
  for (int k = 0; k < 3; ++k)
    if (alone[k]) kill_edge(old[k]);

  const std::array<FaceId, 4> pieces = { f, add_face(e16, e56, e15), add_face(e26, e24, e64), add_face(e45, e34, e35) };
  for (int k = 1; k < 4; ++k)
  {
    split_attributes(pieces[k], f);
    join_surfaces(pieces[k], f);
  }
  split_attributes(f, f);

  if ( ! m_split_callouts->callouts.empty())
  {
    const FaceSplitEvent fs = { f, pieces, cuts };
    call(m_split_callouts, SplitEvent{ edge_events, std::span(&fs, 1) });
  }
}

} // namespace gts
