// Copyright (C) 1999-2003 Wagner Toledo Correa, Stéphane Popinet (GTS: stripe.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-11.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/Strip.h"
#include "gts/EHeap.h"
#include "gts/Topology.h"

#include <algorithm>
#include <cfloat>
#include <deque>
#include <stdexcept>

namespace gts {

namespace {

using Heap = EHeap<FaceId>;
constexpr Heap::Handle kNoHandle = ~Heap::Handle(0);

// tri_data_t, by face index.
struct TriData
{
  bool                used = false;
  std::vector<FaceId> neighbors;
  Heap::Handle        pos = kNoHandle;
};

struct Strips
{
  const Geometrium&    g;
  std::vector<TriData> td;
  Heap                 heap;

  explicit Strips(const Surface& s) : g(s.geometrium()), td(g.face_capacity())
  {
    // gts_triangle_neighbors: of any surface. GTS failed an assertion on
    // one that is not in s; see QUIRKS.md.
    for (FaceId f : s.faces())
    {
      td[f.i].neighbors = face_neighbors(g, f);
      for (FaceId t : td[f.i].neighbors)
        if ( ! s.contains(t)) throw std::logic_error("gts::strips: a neighbouring face is not in the surface");
    }
    for (FaceId f : s.faces()) td[f.i].pos = heap.insert(f, double(unused_neighbors2(f).size()));
  }

  // tri_data_unused_neighbors2: the unused neighbours of f and their unused
  // neighbours other than f, each once.
  std::vector<FaceId> unused_neighbors2(FaceId f) const
  {
    std::vector<FaceId> h;
    auto add = [&](FaceId t) { if (std::find(h.begin(), h.end(), t) == h.end()) h.push_back(t); };
    for (FaceId t2 : td[f.i].neighbors)
      if ( ! td[t2.i].used)
      {
        add(t2);
        for (FaceId t3 : td[t2.i].neighbors)
          if (t3 != f && ! td[t3.i].used) add(t3);
      }
    return h;
  }

  // heap_remove: marks f used and lowers the keys of the faces near it.
  void remove(FaceId f)
  {
    TriData& d = td[f.i];
    d.used = true;
    heap.remove(d.pos);
    d.pos = kNoHandle;
    for (FaceId t : unused_neighbors2(f))
    {
      const double k = double(unused_neighbors2(t).size());
      if (k != heap.key(td[t.i].pos)) heap.decrease_key(td[t.i].pos, k);
    }
  }

  FaceId find_min_neighbor(FaceId f) const
  {
    FaceId min_neighbor;
    double min_key = DBL_MAX;
    for (FaceId t2 : td[f.i].neighbors)
    {
      if (td[t2.i].used) continue;
      const double k = heap.key(td[t2.i].pos);
      if (k < min_key)
      {
        min_key = k;
        min_neighbor = t2;
      }
    }
    return min_neighbor;
  }

  using V3 = std::array<VertexId, 3>;

  // vertices_match: rotates w until its vertices match the given ones,
  // where given; invalid ones match anything.
  static bool vertices_match(VertexId a, VertexId b, VertexId c, V3& w)
  {
    auto match = [&] { return ( ! a.valid() || a == w[0]) && ( ! b.valid() || b == w[1]) && ( ! c.valid() || c == w[2]); };
    for (int i = 0; i < 2; i++)
    {
      if (match()) return true;
      std::rotate(w.begin(), w.begin() + 1, w.end());
    }
    return match();
  }

  FaceId find_neighbor_forward(FaceId f, V3& v, bool left_turn) const
  {
    for (FaceId t2 : td[f.i].neighbors)
    {
      if (t2 == f || td[t2.i].used) continue;
      V3 w = g.vertices_of(t2);
      if (left_turn ? ! vertices_match(v[0], v[2], VertexId(), w) : ! vertices_match(v[2], v[1], VertexId(), w))
        continue;
      v = w;
      return t2;
    }
    return {};
  }

  FaceId find_neighbor_backward(FaceId f, V3& v, bool left_turn) const
  {
    for (FaceId t2 : td[f.i].neighbors)
    {
      if (t2 == f || td[t2.i].used) continue;
      V3 w = g.vertices_of(t2);
      if (left_turn ? ! vertices_match(VertexId(), v[1], v[0], w) : ! vertices_match(v[0], VertexId(), v[1], w))
        continue;
      v = w;
      return t2;
    }
    return {};
  }

  static bool one_of(VertexId x, const V3& v) { return x == v[0] || x == v[1] || x == v[2]; }

  // match_vertex: rotates v until x is first.
  static void match_vertex(VertexId x, V3& v)
  {
    while (v[0] != x) std::rotate(v.begin(), v.begin() + 1, v.end());
  }

  // find_right_turn: v and w, which share an edge, rotated so that v ends
  // at w's first vertex; whether w turns right from v.
  static bool find_right_turn(V3& v, V3& w)
  {
    VertexId u;
    for (VertexId x : v)
      if ( ! one_of(x, w)) { u = x; break; }
    if ( ! u.valid()) return false;
    match_vertex(u, v);
    if ( ! one_of(v[2], w)) return false;
    match_vertex(v[2], w);
    return w[1] == v[1];
  }
};

} // namespace

std::vector<std::vector<FaceId>> strips(const Surface& s)
{
  Strips h(s);
  std::vector<std::vector<FaceId>> result;
  while ( ! h.heap.empty())
  {
    const FaceId t1 = h.heap.top()->first;
    h.remove(t1);
    // GTS's list, front first; built by prepending, then reversed.
    std::deque<FaceId> strip{ t1 };

    if (const FaceId t2 = h.find_min_neighbor(t1); t2.valid())
    {
      Strips::V3 v = h.g.vertices_of(t1), w = h.g.vertices_of(t2);
      if (Strips::find_right_turn(v, w))
      {
        h.remove(t2);
        strip.push_front(t2);

        // grow_strip_forward
        bool left_turn = true;
        FaceId t = t2;
        while ((t = h.find_neighbor_forward(t, w, left_turn)).valid())
        {
          h.remove(t);
          strip.push_front(t);
          left_turn = ! left_turn;
        }
        std::reverse(strip.begin(), strip.end());

        // grow_strip_backward, two faces at a time, so that the turns
        // keep alternating.
        t = t1;
        FaceId ta;
        while ((ta = h.find_neighbor_backward(t, v, false)).valid() &&
               (t = h.find_neighbor_backward(ta, v, true)).valid())
        {
          h.remove(ta);
          h.remove(t);
          strip.push_front(ta);
          strip.push_front(t);
        }
      }
    }
    result.emplace_back(strip.begin(), strip.end());
  }
  return result;
}

} // namespace gts
