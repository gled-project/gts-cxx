// Progressive and hierarchical surfaces, after Hoppe, "Progressive meshes"
// (1996) and "View-dependent refinement of progressive meshes" (1997): a
// surface coarsened by edge collapses that are recorded as vertex splits, so
// that it can be refined again one vertex at a time (PSurface), or in any
// order that the splits allow (HSurface).
//
// A collapse parks the faces on the collapsed edge and the edge's two
// vertices, and an expansion parks the split's vertex; parked elements keep
// their indices and attributes, and come back unchanged. The edges a split
// removes and makes are new each time, as in GTS, with default attributes.
// The surface must not change except through its PSurface or HSurface;
// collect() may be called, and the splits follow the renumbering.
//
// Copyright (C) 1999 Stéphane Popinet (GTS: split.c, psurface.c, hsurface.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-11.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "gts/Coarsen.h"
#include "gts/EHeap.h"
#include "gts/io.h"

#include <expected>
#include <functional>
#include <iosfwd>
#include <memory>
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>

namespace gts {

// GtsSplit: the collapse of the edge between v1 and v2 into v.
struct VertexSplit
{
  VertexId v, v1, v2;
  // The splits whose v is v1 and v2, and the one whose v1 or v2 is v, as
  // indices into the list of splits; -1 for none.
  int child1 = -1, child2 = -1, parent = -1;

  // GtsSplitCFace and CFace: a face on the collapsed edge.
  struct CFace
  {
    FaceId   f;          // invalid until first expanded, after reading
    FaceId   t;          // a face on the edge that replaced f's edges, to find it again
    unsigned flags = 0;  // GTS's CFACE_ flags: orientation, where the edge is in t, edges kept
    std::vector<FaceId> a1, a2; // the faces the collapse moved off f's edges at v1 and v2
  };
  std::vector<CFace> cfaces;
};

// Whether collapse() can collapse vs now: v1 and v2 are joined by an edge
// with faces, which PSurface would collapse (Geometrium::collapsible() and
// no contact), and each of these faces has a neighbour on one of its other
// two edges, which the expansion finds it by. Splits taken out of order, as
// a hierarchical surface takes them, may fail these tests; GTS made no
// test, and failed an assertion half way through the collapse where the
// last one fails.
bool collapsible(const Geometrium& g, const VertexSplit& vs);
// gts_split_collapse: collapses the edge between vs.v1 and vs.v2 into vs.v,
// which must be parked or without edges. on_removed_edge is called with
// each edge removed besides the collapsed one. Throws std::logic_error,
// changing nothing, when vs is not collapsible().
void collapse(Geometrium& g, VertexSplit& vs, const std::function<void(EdgeId)>& on_removed_edge = {});
// gts_split_expand: undoes collapse(); the faces come back into s. Throws
// std::logic_error, before changing anything, when vs is not collapsed.
void expand(Geometrium& g, VertexSplit& vs, Surface& s);
// gts_split_height: the height of the tree of splits below splits[root].
unsigned split_height(const std::vector<VertexSplit>& splits, int root);

class PSurface
{
public:
  // gts_psurface_new: coarsens s as coarsen() does, recording each collapse
  // as a split, and stays at the coarsest level. Edges whose collapse would
  // join two sets of faces only at that edge are not collapsed either.
  PSurface(Surface& s, const CoarsenStop& stop, const EdgeCost& cost = {},
           const CollapsePosition& position = {}, double min_angle = 0.);
  PSurface(const PSurface&) = delete;
  PSurface& operator=(const PSurface&) = delete;

  // gts_psurface_open, gts_psurface_read_vertex and gts_psurface_close: a
  // progressive surface written by write(), or by GTS, read into a new
  // surface of g at its finest level.
  static std::expected<std::unique_ptr<PSurface>, ReadError> read(Geometrium& g, std::istream& in);

  Surface&       surface()       { return *m_s; }
  const Surface& surface() const { return *m_s; }
  // The splits in the order of their collapses, the finest first.
  const std::vector<VertexSplit>& splits() const { return m_splits; }

  // gts_psurface_min_vertex_number, gts_psurface_max_vertex_number and
  // gts_psurface_get_vertex_number.
  std::size_t min_vertex_number() const { return m_min; }
  std::size_t max_vertex_number() const { return m_min + m_splits.size(); }
  std::size_t vertex_number() const     { return m_min + m_splits.size() - m_pos; }

  // gts_psurface_add_vertex and gts_psurface_remove_vertex: expand the next
  // split or collapse the last one; null at the finest or coarsest level.
  const VertexSplit* add_vertex();
  const VertexSplit* remove_vertex();
  // gts_psurface_set_vertex_number: as close to n vertices as there are.
  void set_vertex_number(std::size_t n);

  // gts_psurface_write: collapses s completely and writes it, then expands
  // each split while writing it, so that s ends at its finest level.
  void write(std::ostream& out);

private:
  friend class HSurface;
  PSurface() = default;
  void follow_renumbering();

  Surface*                 m_s = nullptr;
  std::vector<VertexSplit> m_splits;
  std::size_t              m_pos = 0; // the splits before m_pos are collapsed
  std::size_t              m_min = 0;
  Subscription             m_renumber;
};

class HSurface
{
public:
  // The key of a split in the heaps of splits to expand and to collapse,
  // smallest on top. It may change, as with the view; update_keys() takes
  // the new values.
  using SplitKey = std::function<double(const VertexSplit&)>;

  // gts_hsurface_new: takes the splits of ps, which is left without any, as
  // GTS destroyed it, and expands them all.
  HSurface(PSurface& ps, SplitKey expand_key, SplitKey collapse_key);
  HSurface(const HSurface&) = delete;
  HSurface& operator=(const HSurface&) = delete;

  Surface&       surface()       { return *m_s; }
  const Surface& surface() const { return *m_s; }
  const std::vector<VertexSplit>& splits() const { return m_splits; }
  // The splits that are not below another, in the order of GTS's list.
  const std::vector<int>& roots() const { return m_roots; }
  std::size_t vertex_number() const { return m_nvertex; }
  bool expanded(int i) const { return m_expanded[i]; }

  // The split on top of the heap of those that can be expanded, whose
  // vertex is there, and of those that can be collapsed, whose two
  // vertices are there; with its key.
  std::optional<std::pair<int, double>> top_expandable() const  { return m_expandable.top(); }
  std::optional<std::pair<int, double>> top_collapsable() const { return m_collapsable.top(); }
  // gts_eheap_update of both heaps.
  void update_keys();

  // Whether split i is in the heap of those to expand, and whether it is
  // in that of those to collapse and collapsible() now.
  bool can_expand(int i) const   { return ! m_expanded[i] && m_index[i] != kNone; }
  bool can_collapse(int i) const;

  // gts_hsplit_expand: expands split i, which must be expandable. Where
  // faces that it gives back to v1 and v2 are collapsed, it expands their
  // splits first, as force_expand() does; GTS kept them collapsed.
  void expand(int i);
  // gts_hsplit_force_expand: expands split i, and first its parent and the
  // splits of its collapsed neighbours, as far as needed.
  void force_expand(int i);
  // gts_hsplit_collapse: collapses split i, which must be collapsable.
  // Returns false, changing nothing, when the split is not collapsible()
  // now, as it can happen at the boundary of an open surface.
  bool collapse(int i);

  // gts_hsurface_traverse: calls f with each split from the roots, before
  // (pre_order) or after its children, down to depth (all for -1), until f
  // returns true.
  void traverse(bool pre_order, int depth, const std::function<bool(int)>& f) const;
  // gts_hsurface_foreach: expands all splits, then collapses them in their
  // order, calling f before (pre_order) or after each collapse, until f
  // returns true.
  void foreach(bool pre_order, const std::function<bool(int)>& f);
  // gts_hsurface_height
  unsigned height() const;

private:
  using Heap = EHeap<int>;
  static constexpr Heap::Handle kNone = ~Heap::Handle(0);

  void insert(Heap& h, int i) { m_index[i] = h.insert(i); }
  void remove(Heap& h, int i) { h.remove(m_index[i]); m_index[i] = kNone; }
  void record_owner(int i, bool collapsed);

  Surface*                 m_s;
  std::vector<VertexSplit> m_splits;
  std::vector<int>         m_roots;
  std::vector<unsigned>    m_nchild;   // children whose vertex is there, as GtsHSplit's nchild
  std::vector<bool>        m_expanded;
  std::vector<Heap::Handle> m_index;
  std::unordered_map<Index, int> m_owner; // the split holding each collapsed face
  SplitKey                 m_expand_key, m_collapse_key;
  Heap                     m_expandable, m_collapsable;
  std::size_t              m_nvertex = 0;
  Subscription             m_renumber;
};

} // namespace gts
