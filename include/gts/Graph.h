// Graphs, as GTS's GtsGraph: nodes joined by edges, with weights; the dual
// graph of a surface, whose nodes are its faces; progressive graphs, coarsened
// by collapsing the edges of maximal matchings; and partitions of a graph, by
// growing bubbles or by multilevel bisection refined with Kernighan and Lin.
//
// A Graph owns the nodes and edges. A NodeSet is a set of its nodes, as a
// GtsGraph is a container of nodes, and a node may be in several sets: the
// parts of a partition are sets of the nodes of the whole. Weights are
// floats, as in GTS. GTS visited the nodes of a set in the order of a hash
// table keyed by pointer; gts-cxx visits them in the order they were added,
// so partitions differ from GTS's in which of equal choices they take.
//
// Copyright (C) 1999 Stéphane Popinet (GTS: graph.c, pgraph.c, partition.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-11.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "gts/Geometrium.h"
#include "gts/io.h"

#include <array>
#include <deque>
#include <expected>
#include <functional>
#include <iosfwd>
#include <map>
#include <memory>
#include <optional>
#include <ranges>
#include <span>
#include <vector>

namespace gts {

using GNodeId = Id<struct GNodeTag>;
using GEdgeId = Id<struct GEdgeTag>;

//==============================================================================
// Graph and NodeSet
//==============================================================================

class Graph
{
public:
  // gts_gnode_new and gts_wgnode_new. data is the user's, such as the index
  // of a face for the dual graph of a surface (GtsFNode).
  GNodeId add_node(float weight = 1.f, Index data = kNone);
  // gts_gedge_new and gts_wgedge_new: a and b may be the same node.
  GEdgeId add_edge(GNodeId a, GNodeId b, float weight = 1.f, Index data = kNone);

  std::size_t node_capacity() const { return m_nodes.size(); }
  std::size_t edge_capacity() const { return m_edges.size(); }

  // gts_gnode_weight and gts_gedge_weight.
  float weight(GNodeId n) const { return m_nodes[n.i].weight; }
  float weight(GEdgeId e) const { return m_edges[e.i].weight; }
  Index data(GNodeId n) const   { return m_nodes[n.i].data; }
  Index data(GEdgeId e) const   { return m_edges[e.i].data; }
  const std::array<GNodeId, 2>& ends(GEdgeId e) const { return m_edges[e.i].n; }
  bool dead(GEdgeId e) const { return m_edges[e.i].dead; }

  // GTS_GNODE_NEIGHBOR: the other end of e from n.
  GNodeId neighbor(GNodeId n, GEdgeId e) const
  {
    const auto& x = m_edges[e.i].n;
    return x[0] == n ? x[1] : x[1] == n ? x[0] : GNodeId();
  }
  // The edges of n, the newest first, as GTS's list.
  auto edges_of(GNodeId n) const { return m_nodes[n.i].edges | std::views::reverse; }

private:
  friend class NodeSet;
  friend struct GraphAccess;

  struct Node
  {
    float                weight = 1.f;
    Index                data = kNone;
    std::vector<GEdgeId> edges;    // the oldest first
    unsigned             level = 0; // of the last traversal
  };
  struct Edge
  {
    std::array<GNodeId, 2> n;
    float                  weight = 1.f;
    Index                  data = kNone;
    GNodeId                disconnected; // the node a collapse took it from
    bool                   dead = false;
  };

  std::vector<Node> m_nodes;
  std::vector<Edge> m_edges;
};

// GtsGraph and GtsWGraph: a set of the nodes of a graph, in the order they
// were added, and the sum of their weights, kept as nodes come and go.
class NodeSet
{
public:
  explicit NodeSet(const Graph& g) : m_g(&g) {}

  const Graph& graph() const { return *m_g; }
  bool contains(GNodeId n) const { return n.i < m_pos.size() && m_pos[n.i] != 0; }
  // gts_container_add and gts_container_remove; adding a node of the set or
  // removing one that is not does nothing.
  void add(GNodeId n);
  void remove(GNodeId n);
  std::size_t size() const { return m_size; }
  // gts_graph_weight
  float weight() const { return m_weight; }
  // The nodes, in the order they were added.
  std::vector<GNodeId> nodes() const;

private:
  const Graph*         m_g;
  std::vector<GNodeId> m_order;   // with the removed ones, until compacted
  std::vector<Index>   m_pos;     // 1 + the place in m_order, by node index; 0 for out
  std::size_t          m_size = 0;
  float                m_weight = 0.f;
};

// gts_gnode_degree: the edges of n to nodes of g, or all of them for null.
unsigned degree(const Graph& gr, GNodeId n, const NodeSet* g = nullptr);
// gts_gnode_move_cost: the weight of n's edges into src, less that of its
// edges into dst, which moving n from src to dst saves.
float move_cost(GNodeId n, const NodeSet& src, const NodeSet& dst);
// gts_graph_foreach_edge: each edge of a node of g once, also those to
// nodes outside g, in the order of the nodes and of their lists.
std::vector<GEdgeId> edges_of(const NodeSet& g);
// gts_graph_edges_cut and gts_graph_edges_cut_weight: the edges from
// nodes of g to nodes outside it.
unsigned edges_cut(const NodeSet& g);
float    edges_cut_weight(const NodeSet& g);
// gts_wgraph_weight_max: the largest weight of a node of g.
float    weight_max(const NodeSet& g);
// gts_graph_distance_sum: the sum of the distances, in edges, from center
// to the nodes of g it reaches.
unsigned distance_sum(const NodeSet& g, GNodeId center);
// gts_graph_farthest: the node of g farthest from all of seeds, by
// breadth-first traversals from each.
GNodeId  farthest(const NodeSet& g, std::span<const GNodeId> seeds);

// gts_graph_traverse_*: the nodes of g in breadth-first order from n, each
// level 1 + its distance; with reinit, the levels of the other nodes of g
// are reset first. Traversals may share levels, as farthest()'s do.
class GraphTraverse
{
public:
  GraphTraverse(const NodeSet& g, GNodeId n, bool reinit = true);
  // gts_graph_traverse_next and gts_graph_traverse_what_next.
  GNodeId next();
  GNodeId what_next() const { return m_queue.empty() ? GNodeId() : m_queue.front(); }
  unsigned level(GNodeId n) const;

private:
  const NodeSet&      m_g;
  std::deque<GNodeId> m_queue;
};

// gts_graph_print_stats
void print_stats(const NodeSet& g, std::ostream& out);
// gts_graph_write: "nn ne", the nodes, then each edge "i j", from 1, as
// GTS's GtsGraph with plain nodes and edges.
void write_graph(const NodeSet& g, std::ostream& out);
// gts_graph_write_dot
void write_dot(const NodeSet& g, std::ostream& out);
// gts_graph_read: a graph written by write_graph(), into gr, as a new set.
// As GTS, it reads the lines of plain nodes wrongly; see QUIRKS.md.
std::expected<NodeSet, ReadError> read_graph(Graph& gr, std::istream& in);
// gts_graph_read_jostle: a graph in the format of Jostle and the Graph
// Partitioning Archive; each node's data is its line number. As GTS, it
// loops for ever on a graph with an edge; see QUIRKS.md.
std::expected<NodeSet, ReadError> read_jostle(Graph& gr, std::istream& in);

// gts_surface_graph_new: the dual graph of s, a node for each face, with
// the face's index as data, and an edge for each pair of faces of s on an
// edge, with the edge's index as data.
NodeSet surface_graph(Graph& gr, const Surface& s);
// gts_segments_graph_new: a node for each vertex of edges, with the
// vertex's index as data, and an edge for each edge.
NodeSet segments_graph(Graph& gr, const Geometrium& geo, std::span<const EdgeId> edges);
// gts_surface_graph_surface: a new surface of geo with the faces of the
// nodes of g, those of a surface graph.
Surface& graph_surface(const NodeSet& g, Geometrium& geo);

//==============================================================================
// Progressive graphs, from pgraph.c
//==============================================================================

// GtsGNodeSplit: the collapse of n1 and n2 into n.
struct GNodeSplit
{
  GNodeId n, n1, n2;
};

// gts_gnode_split_collapse: replaces n1 and n2 by n in g. An edge from n1
// or n2 to a node joined to both gives way to one edge from n with their
// summed weight.
void collapse(Graph& gr, const GNodeSplit& ns, NodeSet& g);
// gts_gnode_split_expand: undoes collapse().
void expand(Graph& gr, const GNodeSplit& ns, NodeSet& g);

// GtsPGraph: g coarsened by collapsing the edges of maximal matchings, each
// matching a level, down to min nodes.
class PGraph
{
public:
  PGraph(Graph& gr, NodeSet& g, std::size_t min);

  const std::vector<GNodeSplit>& splits() const { return m_splits; }
  // gts_pgraph_add_node and gts_pgraph_remove_node: the split expanded or
  // collapsed; null at the end.
  const GNodeSplit* add_node();
  const GNodeSplit* remove_node();
  // gts_pgraph_set_node_number, get, min and max.
  void set_node_number(std::size_t n);
  std::size_t node_number() const     { return m_min + m_splits.size() - m_pos; }
  std::size_t min_node_number() const { return m_min; }
  std::size_t max_node_number() const { return m_min + m_splits.size(); }
  // The number of levels, matchings, from the finest.
  std::size_t levels() const { return m_levels.size(); }
  // gts_pgraph_down: expands the splits of one level, calling f with each;
  // false when at the finest level.
  bool down(const std::function<void(const GNodeSplit&)>& f = {});

private:
  Graph*                  m_gr;
  NodeSet*                m_g;
  std::vector<GNodeSplit> m_splits;
  std::vector<std::size_t> m_levels; // the size of g before each level
  std::size_t             m_pos = 0, m_min = 0, m_level = 0;
};

//==============================================================================
// Partitions, from partition.c
//==============================================================================

// gts_graph_partition_edges_cut, gts_graph_partition_edges_cut_weight and
// gts_graph_partition_balance: the edges between the parts, and the
// difference between the heaviest and the lightest part.
unsigned partition_edges_cut(std::span<const NodeSet> parts);
float    partition_edges_cut_weight(std::span<const NodeSet> parts);
float    partition_balance(std::span<const NodeSet> parts);
// gts_graph_partition_print_stats
void print_partition_stats(std::span<const NodeSet> parts, std::ostream& out);

// gts_graph_bubble_partition: np parts grown in breadth-first bubbles from
// seeds as far apart as can be, then niter times from the centres of the
// bubbles. step, when given, is called after each iteration.
std::vector<NodeSet> bubble_partition(const NodeSet& g, unsigned np, unsigned niter,
                                      const std::function<void(std::span<const NodeSet>)>& step = {});

// GtsGraphBisection: g in two, g1 and g2, and the nodes of each with a
// neighbour in the other, mapped to that neighbour.
struct GraphBisection
{
  const NodeSet* g;
  NodeSet g1, g2;
  std::map<Index, GNodeId> bg1, bg2;
};

// gts_graph_bisection_check: whether bg1 and bg2 are the boundary nodes.
bool check(const GraphBisection& bg);
// gts_graph_ggg_bisection: grows g1 from ntry seeds of the smallest degree
// by the greedy graph growing of Karypis and Kumar, up to half the weight,
// and keeps the cheapest balanced one.
GraphBisection ggg_bisection(const NodeSet& g, unsigned ntry);
// gts_graph_bfgg_bisection: the same, growing breadth-first.
GraphBisection bfgg_bisection(const NodeSet& g, unsigned ntry);
// gts_graph_bisection_kl_refine: Kernighan-Lin moves from the heavier part,
// stopping after mmax moves that gain nothing, which are undone. Returns
// the gain, negative for a better cut.
double kl_refine(GraphBisection& bg, unsigned mmax);
// gts_graph_bisection_bkl_refine: the same over the boundary nodes, within
// imbalance, a fraction of the weight of g.
double bkl_refine(GraphBisection& bg, unsigned mmax, float imbalance);
// gts_graph_bisection_new: the multilevel bisection of g: coarsened to nmin
// nodes, bisected by ggg_bisection(), then refined with bkl_refine() at
// each level on the way back.
GraphBisection bisection(Graph& gr, NodeSet& g, unsigned ntry, unsigned mmax, unsigned nmin, float imbalance);
// gts_graph_recursive_bisection: 2^n parts by bisection() of the parts,
// in the order of GTS's list.
std::vector<NodeSet> recursive_bisection(Graph& gr, NodeSet& g, unsigned n, unsigned ntry, unsigned mmax,
                                         unsigned nmin, float imbalance);

} // namespace gts
