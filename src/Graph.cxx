// Copyright (C) 1999 Stéphane Popinet (GTS: graph.c, pgraph.c, partition.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-11.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/Graph.h"
#include "gts/EHeap.h"
#include "gts/Stats.h"
#include "io_detail.h"

#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <format>
#include <istream>
#include <iterator>
#include <ostream>
#include <stdexcept>

namespace gts {

// The changes to nodes and edges that progressive graphs make.
struct GraphAccess
{
  static Graph::Node& node(Graph& g, GNodeId n) { return g.m_nodes[n.i]; }
  static Graph::Edge& edge(Graph& g, GEdgeId e) { return g.m_edges[e.i]; }
  static const Graph::Node& node(const Graph& g, GNodeId n) { return g.m_nodes[n.i]; }

  // g_slist_remove: e out of n's list, the newest occurrence.
  static void unlist(Graph& g, GNodeId n, GEdgeId e)
  {
    auto& v = g.m_nodes[n.i].edges;
    for (auto it = v.rbegin(); it != v.rend(); ++it)
      if (*it == e) { v.erase(std::next(it).base()); return; }
  }
  // gts_container_add on a node: e newest in n's list.
  static void list(Graph& g, GNodeId n, GEdgeId e) { g.m_nodes[n.i].edges.push_back(e); }
  static unsigned& level(const Graph& g, GNodeId n) { return const_cast<Graph&>(g).m_nodes[n.i].level; }
};

//==============================================================================
// Graph and NodeSet
//==============================================================================

GNodeId Graph::add_node(float weight, Index data)
{
  m_nodes.push_back({ weight, data, {}, 0 });
  return GNodeId(Index(m_nodes.size() - 1));
}

GEdgeId Graph::add_edge(GNodeId a, GNodeId b, float weight, Index data)
{
  const GEdgeId e(Index(m_edges.size()));
  m_edges.push_back({ { a, b }, weight, data, {}, false });
  m_nodes[a.i].edges.push_back(e);
  if (a != b) m_nodes[b.i].edges.push_back(e);
  return e;
}

void NodeSet::add(GNodeId n)
{
  if (contains(n)) return;
  if (m_pos.size() <= n.i) m_pos.resize(n.i + 1, 0);
  m_order.push_back(n);
  m_pos[n.i] = Index(m_order.size());
  ++m_size;
  m_weight += m_g->weight(n);
}

void NodeSet::remove(GNodeId n)
{
  if ( ! contains(n)) return;
  m_pos[n.i] = 0;
  --m_size;
  m_weight -= m_g->weight(n);
  if (m_order.size() > 2 * m_size + 16)
  {
    std::vector<GNodeId> order;
    for (std::size_t i = 0; i < m_order.size(); ++i)
      if (m_pos[m_order[i].i] == i + 1) order.push_back(m_order[i]);
    m_order.swap(order);
    for (std::size_t i = 0; i < m_order.size(); ++i) m_pos[m_order[i].i] = Index(i + 1);
  }
}

std::vector<GNodeId> NodeSet::nodes() const
{
  std::vector<GNodeId> r;
  r.reserve(m_size);
  for (std::size_t i = 0; i < m_order.size(); ++i)
    if (m_pos[m_order[i].i] == i + 1) r.push_back(m_order[i]);
  return r;
}

namespace {

bool in(const NodeSet* g, GNodeId n) { return n.valid() && ( ! g || g->contains(n)); }

} // namespace

unsigned degree(const Graph& gr, GNodeId n, const NodeSet* g)
{
  unsigned nn = 0;
  for (GEdgeId e : gr.edges_of(n))
    if (in(g, gr.neighbor(n, e))) ++nn;
  return nn;
}

float move_cost(GNodeId n, const NodeSet& src, const NodeSet& dst)
{
  const Graph& gr = src.graph();
  float cost = 0.f;
  for (GEdgeId e : gr.edges_of(n))
  {
    const GNodeId m = gr.neighbor(n, e);
    if (in(&src, m))      cost += gr.weight(e);
    else if (in(&dst, m)) cost -= gr.weight(e);
  }
  return cost;
}

std::vector<GEdgeId> edges_of(const NodeSet& g)
{
  const Graph& gr = g.graph();
  std::vector<bool> seen(gr.edge_capacity());
  std::vector<GEdgeId> r;
  for (GNodeId n : g.nodes())
    for (GEdgeId e : gr.edges_of(n))
      if ( ! seen[e.i]) { seen[e.i] = true; r.push_back(e); }
  return r;
}

unsigned edges_cut(const NodeSet& g)
{
  const Graph& gr = g.graph();
  unsigned cuts = 0;
  for (GNodeId n : g.nodes())
    for (GEdgeId e : gr.edges_of(n))
      if (const GNodeId m = gr.neighbor(n, e); m.valid() && ! g.contains(m)) ++cuts;
  return cuts;
}

float edges_cut_weight(const NodeSet& g)
{
  const Graph& gr = g.graph();
  float weight = 0.f;
  for (GNodeId n : g.nodes())
    for (GEdgeId e : gr.edges_of(n))
      if (const GNodeId m = gr.neighbor(n, e); m.valid() && ! g.contains(m)) weight += gr.weight(e);
  return weight;
}

float weight_max(const NodeSet& g)
{
  float wmax = -FLT_MAX;
  for (GNodeId n : g.nodes()) wmax = std::max(wmax, g.graph().weight(n));
  return wmax;
}

GraphTraverse::GraphTraverse(const NodeSet& g, GNodeId n, bool reinit) : m_g(g)
{
  if ( ! g.contains(n)) throw std::invalid_argument("gts::GraphTraverse: the node is not in the set");
  if (reinit)
    for (GNodeId m : g.nodes()) GraphAccess::level(g.graph(), m) = 0;
  GraphAccess::level(g.graph(), n) = 1;
  m_queue.push_back(n);
}

GNodeId GraphTraverse::next()
{
  if (m_queue.empty()) return {};
  const GNodeId u = m_queue.front();
  m_queue.pop_front();
  const Graph& gr = m_g.graph();
  for (GEdgeId e : gr.edges_of(u))
  {
    const GNodeId n = gr.neighbor(u, e);
    if ( ! in(&m_g, n)) continue;
    unsigned& l = GraphAccess::level(gr, n);
    if (l == 0)
    {
      l = GraphAccess::level(gr, u) + 1;
      m_queue.push_back(n);
    }
  }
  return u;
}

unsigned GraphTraverse::level(GNodeId n) const { return GraphAccess::level(m_g.graph(), n); }

unsigned distance_sum(const NodeSet& g, GNodeId center)
{
  GraphTraverse t(g, center, true);
  unsigned sum = 0;
  for (GNodeId n = t.next(); n.valid(); n = t.next()) sum += t.level(n) - 1;
  return sum;
}

GNodeId farthest(const NodeSet& g, std::span<const GNodeId> seeds)
{
  std::vector<GraphTraverse> ts;
  ts.reserve(seeds.size());
  bool reinit = true;
  for (GNodeId s : seeds)
  {
    ts.emplace_back(g, s, reinit);
    reinit = false;
  }
  GNodeId far;
  bool changed = true;
  unsigned level = 1;
  while (changed)
  {
    changed = false;
    for (auto& t : ts)
    {
      GNodeId n;
      while ((n = t.what_next()).valid() && t.level(n) == level)
      {
        changed = true;
        far = n;
        t.next();
      }
    }
    ++level;
  }
  return far;
}

void print_stats(const NodeSet& g, std::ostream& out)
{
  out << std::format("# nodes: {} weight: {:g}\n", g.size(), g.weight());
  out << "#   degree: ";
  Range degree;
  for (GNodeId n : g.nodes()) degree.add(gts::degree(g.graph(), n, &g));
  degree.update();
  out << to_string(degree) << '\n';
  out << std::format("#   edges cut: {} edges cut weight: {:g}\n", edges_cut(g), edges_cut_weight(g));
}

void write_graph(const NodeSet& g, std::ostream& out)
{
  const Graph& gr = g.graph();
  const std::vector<GNodeId> nodes = g.nodes();
  std::vector<Index> num(gr.node_capacity(), 0);
  for (std::size_t i = 0; i < nodes.size(); ++i) num[nodes[i].i] = Index(i + 1);
  std::vector<GEdgeId> edges;
  for (GEdgeId e : edges_of(g))
    if (g.contains(gr.ends(e)[0]) && g.contains(gr.ends(e)[1])) edges.push_back(e);
  out << std::format("{} {} GtsGraph GtsGNode GtsGEdge\n", nodes.size(), edges.size());
  for (std::size_t i = 0; i < nodes.size(); ++i) out << '\n';
  for (GEdgeId e : edges) out << std::format("{} {}\n", num[gr.ends(e)[0].i], num[gr.ends(e)[1].i]);
}

void write_dot(const NodeSet& g, std::ostream& out)
{
  const Graph& gr = g.graph();
  const std::vector<GNodeId> nodes = g.nodes();
  std::vector<Index> num(gr.node_capacity(), 0);
  out << "digraph \"gts\" {\n";
  for (std::size_t i = 0; i < nodes.size(); ++i)
  {
    num[nodes[i].i] = Index(i + 1);
    out << std::format("  n{}", i + 1);
    if (gr.data(nodes[i]) != kNone) out << std::format(" [label=\"{}\",]", gr.data(nodes[i]));
    out << ";\n";
  }
  for (GEdgeId e : edges_of(g))
  {
    const auto& n = gr.ends(e);
    if ( ! g.contains(n[0]) || ! g.contains(n[1])) continue;
    out << std::format("  n{} -> n{}", num[n[0].i], num[n[1].i]);
    if (gr.data(e) != kNone) out << std::format(" [label=\"{}\",]", gr.data(e));
    out << ";\n";
  }
  out << "}\n";
}

namespace {

// The rest of the line and its newline, but no blank line after it.
void end_line(detail::Tokenizer& f)
{
  while ( ! f.failed() && f.type() != detail::Tokenizer::None && f.ch() != '\n') f.next();
  if (f.ch() == '\n') f.next();
}

} // namespace

// GTS skipped every blank line after the header, which took the empty lines
// of its own nodes for one, and then read edges as nodes; gts-cxx reads a
// line for each node. See QUIRKS.md.
std::expected<NodeSet, ReadError> read_graph(Graph& gr, std::istream& in)
{
  using T = detail::Tokenizer;
  T f(std::string(std::istreambuf_iterator<char>(in), {}));
  while (f.ch() == '\n') f.next();
  auto count = [&](const char* what) -> long {
    if (f.type() != T::Int) { f.fail(std::format("expecting an integer ({})", what)); return 0; }
    const long n = std::atol(f.token().c_str());
    f.next();
    return n;
  };
  const long nn = count("number of nodes");
  const long ne = f.failed() ? 0 : count("number of edges");
  if (f.failed()) return std::unexpected(f.error());
  end_line(f);

  NodeSet g(gr);
  std::vector<GNodeId> nodes;
  for (long n = 0; n < nn && ! f.failed(); ++n)
  {
    nodes.push_back(gr.add_node());
    g.add(nodes.back());
    end_line(f);
  }
  for (long n = 0; n < ne && ! f.failed(); ++n)
  {
    long idx[2];
    for (int k = 0; k < 2 && ! f.failed(); ++k)
    {
      if (f.type() != T::Int) { f.fail(std::format("expecting an integer ({} node index)", k ? "second" : "first")); break; }
      idx[k] = std::atol(f.token().c_str());
      if (idx[k] <= 0 || idx[k] > long(nodes.size()))
        f.fail(std::format("node index `{}' is out of range `[1,{}]'", idx[k], nodes.size()));
      f.next();
    }
    if (f.failed()) break;
    gr.add_edge(nodes[idx[0] - 1], nodes[idx[1] - 1]);
    end_line(f);
  }
  if (f.failed()) return std::unexpected(f.error());
  return g;
}

std::expected<NodeSet, ReadError> read_jostle(Graph& gr, std::istream& in)
{
  using T = detail::Tokenizer;
  T f(std::string(std::istreambuf_iterator<char>(in), {}));
  while (f.ch() == '\n') f.next();
  if (f.type() != T::Int) { f.fail("expecting an integer (number of nodes)"); return std::unexpected(f.error()); }
  const long nn = std::atol(f.token().c_str());
  f.next();
  if (f.type() != T::Int) { f.fail("expecting an integer (number of edges)"); return std::unexpected(f.error()); }
  long ne = std::atol(f.token().c_str());
  f.skip_line();

  NodeSet g(gr);
  std::vector<GNodeId> nodes;
  long n = 0;
  // GTS did not move past a neighbour numbered after the node, and looped
  // for ever there, and took a node without neighbours for an error; see
  // QUIRKS.md. A line may list no neighbours.
  while (n < nn && ! f.failed() && f.type() != T::None)
  {
    // GtsNGNode: the line of the node is its id.
    const GNodeId node = gr.add_node(1.f, f.line());
    nodes.push_back(node);
    ++n;
    g.add(node);
    while ( ! f.failed() && f.type() == T::Int)
    {
      const long i = std::atol(f.token().c_str());
      if (i <= 0 || i > nn) f.fail(std::format("node index `{}' is out of range `[1,{}]'", i, nn));
      else if (i == n)      f.fail(std::format("node index `{}' references itself", i));
      else
      {
        if (i < n)
        {
          gr.add_edge(node, nodes[i - 1]);
          --ne;
        }
        f.next();
      }
    }
    if ( ! f.failed() && f.type() != T::None && f.ch() != '\n') f.fail("expecting an integer (node index)");
    if (f.ch() == '\n') f.next();
  }
  if ( ! f.failed())
  {
    if (n != nn)     f.fail(std::format("only `{}' nodes read out of `{}'", n, nn));
    else if (ne > 0) f.fail(std::format("`{}' unallocated edges remaining", ne));
  }
  if (f.failed()) return std::unexpected(f.error());
  return g;
}

NodeSet surface_graph(Graph& gr, const Surface& s)
{
  const Geometrium& geo = s.geometrium();
  NodeSet g(gr);
  std::vector<GNodeId> node(geo.face_capacity());
  for (FaceId f : s.faces())
  {
    node[f.i] = gr.add_node(1.f, f.i);
    g.add(node[f.i]);
  }
  for (EdgeId e : s.edges())
  {
    std::vector<FaceId> faces;
    for (FaceId f : geo.faces_of(e)) if (s.contains(f)) faces.push_back(f);
    for (std::size_t i = 0; i < faces.size(); ++i)
      for (std::size_t j = i + 1; j < faces.size(); ++j)
        gr.add_edge(node[faces[i].i], node[faces[j].i], 1.f, e.i);
  }
  return g;
}

NodeSet segments_graph(Graph& gr, const Geometrium& geo, std::span<const EdgeId> edges)
{
  NodeSet g(gr);
  std::vector<GNodeId> node(geo.vertex_capacity());
  for (EdgeId e : edges)
  {
    const auto& v = geo.edge(e).v;
    for (VertexId x : v)
      if ( ! node[x.i].valid())
      {
        node[x.i] = gr.add_node(1.f, x.i);
        g.add(node[x.i]);
      }
    gr.add_edge(node[v[0].i], node[v[1].i], 1.f, e.i);
  }
  return g;
}

Surface& graph_surface(const NodeSet& g, Geometrium& geo)
{
  Surface& s = geo.new_surface();
  for (GNodeId n : g.nodes())
    if (g.graph().data(n) != kNone) s.add_face(FaceId(g.graph().data(n)));
  return s;
}

//==============================================================================
// Progressive graphs, from pgraph.c
//==============================================================================

void collapse(Graph& gr, const GNodeSplit& ns, NodeSet& g)
{
  using A = GraphAccess;
  const GNodeId n = ns.n, n1 = ns.n1, n2 = ns.n2;
  if ( ! gr.edges_of(n).empty()) throw std::logic_error("gts::collapse: the node of the split has edges");

  // Triangles n1 (e13) n3 (e32) n2: one edge from n to n3 takes the place
  // of both, which are disconnected from n3.
  std::vector<GEdgeId> items1(gr.edges_of(n1).begin(), gr.edges_of(n1).end());
  for (GEdgeId e13 : items1)
  {
    const GNodeId n3 = gr.neighbor(n1, e13);
    if (n3 == n2 || ! n3.valid()) continue;
    std::vector<GEdgeId> items3(gr.edges_of(n3).begin(), gr.edges_of(n3).end());
    for (GEdgeId e32 : items3)
      if (gr.neighbor(n3, e32) == n2)
      {
        gr.add_edge(n, n3, gr.weight(e13) + gr.weight(e32));
        A::edge(gr, e13).disconnected = n3;
        A::edge(gr, e32).disconnected = n3;
        A::unlist(gr, n3, e32);
      }
    if (A::edge(gr, e13).disconnected == n3) A::unlist(gr, n3, e13);
  }

  // The other edges move to n.
  for (GNodeId x : { n1, n2 })
  {
    std::vector<GEdgeId> items(gr.edges_of(x).begin(), gr.edges_of(x).end());
    for (GEdgeId e : items)
    {
      auto& E = A::edge(gr, e);
      const bool connects = (E.n[0] == n1 && E.n[1] == n2) || (E.n[0] == n2 && E.n[1] == n1);
      if (E.disconnected.valid() || connects) continue;
      if (E.n[0] == n1 || E.n[0] == n2)      E.n[0] = n;
      else if (E.n[1] == n1 || E.n[1] == n2) E.n[1] = n;
      else throw std::logic_error("gts::collapse: an edge of the node is not on it");
      A::list(gr, n, e);
    }
  }
  g.remove(n1);
  g.remove(n2);
  g.add(n);
}

void expand(Graph& gr, const GNodeSplit& ns, NodeSet& g)
{
  using A = GraphAccess;
  const GNodeId n = ns.n;
  if ( ! g.contains(n)) throw std::logic_error("gts::expand: the node of the split is not in the set");

  // restore_edge
  for (auto [x, y] : { std::pair{ ns.n1, ns.n2 }, std::pair{ ns.n2, ns.n1 } })
  {
    std::vector<GEdgeId> items(gr.edges_of(x).begin(), gr.edges_of(x).end());
    for (GEdgeId e : items)
    {
      auto& E = A::edge(gr, e);
      if (const GNodeId n3 = E.disconnected; n3.valid())
      {
        E.disconnected = GNodeId();
        A::list(gr, n3, e);
        continue;
      }
      if ((E.n[0] == x && E.n[1] == y) || (E.n[0] == y && E.n[1] == x)) continue;
      if (E.n[0] == n)      E.n[0] = x;
      else if (E.n[1] == n) E.n[1] = x;
      else throw std::logic_error("gts::expand: an edge of the node is not on it");
      A::unlist(gr, n, e);
    }
  }
  // The edges of n left are those that replaced triangles.
  std::vector<GEdgeId> items(gr.edges_of(n).begin(), gr.edges_of(n).end());
  for (GEdgeId e : items)
  {
    auto& E = A::edge(gr, e);
    A::unlist(gr, E.n[0], e);
    if (E.n[1] != E.n[0]) A::unlist(gr, E.n[1], e);
    E.dead = true;
  }
  g.remove(n);
  g.add(ns.n1);
  g.add(ns.n2);
}

namespace {

// maximal_matching: for each node of g not yet matched, its heaviest edge
// to a node of g not yet matched; in the order of GTS's list.
std::vector<GEdgeId> maximal_matching(const Graph& gr, const NodeSet& g)
{
  std::vector<GEdgeId> list;
  std::vector<bool> matched(gr.node_capacity());
  for (GNodeId n : g.nodes())
  {
    if (matched[n.i]) continue;
    float wmax = -FLT_MAX;
    GEdgeId emax;
    for (GEdgeId e : gr.edges_of(n))
    {
      const GNodeId n1 = gr.neighbor(n, e);
      if (n1.valid() && ! matched[n1.i] && gr.weight(e) > wmax && g.contains(n1))
      {
        emax = e;
        wmax = gr.weight(e);
      }
    }
    if (emax.valid())
    {
      matched[gr.neighbor(n, emax).i] = true;
      matched[n.i] = true;
      list.push_back(emax);
    }
  }
  std::reverse(list.begin(), list.end());
  return list;
}

} // namespace

PGraph::PGraph(Graph& gr, NodeSet& g, std::size_t min) : m_gr(&gr), m_g(&g)
{
  std::vector<GEdgeId> matching;
  while (g.size() > min && ! (matching = maximal_matching(gr, g)).empty())
  {
    m_levels.push_back(g.size());
    for (GEdgeId e : matching)
    {
      if (g.size() <= min) break;
      const auto [a, b] = gr.ends(e);
      const GNodeSplit ns{ gr.add_node(gr.weight(a) + gr.weight(b)), a, b };
      collapse(gr, ns, g);
      m_splits.push_back(ns);
    }
  }
  m_pos = m_splits.size();
  m_min = g.size();
  m_level = m_levels.size();
}

const GNodeSplit* PGraph::add_node()
{
  if (m_pos == 0) return nullptr;
  const GNodeSplit& ns = m_splits[--m_pos];
  expand(*m_gr, ns, *m_g);
  return &ns;
}

const GNodeSplit* PGraph::remove_node()
{
  if (m_pos == m_splits.size()) return nullptr;
  const GNodeSplit& ns = m_splits[m_pos++];
  collapse(*m_gr, ns, *m_g);
  return &ns;
}

void PGraph::set_node_number(std::size_t n)
{
  n = std::clamp(n, min_node_number(), max_node_number());
  const std::size_t pos = m_min + m_splits.size() - n;
  while (m_pos > pos && add_node()) {}
  while (m_pos < pos && remove_node()) {}
}

bool PGraph::down(const std::function<void(const GNodeSplit&)>& f)
{
  if (m_level == 0) return false;
  const std::size_t size = m_levels[--m_level];
  while (m_g->size() < size)
  {
    const GNodeSplit* ns = add_node();
    if ( ! ns) throw std::logic_error("gts::PGraph::down: no split left");
    if (f) f(*ns);
  }
  return true;
}

//==============================================================================
// Partitions, from partition.c
//==============================================================================

unsigned partition_edges_cut(std::span<const NodeSet> parts)
{
  unsigned cuts = 0;
  for (const NodeSet& p : parts) cuts += edges_cut(p);
  return cuts / 2;
}

float partition_edges_cut_weight(std::span<const NodeSet> parts)
{
  float weight = 0.f;
  for (const NodeSet& p : parts) weight += edges_cut_weight(p);
  return weight / 2.f;
}

float partition_balance(std::span<const NodeSet> parts)
{
  float wmin = FLT_MAX, wmax = -FLT_MAX;
  for (const NodeSet& p : parts)
  {
    wmin = std::min(wmin, p.weight());
    wmax = std::max(wmax, p.weight());
  }
  return wmax - wmin;
}

void print_partition_stats(std::span<const NodeSet> parts, std::ostream& out)
{
  Range weight;
  for (const NodeSet& p : parts) weight.add(p.weight());
  weight.update();
  out << std::format("# parts: {}\n#   edge cuts: {:5} edge cuts weight: {:5g}\n#   weight: ",
                     parts.size(), partition_edges_cut(parts), partition_edges_cut_weight(parts));
  out << to_string(weight) << '\n';
}

namespace {

// partition_update: each part grows from its seed, breadth first, the
// lightest part first, until no part can grow.
void partition_update(std::vector<NodeSet>& list, const std::vector<GNodeId>& seeds, const NodeSet& g)
{
  std::vector<GraphTraverse> ts;
  ts.reserve(list.size());
  bool reinit = true;
  for (GNodeId s : seeds)
  {
    ts.emplace_back(g, s, reinit);
    reinit = false;
  }
  // GtsHeap with graph_comp_weight: the lightest part on top.
  EHeap<int> size_heap([&](const int& i) { return double(list[i].weight()); });
  for (int i = 0; i < int(list.size()); ++i) size_heap.insert(i);
  while (auto top = size_heap.remove_top())
  {
    const int i = top->first;
    const GNodeId n = ts[i].next();
    if (n.valid())
    {
      list[i].add(n);
      size_heap.insert(i);
    }
  }
}

// graph_new_seed: the neighbour of seed in g, or seed, with the smallest
// sum of distances to the nodes of g.
GNodeId graph_new_seed(const NodeSet& g, GNodeId seed)
{
  const Graph& gr = g.graph();
  unsigned sum = distance_sum(g, seed);
  GNodeId new_seed = seed;
  for (GEdgeId e : gr.edges_of(seed))
  {
    const GNodeId n = gr.neighbor(seed, e);
    if ( ! in(&g, n)) continue;
    const unsigned sum1 = distance_sum(g, n);
    if (sum1 < sum)
    {
      sum = sum1;
      new_seed = n;
    }
  }
  return new_seed;
}

} // namespace

std::vector<NodeSet> bubble_partition(const NodeSet& g, unsigned np, unsigned niter,
                                      const std::function<void(std::span<const NodeSet>)>& step)
{
  if (np == 0) throw std::invalid_argument("gts::bubble_partition: no parts");
  const Graph& gr = g.graph();
  // find_smallest_degree
  GNodeId seed;
  unsigned min = INT_MAX / 2 - 1;
  for (GNodeId n : g.nodes())
    if (const unsigned d = degree(gr, n, &g); d < min) { min = d; seed = n; }
  if ( ! seed.valid()) return {};

  // GTS's lists, built by prepending: list and seeds are newest first.
  std::vector<NodeSet> list;
  std::vector<GNodeId> part_seeds;
  std::vector<GNodeId> seeds;
  auto new_part = [&](GNodeId s) {
    NodeSet p(gr);
    p.add(s);
    list.insert(list.begin(), std::move(p));
    part_seeds.insert(part_seeds.begin(), s);
  };
  new_part(seed);
  seeds.insert(seeds.begin(), seed);
  while (--np && seed.valid())
    if ((seed = farthest(g, seeds)).valid())
    {
      new_part(seed);
      seeds.insert(seeds.begin(), seed);
    }
  partition_update(list, part_seeds, g);

  bool changed = true;
  while (changed && niter--)
  {
    changed = false;
    for (std::size_t i = 0; i < list.size(); ++i)
    {
      const GNodeId new_seed = graph_new_seed(list[i], part_seeds[i]);
      if (new_seed != part_seeds[i])
      {
        changed = true;
        part_seeds[i] = new_seed;
      }
    }
    if (changed)
    {
      for (std::size_t i = 0; i < list.size(); ++i)
      {
        NodeSet p(gr);
        p.add(part_seeds[i]);
        list[i] = std::move(p);
      }
      partition_update(list, part_seeds, g);
      if (step) step(list);
    }
  }
  return list;
}

//==============================================================================
// Bisection
//==============================================================================

namespace {

// boundary_node1 and boundary_node2: each node of a with a neighbour in b,
// mapped to the first one found.
void boundary(const NodeSet& a, const NodeSet& b, std::map<Index, GNodeId>& bg)
{
  const Graph& gr = a.graph();
  for (GNodeId n : a.nodes())
    for (GEdgeId e : gr.edges_of(n))
      if (const GNodeId n1 = gr.neighbor(n, e); in(&b, n1))
      {
        bg.emplace(n.i, n1);
        break;
      }
}

// add_unused: the nodes of g that the growing did not take go to g2.
void add_unused(const NodeSet& g, const std::vector<bool>& taken, NodeSet& g2)
{
  for (GNodeId n : g.nodes())
    if ( ! taken[n.i]) g2.add(n);
}

// The seeds by smallest degree in g (degree_cost).
EHeap<GNodeId> degree_heap(const NodeSet& g)
{
  const Graph& gr = g.graph();
  EHeap<GNodeId> h([&gr, &g](const GNodeId& n) { return double(degree(gr, n, &g)); });
  h.freeze();
  for (GNodeId n : g.nodes()) h.insert(n);
  h.thaw();
  return h;
}

} // namespace

bool check(const GraphBisection& bg)
{
  const Graph& gr = bg.g->graph();
  auto check_side = [&](const NodeSet& a, const NodeSet& b, const std::map<Index, GNodeId>& m) {
    bool ok = true;
    std::size_t nb = 0;
    for (GNodeId n : a.nodes())
    {
      const unsigned nn = degree(gr, n, &b);
      if (nn > 0) ++nb;
      if ((nn > 0) != m.contains(n.i)) ok = false;
    }
    return ok && m.size() == nb;
  };
  return check_side(bg.g1, bg.g2, bg.bg1) && check_side(bg.g2, bg.g1, bg.bg2);
}

GraphBisection ggg_bisection(const NodeSet& g, unsigned ntry)
{
  const Graph& gr = g.graph();
  const float size = g.weight() / 2.f;
  const float smin = 0.9f * size;
  float bestcost = FLT_MAX;
  std::optional<NodeSet> bestg1, bestg2;
  bool balanced = false;

  EHeap<GNodeId> seeds = degree_heap(g);
  std::optional<std::pair<GNodeId, double>> top;
  while (ntry && (top = seeds.remove_top()))
  {
    const GNodeId seed = top->first;
    NodeSet g1(gr), g2(gr);
    // node_cost: the weight of the edges to g outside g1, less that into g1.
    auto node_cost = [&](const GNodeId& n) {
      double cost = 0.;
      for (GEdgeId e : gr.edges_of(n))
      {
        const GNodeId n1 = gr.neighbor(n, e);
        if (in(&g, n1)) cost += g1.contains(n1) ? -double(gr.weight(e)) : double(gr.weight(e));
      }
      return cost;
    };
    EHeap<GNodeId> heap(node_cost);
    // The heap handle of each node, or taken for a node of g1.
    using H = EHeap<GNodeId>::Handle;
    constexpr H kOut = ~H(0), kTaken = ~H(0) - 1;
    std::vector<H> handle(gr.node_capacity(), kOut);
    std::vector<bool> taken(gr.node_capacity());
    auto add_neighbors = [&](GNodeId n) {
      for (GEdgeId e : gr.edges_of(n))
      {
        const GNodeId m = gr.neighbor(n, e);
        if ( ! in(&g, m) || handle[m.i] == kTaken) continue;
        if (handle[m.i] != kOut) heap.remove(handle[m.i]);
        handle[m.i] = heap.insert(m);
      }
    };

    g1.add(seed);
    handle[seed.i] = kTaken;
    taken[seed.i] = true;
    add_neighbors(seed);
    while (auto t = heap.remove_top())
    {
      const GNodeId n = t->first;
      if (g1.weight() + gr.weight(n) <= size)
      {
        g1.add(n);
        handle[n.i] = kTaken;
        taken[n.i] = true;
        add_neighbors(n);
      }
      else
        handle[n.i] = kOut;
    }
    add_unused(g, taken, g2);

    const float cost = edges_cut_weight(g1);
    if ( ! bestg1 || ( ! balanced && g1.weight() >= smin) || (cost < bestcost && g1.weight() >= smin))
    {
      if (bestg1) bestcost = cost;
      if (g1.weight() >= smin) balanced = true;
      bestg1 = std::move(g1);
      bestg2 = std::move(g2);
    }
    --ntry;
  }
  if ( ! bestg1) throw std::invalid_argument("gts::ggg_bisection: no seed to grow from");
  GraphBisection bg{ &g, std::move(*bestg1), std::move(*bestg2), {}, {} };
  boundary(bg.g1, bg.g2, bg.bg1);
  boundary(bg.g2, bg.g1, bg.bg2);
  return bg;
}

GraphBisection bfgg_bisection(const NodeSet& g, unsigned ntry)
{
  const Graph& gr = g.graph();
  const float size = g.weight() / 2.f;
  const float smin = 0.9f * size;
  float bestcost = FLT_MAX;
  std::optional<NodeSet> bestg1, bestg2;

  EHeap<GNodeId> seeds = degree_heap(g);
  std::optional<std::pair<GNodeId, double>> top;
  while (ntry && (top = seeds.remove_top()))
  {
    NodeSet g1(gr), g2(gr);
    std::vector<bool> taken(gr.node_capacity());
    GraphTraverse t(g, top->first, true);
    for (GNodeId n = t.next(); n.valid(); n = t.next())
      if (g1.weight() + gr.weight(n) <= size)
      {
        g1.add(n);
        taken[n.i] = true;
      }
    add_unused(g, taken, g2);

    const float cost = edges_cut_weight(g1);
    if ( ! bestg1 || (cost < bestcost && g1.weight() >= smin))
    {
      if (bestg1) bestcost = cost;
      bestg1 = std::move(g1);
      bestg2 = std::move(g2);
    }
    --ntry;
  }
  if ( ! bestg1) throw std::invalid_argument("gts::bfgg_bisection: no seed to grow from");
  GraphBisection bg{ &g, std::move(*bestg1), std::move(*bestg2), {}, {} };
  boundary(bg.g1, bg.g2, bg.bg1);
  boundary(bg.g2, bg.g1, bg.bg2);
  return bg;
}

double kl_refine(GraphBisection& bg, unsigned mmax)
{
  if (mmax == 0) throw std::invalid_argument("gts::kl_refine: mmax must be positive");
  const Graph& gr = bg.g->graph();
  using H = EHeap<GNodeId>::Handle;
  constexpr H kOut = ~H(0);
  std::vector<H> handle(gr.node_capacity(), kOut);
  EHeap<GNodeId> h1([&](const GNodeId& n) { return double(move_cost(n, bg.g1, bg.g2)); });
  EHeap<GNodeId> h2([&](const GNodeId& n) { return double(move_cost(n, bg.g2, bg.g1)); });
  h1.freeze();
  for (GNodeId n : bg.g1.nodes()) handle[n.i] = h1.insert(n);
  h1.thaw();
  h2.freeze();
  for (GNodeId n : bg.g2.nodes()) handle[n.i] = h2.insert(n);
  h2.thaw();

  std::vector<GNodeId> moves;
  moves.reserve(mmax);
  double bestcost = 0., totalcost = 0.;
  double best_balance = std::fabs(double(bg.g1.weight() - bg.g2.weight()));
  GNodeId n;
  do
  {
    NodeSet *g1, *g2;
    std::optional<std::pair<GNodeId, double>> top;
    if (bg.g1.weight() > bg.g2.weight()) { top = h1.remove_top(); g1 = &bg.g1; g2 = &bg.g2; }
    else                                 { top = h2.remove_top(); g1 = &bg.g2; g2 = &bg.g1; }
    n = top ? top->first : GNodeId();
    if (n.valid())
    {
      const double cost = top->second;
      handle[n.i] = kOut;
      g2->add(n);
      g1->remove(n);

      totalcost += cost;
      if (totalcost < bestcost)
      {
        bestcost = totalcost;
        moves.clear();
      }
      else if (totalcost == bestcost)
      {
        const double balance = std::fabs(double(g1->weight() - g2->weight()));
        if (balance < best_balance)
        {
          best_balance = balance;
          moves.clear();
        }
      }
      else
        moves.push_back(n);

      for (GEdgeId e : gr.edges_of(n))
      {
        const GNodeId n1 = gr.neighbor(n, e);
        if (n1.valid() && handle[n1.i] != kOut && bg.g->contains(n1))
        {
          EHeap<GNodeId>& h = bg.g1.contains(n1) ? h1 : h2;
          h.remove(handle[n1.i]);
          handle[n1.i] = h.insert(n1);
        }
      }
    }
  } while (n.valid() && moves.size() < mmax);

  // Undo the last moves.
  for (GNodeId m : moves)
  {
    NodeSet& g1 = bg.g1.contains(m) ? bg.g1 : bg.g2;
    NodeSet& g2 = &g1 == &bg.g1 ? bg.g2 : bg.g1;
    g2.add(m);
    g1.remove(m);
  }
  return bestcost;
}

namespace {

// update_neighbors in partition.c: the boundary maps and the heaps after n
// moved. handle marks a node moved in this refinement with kMoved.
void update_neighbors(GNodeId n, GraphBisection& bg, EHeap<GNodeId>* h1, EHeap<GNodeId>* h2,
                      std::vector<EHeap<GNodeId>::Handle>& handle)
{
  using H = EHeap<GNodeId>::Handle;
  constexpr H kOut = ~H(0), kMoved = ~H(0) - 1;
  const Graph& gr = bg.g->graph();
  for (GEdgeId e : gr.edges_of(n))
  {
    const GNodeId n1 = gr.neighbor(n, e);
    if ( ! in(bg.g, n1)) continue;
    const bool in1 = bg.g1.contains(n1);
    EHeap<GNodeId>* h = in1 ? h1 : h2;
    const NodeSet& g2 = in1 ? bg.g2 : bg.g1;
    auto& b1 = in1 ? bg.bg1 : bg.bg2;
    b1.erase(n1.i);
    if (h && handle[n1.i] != kOut && handle[n1.i] != kMoved)
    {
      h->remove(handle[n1.i]);
      handle[n1.i] = kOut;
    }
    if (degree(gr, n1, &g2))
    {
      b1.emplace(n1.i, n1);
      if (h && handle[n1.i] != kMoved) handle[n1.i] = h->insert(n1);
    }
  }
}

} // namespace

double bkl_refine(GraphBisection& bg, unsigned mmax, float imbalance)
{
  if (mmax == 0) throw std::invalid_argument("gts::bkl_refine: mmax must be positive");
  if ( ! (imbalance >= 0.f && imbalance <= 1.f)) throw std::invalid_argument("gts::bkl_refine: imbalance must be in [0, 1]");
  const Graph& gr = bg.g->graph();
  using H = EHeap<GNodeId>::Handle;
  constexpr H kOut = ~H(0), kMoved = ~H(0) - 1;
  std::vector<H> handle(gr.node_capacity(), kOut);
  EHeap<GNodeId> h1([&](const GNodeId& n) { return double(move_cost(n, bg.g1, bg.g2)); });
  EHeap<GNodeId> h2([&](const GNodeId& n) { return double(move_cost(n, bg.g2, bg.g1)); });
  h1.freeze();
  for (const auto& [i, m] : bg.bg1) handle[i] = h1.insert(GNodeId(i));
  h1.thaw();
  h2.freeze();
  for (const auto& [i, m] : bg.bg2) handle[i] = h2.insert(GNodeId(i));
  h2.thaw();

  std::vector<GNodeId> moves;
  moves.reserve(mmax);
  const double imb = double(imbalance * bg.g->weight());
  double bestcost = 0., totalcost = 0.;
  double best_balance = std::fabs(double(bg.g1.weight() - bg.g2.weight()));
  bool balanced = best_balance <= imb;
  GNodeId n;
  do
  {
    NodeSet *g1, *g2;
    std::map<Index, GNodeId> *b1, *b2;
    std::optional<std::pair<GNodeId, double>> top;
    if (bg.g1.weight() > bg.g2.weight())
    {
      top = h1.remove_top(); g1 = &bg.g1; g2 = &bg.g2; b1 = &bg.bg1; b2 = &bg.bg2;
    }
    else
    {
      top = h2.remove_top(); g1 = &bg.g2; g2 = &bg.g1; b1 = &bg.bg2; b2 = &bg.bg1;
    }
    n = top ? top->first : GNodeId();
    if (n.valid())
    {
      const double cost = top->second;
      handle[n.i] = kMoved;
      g2->add(n);
      g1->remove(n);
      b1->erase(n.i);
      if (degree(gr, n, g1)) b2->emplace(n.i, n);
      update_neighbors(n, bg, &h1, &h2, handle);

      totalcost += cost;
      const double balance = std::fabs(double(g1->weight() - g2->weight()));
      if ( ! balanced && balance <= imb)
      {
        bestcost = totalcost;
        best_balance = balance;
        balanced = true;
        moves.clear();
      }
      else if (totalcost < bestcost && (balance < best_balance || balance <= imb))
      {
        bestcost = totalcost;
        best_balance = balance;
        moves.clear();
      }
      else if (totalcost == bestcost && balance < best_balance)
      {
        best_balance = balance;
        moves.clear();
      }
      else
        moves.push_back(n);
    }
  } while (n.valid() && moves.size() < mmax);

  // Undo the last moves.
  for (GNodeId m : moves)
  {
    const bool in1 = bg.g1.contains(m);
    NodeSet& g1 = in1 ? bg.g1 : bg.g2;
    NodeSet& g2 = in1 ? bg.g2 : bg.g1;
    auto& b1 = in1 ? bg.bg1 : bg.bg2;
    auto& b2 = in1 ? bg.bg2 : bg.bg1;
    g2.add(m);
    g1.remove(m);
    b1.erase(m.i);
    if (degree(gr, m, &g1)) b2.emplace(m.i, m);
    update_neighbors(m, bg, nullptr, nullptr, handle);
  }
  return bestcost;
}

GraphBisection bisection(Graph& gr, NodeSet& g, unsigned ntry, unsigned mmax, unsigned nmin, float imbalance)
{
  PGraph pg(gr, g, nmin);
  GraphBisection bg = ggg_bisection(g, ntry);
  while (bkl_refine(bg, mmax, imbalance) != 0.) {}
  // bisection_children: the nodes of each split take its place on its side.
  auto children = [&](const GNodeSplit& ns) {
    const bool in1 = bg.g1.contains(ns.n);
    NodeSet& side = in1 ? bg.g1 : bg.g2;
    const NodeSet& other = in1 ? bg.g2 : bg.g1;
    auto& bbg = in1 ? bg.bg1 : bg.bg2;
    side.remove(ns.n);
    side.add(ns.n1);
    side.add(ns.n2);
    if (bbg.erase(ns.n.i))
    {
      if (degree(gr, ns.n1, &other) > 0) bbg.emplace(ns.n1.i, ns.n1);
      if (degree(gr, ns.n2, &other) > 0) bbg.emplace(ns.n2.i, ns.n2);
    }
  };
  while (pg.down(children))
    while (bkl_refine(bg, mmax, imbalance) != 0.) {}
  return bg;
}

namespace {

void recursive_bisection(Graph& gr, NodeSet& g, unsigned np, unsigned ntry, unsigned mmax, unsigned nmin,
                         float imbalance, std::vector<NodeSet>& list)
{
  if (np == 0)
  {
    list.insert(list.begin(), std::move(g));
    return;
  }
  GraphBisection bg = bisection(gr, g, ntry, mmax, nmin, imbalance);
  NodeSet g1 = std::move(bg.g1), g2 = std::move(bg.g2);
  recursive_bisection(gr, g1, np - 1, ntry, mmax, nmin, imbalance, list);
  recursive_bisection(gr, g2, np - 1, ntry, mmax, nmin, imbalance, list);
}

} // namespace

std::vector<NodeSet> recursive_bisection(Graph& gr, NodeSet& g, unsigned n, unsigned ntry, unsigned mmax,
                                         unsigned nmin, float imbalance)
{
  if (n == 0) throw std::invalid_argument("gts::recursive_bisection: n must be positive");
  std::vector<NodeSet> list;
  GraphBisection bg = bisection(gr, g, ntry, mmax, nmin, imbalance);
  NodeSet g1 = std::move(bg.g1), g2 = std::move(bg.g2);
  recursive_bisection(gr, g1, n - 1, ntry, mmax, nmin, imbalance, list);
  recursive_bisection(gr, g2, n - 1, ntry, mmax, nmin, imbalance, list);
  return list;
}

} // namespace gts
