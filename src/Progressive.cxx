// Copyright (C) 1999 Stéphane Popinet (GTS: split.c, psurface.c, hsurface.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-11.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/Progressive.h"
#include "gts/Stats.h"
#include "gts/Topology.h"
#include "io_detail.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <format>
#include <istream>
#include <iterator>
#include <ostream>
#include <stdexcept>

namespace gts {

// GTS rewires the faces of a collapse or an expansion before it moves the
// ends of the edges, so for a while they are not triangles.
struct detail::SplitAccess
{
  static void relink_face(Geometrium& g, FaceId f, const std::array<EdgeId, 3>& e) { g.relink_face(f, e); }
};

//==============================================================================
// Vertex splits, from split.c
//==============================================================================

namespace {

// The CFACE_ flags of split.c.
enum : unsigned
{
  kOrientation = 0x1,  // f ran v1 -> v2 -> its third vertex
  kVvsDirect   = 0x2,  // where the edge to the split's vertex is in t
  kE1          = 0x4,  // the edge to v1 was kept for another collapsed face
  kE2          = 0x8,
  kKeepVvs     = 0x10  // the edge to the split's vertex is kept on expansion
};

bool uses(const Geometrium& g, EdgeId e, VertexId v)
{
  const Edge& E = g.edge(e);
  return E.v[0] == v || E.v[1] == v;
}

// gts_segments_touch
bool touch(const Geometrium& g, EdgeId a, EdgeId b)
{
  const Edge& A = g.edge(a);
  return uses(g, b, A.v[0]) || uses(g, b, A.v[1]);
}

// TRIANGLE_REPLACE_EDGE: t takes with in place of e; it becomes the newest
// face of with.
void replace_edge_of(Geometrium& g, FaceId t, EdgeId e, EdgeId with)
{
  auto E = g.face(t).e;
  for (EdgeId& x : E)
    if (x == e) { x = with; break; }
  detail::SplitAccess::relink_face(g, t, E);
}

// replace_edge_collapse: moves the faces on e, but those of the collapse in
// progress, to with, and lists them in a. Faces of the collapse in
// progress that stay on e get the flag edge_flag, and cf keeps its edge to
// the split's vertex. Removes e when it is left without faces. Returns the
// last face moved.
FaceId replace_edge_collapse(Geometrium& g, EdgeId e, EdgeId with, unsigned& cf_flags, unsigned edge_flag,
                             std::unordered_map<Index, unsigned>& pending, std::vector<FaceId>& a,
                             const std::function<void(EdgeId)>& on_removed_edge)
{
  std::vector<FaceId> ring;
  for (FaceId t : g.faces_of(e)) ring.push_back(t);
  FaceId rt;
  for (FaceId t : ring)
  {
    if (auto it = pending.find(t.i); it != pending.end())
    {
      it->second = edge_flag;
      cf_flags |= kKeepVvs;
    }
    else
    {
      replace_edge_of(g, t, e, with);
      a.push_back(t);
      rt = t;
    }
  }
  if ( ! g.edge(e).face.valid())
  {
    if (on_removed_edge) on_removed_edge(e);
    g.kill_edge(e);
  }
  return rt;
}

// cface_new
void cface_new(Geometrium& g, FaceId f, EdgeId e, VertexId v1, VertexId v2, VertexSplit& vs, VertexSplit::CFace& cf,
               std::unordered_map<Index, unsigned>& pending, const std::function<void(EdgeId)>& on_removed_edge)
{
  cf = {};
  cf.f = f;
  cf.flags = pending[f.i];
  pending.erase(f.i);

  auto [e1, e2, e3] = g.face(f).e;
  // ROTATE_ORIENT: e1 and e2 the two other edges, in f's order after e.
  if      (e1 == e) { e1 = e2; e2 = e3; }
  else if (e2 == e) { e2 = e1; e1 = e3; }
  if (uses(g, e1, v2))
  {
    cf.flags |= kOrientation; // v1 -> v2 -> w
    std::swap(e1, e2);
  }
  const VertexId w = g.other_vertex(e1, v1);
  g.park(f);

  EdgeId vvs;
  if (cf.flags & (kE1 | kE2))
  {
    vvs = g.find_edge(vs.v, w);
    if ( ! vvs.valid()) throw std::logic_error("gts::collapse: the edge kept for a collapsed face is gone");
  }
  else
    vvs = g.add_edge(w, vs.v);

  const FaceId t1 = replace_edge_collapse(g, e1, vvs, cf.flags, kE1, pending, cf.a1, on_removed_edge);
  const FaceId t2 = replace_edge_collapse(g, e2, vvs, cf.flags, kE2, pending, cf.a2, on_removed_edge);
  cf.t = t1.valid() ? t1 : t2;
  if ( ! cf.t.valid()) throw std::logic_error("gts::collapse: a face on the edge has no neighbour");

  // The flag to find vvs in t again.
  const auto& T = g.face(cf.t).e;
  const EdgeId x = T[0] == vvs ? T[1] : T[1] == vvs ? T[2] : T[0];
  if (uses(g, x, w)) cf.flags |= kVvsDirect;
}

// find_vvs: the edge of t from vs to its other end w.
void find_vvs(const Geometrium& g, VertexId vs, FaceId t, VertexId& w, EdgeId& vvs, bool orientation)
{
  auto [e1, e2, e3] = g.face(t).e;
  if (uses(g, e2, vs))      { const EdgeId tmp = e1; e1 = e2; e2 = e3; e3 = tmp; }
  else if (uses(g, e3, vs)) { const EdgeId tmp = e1; e1 = e3; e3 = e2; e2 = tmp; }
  if (uses(g, e2, vs) || ! touch(g, e1, e2)) { const EdgeId tmp = e1; e1 = e2; e2 = e3; e3 = tmp; }
  vvs = orientation ? e1 : e3;
  w = g.other_vertex(vvs, vs);
}

// The marks of expand.fig: the vertex a face goes to, toggled when it is
// reached a second time.
using Marks = std::unordered_map<Index, VertexId>;

void toggle(Marks& marks, FaceId t, VertexId v)
{
  auto it = marks.find(t.i);
  if (it == marks.end()) marks.emplace(t.i, v);
  else
  {
    if (it->second != v) throw std::logic_error("gts::expand: a face is marked for both vertices");
    marks.erase(it);
  }
}

// replace_edge_expand
void replace_edge_expand(Geometrium& g, EdgeId e, EdgeId with, const std::vector<FaceId>& a, VertexId v, Marks& marks)
{
  for (FaceId t : a)
  {
    replace_edge_of(g, t, e, with);
    toggle(marks, t, v);
  }
}

// cface_expand
void cface_expand(Geometrium& g, VertexSplit::CFace& cf, EdgeId e, VertexId v1, VertexId v2, VertexId vs, Surface& s,
                  Marks& marks)
{
  VertexId w;
  EdgeId vvs;
  find_vvs(g, vs, cf.t, w, vvs, cf.flags & kVvsDirect);

  const EdgeId e1 = cf.flags & kE1 ? g.find_edge(v1, w) : g.add_edge(w, v1);
  const EdgeId e2 = cf.flags & kE2 ? g.find_edge(v2, w) : g.add_edge(w, v2);
  if ( ! e1.valid() || ! e2.valid()) throw std::logic_error("gts::expand: an edge kept for a collapsed face is gone");

  replace_edge_expand(g, vvs, e1, cf.a1, v1, marks);
  replace_edge_expand(g, vvs, e2, cf.a2, v2, marks);
  if ( ! (cf.flags & kKeepVvs))
  {
    if (g.edge(vvs).face.valid()) throw std::logic_error("gts::expand: the edge to the split's vertex still has faces");
    g.kill_edge(vvs);
  }

  const bool orientation = cf.flags & kOrientation;
  const EdgeId a = e, b = orientation ? e2 : e1, c = orientation ? e1 : e2;
  if (cf.f.valid()) g.unpark(cf.f, a, b, c);
  else              cf.f = g.add_face(a, b, c);
  s.add_face(cf.f);
}

} // namespace

bool collapsible(const Geometrium& g, const VertexSplit& vs)
{
  if (g.dead(vs.v1) || g.dead(vs.v2)) return false;
  const EdgeId e = g.find_edge(vs.v1, vs.v2);
  if ( ! e.valid()) return false;
  std::vector<FaceId> ring;
  for (FaceId f : g.faces_of(e)) ring.push_back(f);
  if (ring.empty()) return false;
  // The tests of PSurface's collapses.
  if ( ! g.collapsible(e) || (ring.size() > 2 && edge_contacts(g, e) > 1)) return false;
  auto on_e = [&](FaceId t) { return std::find(ring.begin(), ring.end(), t) != ring.end(); };
  for (FaceId f : ring)
  {
    bool neighbour = false;
    for (EdgeId x : g.face(f).e)
      if (x != e)
        for (FaceId t : g.faces_of(x))
          if ( ! on_e(t)) neighbour = true;
    if ( ! neighbour) return false;
  }
  return true;
}

void collapse(Geometrium& g, VertexSplit& vs, const std::function<void(EdgeId)>& on_removed_edge)
{
  const VertexId v = vs.v, v1 = vs.v1, v2 = vs.v2;
  if ( ! g.parked(v) && g.vertex(v).edge.valid()) throw std::logic_error("gts::collapse: the split's vertex has edges");
  if ( ! collapsible(g, vs)) throw std::logic_error("gts::collapse: the split cannot be collapsed now");
  if (g.parked(v)) g.unpark(v);
  const EdgeId e = g.find_edge(v1, v2);

  std::vector<FaceId> ring;
  for (FaceId f : g.faces_of(e)) ring.push_back(f);

  // The faces of this collapse, with the flags that cface_new() sets on
  // those still waiting (GTS marked them as CFaces at once).
  std::unordered_map<Index, unsigned> pending;
  for (FaceId f : ring) pending.emplace(f.i, 0u);
  vs.cfaces.assign(ring.size(), {});
  for (std::size_t k = 0; k < ring.size(); ++k)
    cface_new(g, ring[k], e, v1, v2, vs, vs.cfaces[k], pending, on_removed_edge);
  g.kill_edge(e);

  // The edges of v1, then those of v2, go to v ahead of its own, each list
  // in its order, as GTS concatenated them.
  for (VertexId x : { v1, v2 })
  {
    std::vector<EdgeId> edges;
    for (EdgeId y : g.edges_of(x)) edges.push_back(y);
    for (auto it = edges.rbegin(); it != edges.rend(); ++it) g.reattach(*it, x, v);
  }
  g.park(v1);
  g.park(v2);
}

void expand(Geometrium& g, VertexSplit& vs, Surface& s)
{
  const VertexId v = vs.v, v1 = vs.v1, v2 = vs.v2;
  if ( ! g.parked(v1) || ! g.parked(v2) || g.dead(v))
    throw std::logic_error("gts::expand: the split is not collapsed");
  for (VertexId x : { v1, v2 }) g.unpark(x);
  const EdgeId e = g.add_edge(v1, v2);
  Marks marks;
  for (auto& cf : vs.cfaces) cface_expand(g, cf, e, v1, v2, v, s, marks);

  // expand.fig: each edge of v whose faces carry a mark goes to the vertex
  // marked, and passes the marks on, round and round until v has no edges.
  while (g.vertex(v).edge.valid())
  {
    std::vector<EdgeId> edges;
    for (EdgeId x : g.edges_of(v)) edges.push_back(x);
    bool changed = false;
    for (EdgeId x : edges)
    {
      VertexId with;
      for (FaceId t : g.faces_of(x))
        if (auto it = marks.find(t.i); it != marks.end()) { with = it->second; break; }
      if ( ! with.valid()) continue;
      for (FaceId t : g.faces_of(x)) toggle(marks, t, with);
      g.reattach(x, v, with);
      changed = true;
    }
    if ( ! changed) throw std::logic_error("gts::expand: edges of the split's vertex that no face leads to");
  }
  g.park(v);
}

unsigned split_height(const std::vector<VertexSplit>& splits, int root)
{
  unsigned height = 0;
  for (int c : { splits[root].child1, splits[root].child2 })
    if (c >= 0) height = std::max(height, split_height(splits, c));
  return height + 1;
}

namespace {

// Applies a renumbering to the elements the splits refer to.
void renumber(std::vector<VertexSplit>& splits, const Renumbering& r)
{
  auto face = [&](FaceId& f) { if (f.valid()) f = r(f); };
  for (VertexSplit& vs : splits)
  {
    vs.v = r(vs.v); vs.v1 = r(vs.v1); vs.v2 = r(vs.v2);
    for (auto& cf : vs.cfaces)
    {
      face(cf.f);
      face(cf.t);
      for (FaceId& x : cf.a1) face(x);
      for (FaceId& x : cf.a2) face(x);
    }
  }
}

} // namespace

//==============================================================================
// Progressive surfaces, from psurface.c
//==============================================================================

void PSurface::follow_renumbering()
{
  m_renumber = m_s->geometrium().on_renumber([this](const Renumbering& r) { renumber(m_splits, r); });
}

PSurface::PSurface(Surface& s, const CoarsenStop& stop, const EdgeCost& cost, const CollapsePosition& position,
                   double min_angle)
  : m_s(&s)
{
  Geometrium& g = s.geometrium();
  const CollapsePosition& pos = position ? position : CollapsePosition(midpoint);
  using Heap = EHeap<EdgeId>;
  auto order = [&](EdgeId e) { return edge_order_key(g, e); };
  Heap heap;
  if (cost)
    heap = Heap([&](EdgeId e) { return cost(g, e); }, order);
  else
    heap = Heap([&](EdgeId e) { const Segment x = g.segment(e); return distance2(x.a, x.b); }, order);

  // The heap handle of each edge; collapses make edges, so it grows.
  constexpr Heap::Handle kOut = Heap::Handle(-1);
  std::vector<Heap::Handle> handle;
  auto handle_of = [&](EdgeId e) -> Heap::Handle& {
    if (handle.size() <= e.i) handle.resize(e.i + 1, kOut);
    return handle[e.i];
  };
  auto insert = [&](EdgeId e) { handle_of(e) = heap.insert(e); };
  auto remove = [&](EdgeId e) {
    Heap::Handle& h = handle_of(e);
    if (h != kOut) { heap.remove(h); h = kOut; }
  };

  double max_cos2 = std::cos(min_angle);
  max_cos2 *= max_cos2;

  heap.freeze();
  for (EdgeId e : s.edges()) insert(e);
  heap.thaw();

  // The split whose v each vertex is, by vertex index.
  std::unordered_map<Index, int> split_of;
  auto split_at = [&](VertexId x) { auto it = split_of.find(x.i); return it == split_of.end() ? -1 : it->second; };

  KeepOrphans keep(g, { .edges = true });
  std::vector<EdgeId> second;
  while (auto top = heap.remove_top())
  {
    const auto [e, c] = *top;
    handle_of(e) = kOut;
    if ( ! (c < DBL_MAX) || stop(c, heap.size() - face_number(s, e)))
      break;

    // edge_collapse: a refused edge goes to the end. So does a non-manifold
    // edge that is a contact edge.
    std::size_t n_faces = 0;
    for (FaceId f : g.faces_of(e)) { (void) f; ++n_faces; }
    if ( ! g.collapsible(e) || (n_faces > 2 && edge_contacts(g, e) > 1))
    {
      handle_of(e) = heap.insert(e, DBL_MAX);
      continue;
    }
    const Vec3 p = pos(g, e);
    if (g.collapse_folds(e, p, max_cos2))
    {
      handle_of(e) = heap.insert(e, DBL_MAX);
      continue;
    }
    const VertexId v1 = g.edge(e).v[0], v2 = g.edge(e).v[1];
    VertexSplit vs;
    vs.v = g.add_merged_vertex(p, v1, v2);
    vs.v1 = v1;
    vs.v2 = v2;
    vs.child1 = split_at(v1);
    vs.child2 = split_at(v2);
    collapse(g, vs, remove);
    const int index = int(m_splits.size());
    for (int ch : { vs.child1, vs.child2 }) if (ch >= 0) m_splits[ch].parent = index;
    split_of[vs.v.i] = index;
    const VertexId v = vs.v;
    m_splits.push_back(std::move(vs));

    // update_2nd_closest_neighbors, as in coarsen().
    second.clear();
    for (EdgeId x : g.edges_of(v))
    {
      const VertexId w = g.other_vertex(x, v);
      for (EdgeId y : g.edges_of(w))
        if (std::find(second.begin(), second.end(), y) == second.end()) second.push_back(y);
    }
    std::reverse(second.begin(), second.end());
    for (EdgeId y : second)
    {
      remove(y);
      insert(y);
    }
  }

  m_pos = m_splits.size();
  m_min = s.vertex_number();
  follow_renumbering();
}

const VertexSplit* PSurface::add_vertex()
{
  if (m_pos == 0) return nullptr;
  VertexSplit& vs = m_splits[--m_pos];
  expand(m_s->geometrium(), vs, *m_s);
  return &vs;
}

const VertexSplit* PSurface::remove_vertex()
{
  if (m_pos == m_splits.size()) return nullptr;
  VertexSplit& vs = m_splits[m_pos++];
  collapse(m_s->geometrium(), vs);
  return &vs;
}

void PSurface::set_vertex_number(std::size_t n)
{
  n = std::clamp(n, min_vertex_number(), max_vertex_number());
  const std::size_t pos = m_min + m_splits.size() - n;
  while (m_pos > pos && add_vertex()) {}
  while (m_pos < pos && remove_vertex()) {}
}

void PSurface::write(std::ostream& out)
{
  Geometrium& g = m_s->geometrium();
  while (remove_vertex()) {}
  write_gts(*m_s, out);

  // The numbers of write_gts(): vertices in index order, faces in the
  // surface's order, from 1. The vertices and faces of the splits follow.
  std::unordered_map<Index, unsigned> vnum, fnum;
  for (FaceId f : m_s->faces())
    for (VertexId v : g.vertices_of(f)) vnum[v.i] = 0;
  std::vector<Index> vs_sorted;
  for (const auto& [i, n] : vnum) vs_sorted.push_back(i);
  std::sort(vs_sorted.begin(), vs_sorted.end());
  unsigned nv = 1, nf = 1;
  for (Index i : vs_sorted) vnum[i] = nv++;
  for (FaceId f : m_s->faces()) fnum[f.i] = nf++;

  out << std::format("{}\n", m_splits.size());
  while (m_pos > 0)
  {
    VertexSplit& vs = m_splits[--m_pos];
    out << std::format("{} {}\n", vnum.at(vs.v.i), vs.cfaces.size());
    for (VertexId x : { vs.v1, vs.v2 })
    {
      vnum[x.i] = nv++;
      const Vec3 p = g.vertex(x).p;
      out << std::format("{:.10g} {:.10g} {:.10g}\n", p.x, p.y, p.z);
    }
    for (const auto& cf : vs.cfaces)
    {
      out << std::format("{} {}\n", fnum.at(cf.t.i), cf.flags);
      for (const auto* a : { &cf.a1, &cf.a2 })
      {
        for (FaceId t : *a) out << std::format("{} ", fnum.at(t.i));
        out << '\n';
      }
      fnum[cf.f.i] = nf++;
    }
    expand(g, vs, *m_s);
  }
}

std::expected<std::unique_ptr<PSurface>, ReadError> PSurface::read(Geometrium& g, std::istream& in)
{
  detail::Tokenizer f(std::string(std::istreambuf_iterator<char>(in), {}));
  std::vector<VertexId> vertices;
  std::vector<FaceId> faces;
  auto s = detail::read_gts(g, f, {}, &vertices, &faces);
  if ( ! s) return std::unexpected(s.error());

  std::unique_ptr<PSurface> ps(new PSurface());
  ps->m_s = *s;
  ps->m_min = (*s)->vertex_number();
  std::size_t ns = 0;
  if (f.type() == detail::Tokenizer::Int)
  {
    ns = std::size_t(std::max(0L, std::atol(f.token().c_str())));
    if (ns > 0) f.skip_line();
  }

  // The split each vertex is v1 or v2 of, by index into vertices.
  std::unordered_map<Index, int> owner;
  auto fail = [&](std::string msg) { f.fail(std::move(msg)); };
  auto read_index = [&](std::size_t max, const char* what, const char* kind) -> std::size_t {
    if (f.type() != detail::Tokenizer::Int) { fail(std::format("expecting an integer ({})", what)); return 0; }
    const long i = std::atol(f.token().c_str());
    if (i <= 0 || std::size_t(i) > max) { fail(std::format("{} index `{}' is out of range `[1,{}]'", kind, i, max)); return 0; }
    return std::size_t(i);
  };
  auto read_point = [&](Vec3& p) {
    double* c[] = { &p.x, &p.y, &p.z };
    const char* axis[] = { "x", "y", "z" };
    for (int k = 0; k < 3 && ! f.failed(); ++k)
    {
      if ( ! detail::is_number(f.type())) { fail(std::format("expecting a number ({} coordinate)", axis[k])); return; }
      *c[k] = std::atof(f.token().c_str());
      f.next();
    }
  };
  // The rest of the line and its newline, but not the blank lines after
  // it, which may be lists without faces. GTS skipped them, and failed to
  // read what it wrote for a split at a boundary; see QUIRKS.md.
  auto end_line = [&] {
    while ( ! f.failed() && f.type() != detail::Tokenizer::None && f.ch() != '\n') f.next();
    if (f.ch() == '\n') f.next();
  };
  // A face index must name a face that exists, not one of the split being
  // read, which gets its face when the split is expanded.
  auto read_face = [&]() -> FaceId {
    const std::size_t it = read_index(faces.size(), "face index", "face");
    if (f.failed()) return {};
    if ( ! faces[it - 1].valid()) fail(std::format("face index `{}' names a face of the split being read", it));
    return faces[it - 1];
  };
  // An empty list is an empty line; the end of the file is not one.
  auto read_faces = [&](std::vector<FaceId>& a) {
    if ( ! f.failed() && f.type() == detail::Tokenizer::None) fail("unexpected end of file (list of faces)");
    while ( ! f.failed() && f.type() == detail::Tokenizer::Int)
    {
      const FaceId t = read_face();
      if (f.failed()) return;
      a.push_back(t);
      f.next();
    }
    if ( ! f.failed() && f.ch() != '\n' && f.type() != detail::Tokenizer::None)
      fail("expecting an integer (face index)");
    end_line();
  };

  // gts_psurface_read_vertex, until the end.
  while (ps->m_pos < ns && ! f.failed() && f.type() != detail::Tokenizer::None)
  {
    VertexSplit vs;
    const std::size_t iv = read_index(vertices.size(), "vertex index", "vertex");
    if (f.failed()) break;
    vs.v = vertices[iv - 1];
    f.next();
    if (f.type() != detail::Tokenizer::Int) { fail("expecting an integer (ncf)"); break; }
    const long ncf = std::atol(f.token().c_str());
    if (ncf < 1) { fail(std::format("a split with `{}' faces", ncf)); break; }
    f.next();
    f.skip_line();

    const int index = int(ps->m_splits.size());
    for (VertexId* x : { &vs.v1, &vs.v2 })
    {
      Vec3 p;
      read_point(p);
      if (f.failed()) break;
      *x = g.add_vertex(p);
      g.park(*x);
      owner[x->i] = index;
      vertices.push_back(*x);
      f.skip_line();
    }
    // The faces of the split get their numbers as they come back.
    std::vector<std::size_t> slots;
    for (long k = 0; k < ncf && ! f.failed(); ++k)
    {
      VertexSplit::CFace cf;
      cf.t = read_face();
      if (f.failed()) break;
      f.next();
      if (f.type() != detail::Tokenizer::Int) { fail("expecting an integer (flags)"); break; }
      cf.flags = unsigned(std::atol(f.token().c_str()));
      f.next();
      end_line();
      read_faces(cf.a1);
      if ( ! f.failed()) read_faces(cf.a2);
      if (f.failed()) break;
      slots.push_back(faces.size());
      faces.push_back(FaceId());
      vs.cfaces.push_back(std::move(cf));
    }
    if (f.failed()) break;

    if (auto it = owner.find(vs.v.i); it != owner.end())
    {
      VertexSplit& parent = ps->m_splits[it->second];
      if (parent.v1 == vs.v) parent.child1 = index; else parent.child2 = index;
      vs.parent = it->second;
    }
    ps->m_splits.push_back(std::move(vs));
    ++ps->m_pos;
    VertexSplit& added = ps->m_splits.back();
    expand(g, added, **s);
    for (std::size_t k = 0; k < slots.size(); ++k) faces[slots[k]] = added.cfaces[k].f;
  }
  if ( ! f.failed() && ps->m_pos < ns) fail(std::format("only `{}' splits read out of `{}'", ps->m_pos, ns));
  if (f.failed()) return std::unexpected(f.error());

  // gts_psurface_close: the splits read, last first, so that the finest
  // comes first, as after gts_psurface_new.
  auto& sp = ps->m_splits;
  const int n = int(sp.size());
  std::reverse(sp.begin(), sp.end());
  auto flip = [n](int& i) { if (i >= 0) i = n - 1 - i; };
  for (auto& vs : sp) { flip(vs.child1); flip(vs.child2); flip(vs.parent); }
  ps->m_pos = 0;
  ps->follow_renumbering();
  return ps;
}

//==============================================================================
// Hierarchical surfaces, from hsurface.c
//==============================================================================

HSurface::HSurface(PSurface& ps, SplitKey expand_key, SplitKey collapse_key)
  : m_s(ps.m_s), m_expand_key(std::move(expand_key)), m_collapse_key(std::move(collapse_key)),
    m_expandable([this](const int& i) { return m_expand_key(m_splits[i]); }),
    m_collapsable([this](const int& i) { return m_collapse_key(m_splits[i]); })
{
  while (ps.remove_vertex()) {}
  ps.m_renumber.reset();
  m_splits = std::move(ps.m_splits);
  ps.m_splits.clear();
  ps.m_pos = 0;

  const int n = int(m_splits.size());
  m_nchild.assign(n, 0);
  m_expanded.assign(n, false);
  m_index.assign(n, kNone);
  Geometrium& g = m_s->geometrium();
  for (int i = 0; i < n; ++i) record_owner(i, true);
  for (int i = n - 1; i >= 0; --i)
  {
    VertexSplit& vs = m_splits[i];
    if (vs.parent < 0) m_roots.insert(m_roots.begin(), i);
    m_nchild[i] = (vs.child1 < 0) + (vs.child2 < 0);
    gts::expand(g, vs, *m_s);
    record_owner(i, false);
    m_expanded[i] = true;
    if (m_nchild[i] == 2) insert(m_collapsable, i);
  }
  m_nvertex = m_s->vertex_number();
  m_renumber = g.on_renumber([this](const Renumbering& r) {
    gts::renumber(m_splits, r);
    std::unordered_map<Index, int> owner;
    for (const auto& [f, i] : m_owner) owner.emplace(r.faces[f], i);
    m_owner.swap(owner);
  });
}

void HSurface::record_owner(int i, bool collapsed)
{
  for (const auto& cf : m_splits[i].cfaces)
    if (cf.f.valid())
    {
      if (collapsed) m_owner[cf.f.i] = i;
      else           m_owner.erase(cf.f.i);
    }
}

void HSurface::update_keys()
{
  m_expandable.update();
  m_collapsable.update();
}

bool HSurface::can_collapse(int i) const
{
  return m_expanded[i] && m_nchild[i] == 2 && collapsible(m_s->geometrium(), m_splits[i]);
}

bool HSurface::collapse(int i)
{
  if ( ! m_expanded[i] || m_nchild[i] != 2)
    throw std::logic_error("gts::HSurface::collapse: the split is not collapsable");
  if ( ! collapsible(m_s->geometrium(), m_splits[i])) return false;
  gts::collapse(m_s->geometrium(), m_splits[i]);
  record_owner(i, true);
  --m_nvertex;
  m_nchild[i] = 0;
  m_expanded[i] = false;
  remove(m_collapsable, i);
  insert(m_expandable, i);
  const VertexSplit& vs = m_splits[i];
  for (int c : { vs.child1, vs.child2 }) if (c >= 0) remove(m_expandable, c);
  if (vs.parent >= 0 && ++m_nchild[vs.parent] == 2) insert(m_collapsable, vs.parent);
  return true;
}

void HSurface::expand(int i)
{
  if (m_expanded[i] || m_index[i] == kNone)
    throw std::logic_error("gts::HSurface::expand: the split is not expandable");
  // The splits holding faces that this one gives back.
  for (const auto& cf : m_splits[i].cfaces)
    for (const auto* a : { &cf.a1, &cf.a2 })
      for (FaceId t : *a)
        if (auto it = m_owner.find(t.i); it != m_owner.end()) force_expand(it->second);
  if (m_expanded[i]) return; // expanded on the way

  gts::expand(m_s->geometrium(), m_splits[i], *m_s);
  record_owner(i, false);
  ++m_nvertex;
  m_nchild[i] = 2;
  m_expanded[i] = true;
  remove(m_expandable, i);
  insert(m_collapsable, i);
  const VertexSplit& vs = m_splits[i];
  for (int c : { vs.child1, vs.child2 }) if (c >= 0) insert(m_expandable, c);
  if (vs.parent >= 0 && m_nchild[vs.parent]-- == 2) remove(m_collapsable, vs.parent);
}

void HSurface::force_expand(int i)
{
  if (m_expanded[i]) return;
  const int parent = m_splits[i].parent;
  if (parent >= 0 && ! m_expanded[parent]) force_expand(parent);
  expand(i);
}

void HSurface::traverse(bool pre_order, int depth, const std::function<bool(int)>& f) const
{
  // gts_split_traverse, with or without a depth.
  std::function<bool(int, int)> walk = [&](int i, int d) {
    if (pre_order && f(i)) return true;
    if (depth < 0 || d > 1)
      for (int c : { m_splits[i].child1, m_splits[i].child2 })
        if (c >= 0 && walk(c, d - 1)) return true;
    return ! pre_order && f(i);
  };
  for (int r : m_roots)
    if (walk(r, depth)) return;
}

void HSurface::foreach(bool pre_order, const std::function<bool(int)>& f)
{
  while (auto top = m_expandable.top()) expand(top->first);
  bool stop = false;
  for (int i = 0; i < int(m_splits.size()) && ! stop; ++i)
  {
    if (pre_order)
    {
      stop = f(i);
      if ( ! stop && ! collapse(i)) throw std::logic_error("gts::HSurface::foreach: a split cannot be collapsed");
    }
    else
    {
      if ( ! collapse(i)) throw std::logic_error("gts::HSurface::foreach: a split cannot be collapsed");
      stop = f(i);
    }
  }
}

unsigned HSurface::height() const
{
  unsigned h = 0;
  for (int r : m_roots) h = std::max(h, split_height(m_splits, r));
  return h;
}

} // namespace gts
