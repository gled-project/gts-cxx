// Copyright (C) 1999 Stéphane Popinet (GTS: surface.c, misc.c, point.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/io.h"
#include "gts/Stats.h"
#include "io_detail.h"

#include <cstdlib>
#include <cstring>
#include <format>
#include <fstream>
#include <istream>
#include <iterator>
#include <ostream>
#include <vector>

namespace gts {

using detail::Tokenizer;
using detail::is_number;

//==============================================================================
// Reading: gts_surface_read
//==============================================================================

std::expected<Surface*, ReadError> detail::read_gts(Geometrium& g, Tokenizer& f, std::string name,
                                                    std::vector<VertexId>* vertices_out,
                                                    std::vector<FaceId>* faces_out)
{
  while (f.ch() == '\n') f.next();

  auto read_count = [&](const char* what) -> long {
    if (f.type() != Tokenizer::Int) { f.fail(std::format("expecting an integer ({})", what)); return 0; }
    long n = std::atol(f.token().c_str());
    f.next();
    return n;
  };
  const long nv = read_count("number of vertices");
  const long ne = f.failed() ? 0 : read_count("number of edges");
  const long nf = f.failed() ? 0 : read_count("number of faces");
  if (f.failed()) return std::unexpected(f.error());

  bool binary = false;
  if (f.type() == Tokenizer::String)
  {
    const char* what[] = { "GtsFaceClass", "GtsEdgeClass", "GtsVertexClass" };
    for (const char* w : what)
    {
      f.next();
      if (f.type() != Tokenizer::String) { f.fail(std::format("expecting a string ({})", w)); return std::unexpected(f.error()); }
    }
    binary = f.token() == "GtsVertexBinary";
    if ( ! binary) f.skip_line();
  }
  else
    f.skip_line();

  const std::size_t nv0 = g.vertex_capacity(), ne0 = g.edge_capacity(), nf0 = g.face_capacity();
  Surface& s = g.new_surface(std::move(name));
  if (nf <= 0) return &s;

  std::vector<VertexId> vertices;
  std::vector<EdgeId>   edges;
  vertices.reserve(nv);
  edges.reserve(ne);

  // Vertices: point_read
  while (long(vertices.size()) < nv && ! f.failed())
  {
    Vec3 p;
    if (binary)
    {
      double* c[] = { &p.x, &p.y, &p.z };
      const char* axis[] = { "x", "y", "z" };
      for (int k = 0; k < 3 && ! f.failed(); ++k)
        if ( ! f.read_raw(c[k], sizeof(double)))
          f.fail(std::format("expecting a binary number ({} coordinate)", axis[k]));
    }
    else
    {
      double* c[] = { &p.x, &p.y, &p.z };
      const char* axis[] = { "x", "y", "z" };
      for (int k = 0; k < 3 && ! f.failed(); ++k)
      {
        if ( ! is_number(f.type())) { f.fail(std::format("expecting a number ({} coordinate)", axis[k])); break; }
        *c[k] = std::atof(f.token().c_str());
        f.next();
      }
      if ( ! f.failed()) f.skip_line();
    }
    if ( ! f.failed()) vertices.push_back(g.add_vertex(p));
  }
  if (binary && ! f.failed())
  {
    f.next();
    f.skip_line();
  }

  // Edges
  auto read_index = [&](long max, const char* what, const char* kind) -> long {
    if (f.type() != Tokenizer::Int) { f.fail(std::format("expecting an integer ({})", what)); return 0; }
    const long i = std::atol(f.token().c_str());
    if (i <= 0 || i > max) { f.fail(std::format("{} index `{}' is out of range `[1,{}]'", kind, i, max)); return 0; }
    return i;
  };
  while (long(edges.size()) < ne && ! f.failed())
  {
    const long a = read_index(long(vertices.size()), "first vertex index", "vertex");
    if (f.failed()) break;
    f.next();
    const long b = read_index(long(vertices.size()), "second vertex index", "vertex");
    if (f.failed()) break;
    if (a == b) { f.fail(std::format("edge from vertex `{}' to itself", a)); break; }
    edges.push_back(g.add_edge(vertices[a - 1], vertices[b - 1]));
    f.next();
    f.skip_line();
  }

  // Faces
  long n = 0;
  while (n < nf && ! f.failed())
  {
    long idx[3];
    const char* what[] = { "first edge index", "second edge index", "third edge index" };
    for (int k = 0; k < 3 && ! f.failed(); ++k)
    {
      if (k > 0) f.next();
      idx[k] = read_index(long(edges.size()), what[k], "edge");
    }
    if (f.failed()) break;
    const EdgeId e1 = edges[idx[0] - 1], e2 = edges[idx[1] - 1], e3 = edges[idx[2] - 1];
    if ( ! g.forms_triangle(e1, e2, e3))
    {
      f.fail(std::format("edges `{}', `{}' and `{}' do not form a triangle", idx[0], idx[1], idx[2]));
      break;
    }
    const FaceId t = g.add_face(e1, e2, e3);
    s.add_face(t);
    if (faces_out) faces_out->push_back(t);
    f.next();
    f.skip_line();
    ++n;
  }

  if (f.failed())
  {
    // Nothing older refers to what was read, so it can simply be cut off.
    KeepOrphans keep(g, { .faces = true, .edges = true, .vertices = true });
    g.remove_surface(s);
    g.truncate(nv0, ne0, nf0);
    if (vertices_out) vertices_out->clear();
    if (faces_out) faces_out->clear();
    return std::unexpected(f.error());
  }
  if (vertices_out) *vertices_out = std::move(vertices);
  return &s;
}

std::expected<Surface*, ReadError> read_gts(Geometrium& g, std::istream& in, std::string name)
{
  detail::Tokenizer f(std::string(std::istreambuf_iterator<char>(in), {}));
  return detail::read_gts(g, f, std::move(name));
}

std::expected<Surface*, ReadError> read_gts(Geometrium& g, const std::filesystem::path& path)
{
  std::ifstream in(path, std::ios::binary);
  if ( ! in) return std::unexpected(ReadError{ 0, 0, "cannot open " + path.string() });
  return read_gts(g, in, path.stem().string());
}

//==============================================================================
// Writing: gts_surface_write
//==============================================================================

void write_gts(const Surface& s, std::ostream& out, WriteOptions o)
{
  const Geometrium& g = s.geometrium();

  // Number the vertices and edges of s, in index order, from 1.
  std::vector<Index> vnum(g.vertex_capacity(), 0), enum_(g.edge_capacity(), 0);
  for (FaceId f : s.faces())
    for (EdgeId e : g.face(f).e)
    {
      enum_[e.i] = 1;
      for (VertexId v : g.edge(e).v) vnum[v.i] = 1;
    }
  Index nv = 0, ne = 0;
  for (auto& x : vnum)  if (x) x = ++nv;
  for (auto& x : enum_) if (x) x = ++ne;

  out << std::format("{} {} {} GtsSurface GtsFace GtsEdge GtsVertex\n", nv, ne, s.face_number());
  for (Index i = 0; i < vnum.size(); ++i)
    if (vnum[i])
    {
      const Vec3 p = g.vertex(VertexId(i)).p;
      out << std::format("{:.{}g} {:.{}g} {:.{}g}\n", p.x, o.precision, p.y, o.precision, p.z, o.precision);
    }
  for (Index i = 0; i < enum_.size(); ++i)
    if (enum_[i])
    {
      const Edge& e = g.edge(EdgeId(i));
      out << std::format("{} {}\n", vnum[e.v[0].i], vnum[e.v[1].i]);
    }
  for (FaceId f : s.faces())
  {
    const auto& e = g.face(f).e;
    out << std::format("{} {} {}\n", enum_[e[0].i], enum_[e[1].i], enum_[e[2].i]);
  }
}

bool write_gts(const Surface& s, const std::filesystem::path& path, WriteOptions o)
{
  std::ofstream out(path);
  if ( ! out) return false;
  write_gts(s, out, o);
  return bool(out);
}

//==============================================================================
// Geomview and VTK
//==============================================================================

namespace {

// The vertices of s in index order, numbered from 0; kNone for the others.
std::vector<Index> number_vertices(const Surface& s, Index& n)
{
  const Geometrium& g = s.geometrium();
  std::vector<Index> num(g.vertex_capacity(), kNone);
  for (FaceId f : s.faces())
    for (VertexId v : g.vertices_of(f)) num[v.i] = 0;
  n = 0;
  for (auto& x : num) if (x != kNone) x = n++;
  return num;
}

} // namespace

void write_oogl(const Surface& s, std::ostream& out, const OoglColors& c)
{
  const Geometrium& g = s.geometrium();
  Index nv;
  const std::vector<Index> num = number_vertices(s, nv);
  out << (c.vertex ? "COFF " : "OFF ");
  out << std::format("{} {} {}\n", nv, s.face_number(), s.edge_number());
  for (Index i = 0; i < num.size(); ++i)
  {
    if (num[i] == kNone) continue;
    const Vec3 p = g.vertex(VertexId(i)).p;
    out << std::format("{:g} {:g} {:g}", p.x, p.y, p.z);
    if (c.vertex)
    {
      const Rgb k = c.vertex(VertexId(i));
      out << std::format(" {:g} {:g} {:g} 1.0\n", k.r, k.g, k.b);
    }
    else out << '\n';
  }
  for (FaceId f : s.faces())
  {
    const auto v = g.vertices_of(f);
    out << std::format("3 {} {} {}", num[v[0].i], num[v[1].i], num[v[2].i]);
    if (c.face)
    {
      const Rgb k = c.face(f);
      out << std::format(" {:g} {:g} {:g}\n", k.r, k.g, k.b);
    }
    else out << '\n';
  }
}

void write_oogl_boundary(const Surface& s, std::ostream& out, const OoglColors& c)
{
  const Geometrium& g = s.geometrium();
  out << "LIST {\n";
  for (EdgeId e : s.edges())
  {
    if ( ! boundary_face(g, e, &s).valid()) continue;
    const Segment x = g.segment(e);
    if (c.edge)
    {
      const Rgb k = c.edge(e);
      out << std::format("VECT 1 2 1 2 1 {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} 1.\n",
                         x.a.x, x.a.y, x.a.z, x.b.x, x.b.y, x.b.z, k.r, k.g, k.b);
    }
    else
      out << std::format("VECT 1 2 0 2 0 {:g} {:g} {:g} {:g} {:g} {:g}\n", x.a.x, x.a.y, x.a.z, x.b.x, x.b.y, x.b.z);
  }
  out << "}\n";
}

void write_vtk(const Surface& s, std::ostream& out)
{
  const Geometrium& g = s.geometrium();
  Index nv;
  const std::vector<Index> num = number_vertices(s, nv);
  out << std::format("# vtk DataFile Version 2.0\n"
                     "Generated by GTS\n"
                     "ASCII\n"
                     "DATASET POLYDATA\n"
                     "POINTS {} float\n", nv);
  for (Index i = 0; i < num.size(); ++i)
    if (num[i] != kNone)
    {
      const Vec3 p = g.vertex(VertexId(i)).p;
      out << std::format("{:g} {:g} {:g}\n", p.x, p.y, p.z);
    }
  out << std::format("POLYGONS {} {}\n", s.face_number(), s.face_number()*4);
  for (FaceId f : s.faces())
  {
    const auto v = g.vertices_of(f);
    out << std::format("3 {} {} {}\n", num[v[0].i], num[v[1].i], num[v[2].i]);
  }
}

} // namespace gts
