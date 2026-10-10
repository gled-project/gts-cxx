// Copyright (C) 1999 Stéphane Popinet (GTS: surface.c, misc.c, point.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "gts/io.h"

#include <cstdlib>
#include <cstring>
#include <format>
#include <fstream>
#include <istream>
#include <iterator>
#include <ostream>
#include <vector>

namespace gts {

namespace {

//==============================================================================
// Tokenizer: GtsFile from misc.c, for what surfaces need
//==============================================================================

class Tokenizer
{
public:
  enum Type { None, Int, Float, String, Error, Char }; // Char: a one-character token in kTokens

  explicit Tokenizer(std::string buf) : m_buf(std::move(buf)) { next(); }

  Type               type() const  { return m_type; }
  char               ch() const    { return m_type == Char ? m_token[0] : 0; }
  const std::string& token() const { return m_token; }
  bool               failed() const { return m_type == Error; }
  ReadError          error() const { return m_error; }

  void fail(std::string msg)
  {
    if (m_type == Error) return;
    m_error = { m_line, m_pos, std::move(msg) };
    m_type = Error;
  }

  // gts_file_next_token
  void next()
  {
    if (m_type == Error) return;
    m_token.clear();
    int c;
    if (m_next != 0)
    {
      if (is_token_char(m_next))
      {
        m_line = m_curline; m_pos = m_curpos - 1;
        m_token.push_back(char(m_next));
        m_type = Char;
        m_next = 0;
        return;
      }
      c = m_next;
      m_next = 0;
    }
    else
      c = getc();

    m_type = None;
    bool in_string = false;
    while (c != EOF && ( ! in_string || ! is_delimiter(c)))
    {
      if (in_string)
      {
        if (is_token_char(c)) { m_next = c; break; }
        m_token.push_back(char(c));
      }
      else if ( ! is_delimiter(c))
      {
        in_string = true;
        m_line = m_curline; m_pos = m_curpos - 1;
        m_token.push_back(char(c));
        if (is_token_char(c)) { m_type = Char; break; }
      }
      c = getc();
    }
    if (m_type == Error) return;
    if (m_type == None && ! m_token.empty()) m_type = classify(m_token);
  }

  // gts_file_first_token_after (f, '\n')
  void skip_line()
  {
    while (m_type != Error && m_type != None && ch() != '\n') next();
    while (ch() == '\n') next();
  }

  // gts_file_read: raw bytes, counting lines and positions.
  bool read_raw(void* dst, std::size_t n)
  {
    if (m_type == Error || m_at + n > m_buf.size()) return false;
    std::memcpy(dst, m_buf.data() + m_at, n);
    for (std::size_t i = 0; i < n; ++i)
    {
      ++m_curpos;
      if (m_buf[m_at + i] == '\n') { ++m_curline; m_curpos = 1; }
    }
    m_at += n;
    return true;
  }

private:
  static constexpr const char* kDelimiters = " \t";
  static constexpr const char* kComments   = "#!";
  static constexpr const char* kTokens     = "\n{}()=";

  static bool in(int c, const char* s) { return c != EOF && c != 0 && std::strchr(s, c) != nullptr; }
  static bool is_delimiter(int c)  { return in(c, kDelimiters); }
  static bool is_token_char(int c) { return in(c, kTokens); }

  int raw()
  {
    if (m_at >= m_buf.size()) return EOF;
    return (unsigned char) m_buf[m_at++];
  }

  // gts_file_getc: skips comments, including their newline, and counts
  // lines, positions and braces.
  int getc()
  {
    if (m_type == Error) return EOF;
    int c = raw();
    ++m_curpos;
    while (in(c, kComments))
    {
      while (c != EOF && c != '\n') c = raw();
      if (c == '\n') { ++m_curline; m_curpos = 1; c = raw(); }
    }
    switch (c)
    {
      case '\n': ++m_curline; m_curpos = 1; break;
      case '{':  ++m_scope; break;
      case '}':
        if (m_scope == 0)
        {
          m_line = m_curline; m_pos = m_curpos - 1;
          m_error = { m_line, m_pos, "no matching opening brace" };
          m_type = Error;
          c = EOF;
        }
        else
          --m_scope;
        break;
    }
    return c;
  }

  static Type classify(const std::string& t)
  {
    auto all = [&](const char* set) { return t.find_first_not_of(set) == std::string::npos; };
    if (all("+-"))                return String;
    if (all("+-0123456789"))      return Int;
    if (all("+-eE."))             return String;
    if (all("+-0123456789eE."))   return Float;
    if (t.starts_with("0x") || t.starts_with("-0x") || t.starts_with("+0x"))
    {
      if (all("+-0123456789abcdefx"))   return Int;
      if (all("+-0123456789abcdefx.p")) return Float;
    }
    return String;
  }

  std::string m_buf;
  std::size_t m_at = 0;
  std::string m_token;
  Type        m_type = None;
  int         m_next = 0;
  unsigned    m_line = 1, m_pos = 1, m_curline = 1, m_curpos = 1, m_scope = 0;
  ReadError   m_error;
};

bool is_number(Tokenizer::Type t) { return t == Tokenizer::Int || t == Tokenizer::Float; }

} // namespace

//==============================================================================
// Reading: gts_surface_read
//==============================================================================

std::expected<Surface*, ReadError> read_gts(Geometrium& g, std::istream& in, std::string name)
{
  Tokenizer f(std::string(std::istreambuf_iterator<char>(in), {}));

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
    s.add_face(g.add_face(e1, e2, e3));
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
    return std::unexpected(f.error());
  }
  return &s;
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

} // namespace gts
