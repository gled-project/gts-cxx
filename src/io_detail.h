// The tokenizer of .gts files and the reader of a surface from it, shared by
// io.cxx and the progressive surfaces, whose files start with a surface.
//
// Copyright (C) 1999 Stéphane Popinet (GTS: misc.c, surface.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "gts/io.h"

#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace gts::detail {

// GtsFile from misc.c, for what surfaces need.
class Tokenizer
{
public:
  enum Type { None, Int, Float, String, Error, Char }; // Char: a one-character token in kTokens

  explicit Tokenizer(std::string buf) : m_buf(std::move(buf)) { next(); }

  Type               type() const  { return m_type; }
  char               ch() const    { return m_type == Char ? m_token[0] : 0; }
  const std::string& token() const { return m_token; }
  bool               failed() const { return m_type == Error; }
  unsigned           line() const  { return m_line; } // of the current token
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

  // gts_file_getc: skips comments up to their newline, and counts lines,
  // positions and braces.
  int getc()
  {
    if (m_type == Error) return EOF;
    int c = raw();
    ++m_curpos;
    if (in(c, kComments))
      while (c != EOF && c != '\n') c = raw();
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

inline bool is_number(Tokenizer::Type t) { return t == Tokenizer::Int || t == Tokenizer::Float; }

// gts_surface_read from f into a new surface of g. vertices and faces, when
// given, get the elements in the order of the file. On error, g is left as
// it was.
std::expected<Surface*, ReadError> read_gts(Geometrium& g, Tokenizer& f, std::string name,
                                            std::vector<VertexId>* vertices = nullptr,
                                            std::vector<FaceId>* faces = nullptr);

} // namespace gts::detail
