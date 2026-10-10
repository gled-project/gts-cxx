// Reading and writing surfaces in GTS's .gts format.
//
// The format: a line "nv ne nf", optionally followed by four class names;
// then nv lines of vertex coordinates "x y z", ne lines of edges "v1 v2"
// with 1-based vertex numbers, and nf lines of faces "e1 e2 e3" with 1-based
// edge numbers. '#' and '!' start comments. Extra tokens on a line are
// ignored. A vertex class named "GtsVertexBinary" means the coordinates are
// binary doubles in the machine's byte order.
//
// Copyright (C) 1999 Stéphane Popinet (GTS: surface.c, misc.c, point.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "gts/Geometrium.h"

#include <expected>
#include <filesystem>
#include <iosfwd>
#include <string>

namespace gts {

struct ReadError
{
  unsigned    line = 0, pos = 0; // 1-based, as GTS reports them
  std::string message;
};

// Reads one surface into g as a new surface. On error, g is left as it was.
std::expected<Surface*, ReadError> read_gts(Geometrium& g, std::istream& in, std::string name = {});
std::expected<Surface*, ReadError> read_gts(Geometrium& g, const std::filesystem::path& path);

struct WriteOptions
{
  // Significant digits of coordinates. GTS writes 10, which loses
  // precision; 17 is exact for doubles.
  int precision = 10;
};

// Writes s with its vertices and edges in index order and faces in the
// surface's order, numbered from 1.
void write_gts(const Surface& s, std::ostream& out, WriteOptions o = {});
bool write_gts(const Surface& s, const std::filesystem::path& path, WriteOptions o = {});

} // namespace gts
