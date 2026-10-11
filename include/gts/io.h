// Reading and writing surfaces in GTS's .gts format, and writing them for
// Geomview (OOGL) and VTK.
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
#include <functional>
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

// A colour for the Geomview writers, GtsColor.
struct Rgb { double r = 0, g = 0, b = 0; };

// Colours for the Geomview writers, taking the place of GTS's colour
// method of a class. A vertex colour makes the object COFF.
struct OoglColors
{
  std::function<Rgb(VertexId)> vertex;
  std::function<Rgb(FaceId)>   face;
  std::function<Rgb(EdgeId)>   edge;
};

// gts_surface_write_oogl: s as a Geomview OFF object, vertices in index
// order, numbered from 0, coordinates with 6 significant digits.
void write_oogl(const Surface& s, std::ostream& out, const OoglColors& c = {});
// gts_surface_write_oogl_boundary: the boundary edges of s as a Geomview
// LIST of VECT lines.
void write_oogl_boundary(const Surface& s, std::ostream& out, const OoglColors& c = {});
// gts_surface_write_vtk: s as legacy VTK polygon data, vertices in index
// order.
void write_vtk(const Surface& s, std::ostream& out);

} // namespace gts
