// The elements of a mesh: vertices, edges and faces, and their indices.
//
// The elements are plain structs, linked by 32-bit indices into the arrays of
// their Geometrium. The links are rings stored in the elements themselves:
//
// - A vertex stores the last edge of the ring of edges around it. Each edge
//   stores, for each of its two vertices, the next edge around that vertex.
// - An edge stores the last face of the ring of faces around it. Each face
//   stores, for each of its three edges, the next face around that edge.
//
// A walk starts at the successor of the last element and ends with the last.
// A new element is inserted after the last, which makes it the first, so a
// ring lists the newest element first, in the order of GTS's lists. Any
// number of faces can share an edge, so non-manifold surfaces are possible.
//
// Copyright (C) 2026 Matevz Tadel
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "gts/Geometry.h"

#include <array>
#include <compare>
#include <cstdint>
#include <type_traits>

namespace gts {

using Index = std::uint32_t;
inline constexpr Index kNone = ~Index(0);

// A typed index, so that a vertex index cannot be passed as an edge index.
template <class Tag>
struct Id
{
  Index i = kNone;

  constexpr Id() = default;
  constexpr explicit Id(Index x) : i(x) {}
  constexpr bool valid() const { return i != kNone; }
  friend constexpr auto operator<=>(Id, Id) = default;
};

using VertexId = Id<struct VertexTag>;
using EdgeId   = Id<struct EdgeTag>;
using FaceId   = Id<struct FaceTag>;

// Flags. The low byte is the library's; users take bits from kUserFlag up.
enum : std::uint32_t
{
  kDead     = 1u << 0, // removed; the storage is reclaimed by collect()
  kUserFlag = 1u << 8
};

struct Vertex final
{
  Vec3          p;
  EdgeId        edge;      // the last edge of the ring around this vertex
  std::uint32_t flags = 0;
};

struct Edge final
{
  std::array<VertexId, 2> v;
  std::array<EdgeId, 2>   next;  // next[k]: the next edge around v[k]
  FaceId                  face;  // the last face of the ring around this edge
  std::uint32_t           flags = 0;
};

struct Face final
{
  std::array<EdgeId, 3> e;
  std::array<FaceId, 3> next;     // next[k]: the next face around e[k]
  std::uint16_t         surfaces = 0; // number of surfaces containing the face
  std::uint32_t         flags = 0;
};

static_assert(std::is_trivially_copyable_v<Vertex> && std::is_standard_layout_v<Vertex>);
static_assert(std::is_trivially_copyable_v<Edge>   && std::is_standard_layout_v<Edge>);
static_assert(std::is_trivially_copyable_v<Face>   && std::is_standard_layout_v<Face>);

} // namespace gts
