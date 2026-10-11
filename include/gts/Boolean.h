// Booleans: the intersection of two surfaces, and their union, intersection
// and differences, as GTS's GtsSurfaceInter; and the test for a surface
// that intersects itself.
//
// SurfaceInter finds the curve where the faces of two surfaces cross, cuts
// the faces it crosses along the curve, and triangulates the pieces. The
// faces of the result of boolean() are faces of the two surfaces and the
// new pieces, shared with them. Both surfaces must be in one Geometrium.
//
// Copyright (C) 1999 Stéphane Popinet (GTS: boolean.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "gts/Geometrium.h"

#include <functional>
#include <vector>

namespace gts {

// GtsBooleanOperation: the part of one surface outside or inside the other.
enum class BooleanOp { A_out_B, A_in_B, B_out_A, B_in_A };

class SurfaceInter
{
public:
  // gts_surface_inter_new.
  SurfaceInter(Surface& a, Surface& b);
  // Removes the two surfaces of pieces; pieces in no other surface are
  // unlinked.
  ~SurfaceInter();
  SurfaceInter(const SurfaceInter&) = delete;
  SurfaceInter& operator=(const SurfaceInter&) = delete;

  // The edges of the intersection curve, newest first as in GTS.
  const std::vector<EdgeId>& edges() const { return m_edges; }

  // The pieces of the faces of a and of b that the curve cuts.
  const Surface& pieces_a() const { return *m_pieces[0]; }
  const Surface& pieces_b() const { return *m_pieces[1]; }

  // How many orientation tests of vertices of both surfaces were exactly
  // degenerate, as where the surfaces touch or have coplanar faces, so that
  // Simulation of Simplicity decided them. When there are any, the curve
  // and the pieces depend on the ranking of the points; see QUIRKS.md.
  std::size_t sos_ties() const { return m_sos_ties; }

  // gts_surface_inter_check: whether the curve is consistently oriented
  // and the pieces join the uncut faces, and whether the curve is closed.
  struct Check { bool ok = true; bool closed = false; };
  Check check() const;

  // gts_surface_inter_boolean: adds to result the faces of the part op
  // selects.
  void boolean(Surface& result, BooleanOp op) const;

private:
  Geometrium&         m_geo;
  Surface*            m_surface[2];
  Surface*            m_pieces[2];
  std::vector<EdgeId> m_edges;
  std::size_t         m_sos_ties = 0;
};

// gts_surface_foreach_intersecting_face: calls f(face1, face2) for each pair
// of faces of s that cross, where an edge of the second passes through the
// first and shares no vertex with it. Each pair may come twice. Returns
// whether there is any.
bool foreach_intersecting_face(const Surface& s, const std::function<void(FaceId, FaceId)>& f);

// gts_surface_is_self_intersecting, as a yes or no. With sos_ties, also
// the number of orientation tests that were exactly zero, so that
// Simulation of Simplicity decided them, as in a flat surface.
bool is_self_intersecting(const Surface& s, std::size_t* sos_ties = nullptr);

} // namespace gts
