// Triangle strips for rendering: the faces of a surface in sequences where
// each face shares an edge with the next.
//
// Copyright (C) 1999-2003 Wagner Toledo Correa, Stéphane Popinet (GTS: stripe.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-11.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "gts/Geometrium.h"

#include <vector>

namespace gts {

// gts_surface_strip: every face of s in one strip. A strip starts from the
// face with the fewest unused faces within two steps of it, takes the
// neighbour with the fewest, and grows forward and backward with faces
// that turn left and right in turn, so that a strip of n faces is drawn
// with n + 2 vertices. Throws std::logic_error where a face of another
// surface is a neighbour, where GTS failed an assertion.
std::vector<std::vector<FaceId>> strips(const Surface& s);

} // namespace gts
