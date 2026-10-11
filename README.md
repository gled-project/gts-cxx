# gts-cxx

A rewrite in C++23 of GTS, the GNU Triangulated Surface Library by Stéphane
Popinet (<https://gts.sourceforge.net/>), without glib and without any other
dependency.

**Status: all of GTS is ported**, except what C++ replaces by design (see
below). The library has its core data structures, surface statistics,
reading and writing of `.gts` files, and writing for Geomview and VTK; the
edge heap, tessellation, refinement by edge splitting, generated spheres,
coarsening by edge collapse with the volume-optimized vertex and cost;
iso-surfaces by marching cubes, on cartesian and other structured grids, and
by tetrahedra; Shewchuk's robust predicates, bounding-box trees, booleans,
the test for self-intersection and the cluster grid. It has the distances
between points, faces, edges and surfaces, the test for a point inside a
surface, discrete curvature, the kd-tree, the merging of vertices and edges,
and the removal of duplicate faces; the neighbourhoods of vertices, edges
and faces, connected components, contact vertices and edges; constrained
Delaunay triangulation in the plane with conforming and refinement;
triangle strips; 4 x 4 matrices; progressive and hierarchical surfaces; and
graphs with their partitions. [Gled](https://github.com/gled-project/gled)
uses it through GledGTS.

The port keeps GTS's algorithms and its floating-point arithmetic, so that
its results can be compared with the C library's; redesigns come later.
Errors of GTS are fixed one by one, and the comparison checks the cases
they change.

- `docs/design.md`: the data structures and the primitives.
- `QUIRKS.md`: the errors of GTS that the port fixes and those it still
  has, odd behaviour it keeps, and where it departs from GTS.
- `test/README.md`: how the port is tested against the C library, and the
  results.

## License

LGPL-2.1-or-later; see `COPYING`. GTS is distributed under the GNU Library
General Public License, version 2 or later, which allows this. Files derived
from GTS keep the copyright lines of the original files and state that they
were rewritten in C++.

## Building

```sh
make            # lib/libgts-cxx.so and the test programs
make lib        # the library only
make install PREFIX=/opt/gled   # include/gts/*.h and lib/libgts-cxx.so
make check      # the unit tests, under AddressSanitizer
make compare    # all cases through the C and the C++ library, compared
make ref-check  # all cases through the C reference only
```

A plain Makefile, with `-std=c++23`, tested with GCC 15 and Clang 20.
GCC 15 is the first whose library has `<flat_set>`. gled-builder clones the
repository and runs `make lib` and `make install`. Only
the reference program needs anything installed: the C GTS library and glib,
found through `pkg-config`. The library is `libgts-cxx`, linked with
`-lgts-cxx`, and its namespace is `gts`. Headers are included as
`<gts/...>`; the C library installs `<gts.h>` and `<gtsconfig.h>`, so the
two do not collide.

## Not ported

Left out by design, as C++ replaces them: the object system (`object.c`)
and its classes of vertices, edges and faces, whose data attributes carry
instead (`GtsVertexNormal`, `GtsColorVertex`, the named objects of
`named.c`); containers (`container.c`), the generic heap (`heap.c`) and the
queue (`fifo.c`); the `GtsFile` tokenizer, of which `io` keeps what it
reads; and the vector and point helpers, which `Vec3` and `BBox` cover.
Not ported either: `GtsListFace`, the faces that carry lists of points
through Delaunay insertions for GTS's example `happrox.c`; reading a
progressive surface one split at a time, for streaming (`PSurface::read()`
reads it all); and GtsEHeap's randomized mode.

## Layout

```
COPYING                  LGPL-2.1
QUIRKS.md                fixed and remaining errors of GTS, choices, departures
docs/design.md           the data structures and the primitives
Makefile
include/gts/
  Geometry.h             value types and measures: Vec3, Segment, Triangle, BBox
  Elements.h             Vertex, Edge, Face, their typed indices and flags
  Attribute.h            attributes and their hooks
  Geometrium.h           Geometrium, Surface, rings, collection, renumbering
  Topology.h             neighbourhoods, fans, contacts, components, traversal
  EHeap.h                the edge heap, as GtsEHeap
  Refine.h               tessellation, refinement, spheres
  Coarsen.h              coarsening, the volume-optimized vertex and cost
  Progressive.h          vertex splits, progressive and hierarchical surfaces
  Iso.h                  iso-surfaces: marching cubes and tetrahedra, slices
  Predicates.h           Shewchuk's orient2d, orient3d, incircle, insphere; SoS
  BBTree.h               bounding-box trees of faces or edges, box overlaps
  Boolean.h              intersection of surfaces, booleans, self-intersection
  Cluster.h              the cluster grid: simplification by vertex clustering
  Distance.h             distances, closest points, inside test, surface distance
  Curvature.h            discrete mean and Gaussian curvature, principal directions
  KdTree.h               the kd-tree of vertices, for merging them
  Delaunay.h             constrained Delaunay triangulation, conforming, refinement
  Strip.h                triangle strips
  Matrix.h               4 x 4 matrices and transformations
  Graph.h                graphs, progressive graphs, partitions
  Stats.h                area, volume, closed/manifold/orientable, boundary, centres
  io.h                   .gts reading and writing, Geomview and VTK writing
  all.h                  all of the above
src/                     the implementation
test/unit/test-core.cxx  unit tests
test/README.md           testing against the C library
test/compare/
  gts-ref.cxx            the reference program, linked with the C GTS library
  gts-cxx.cxx            its twin, linked with this library
  cases.txt              the cases
  ties.txt               cases decided by GTS's memory order, compared by distance
  fixed.txt              cases where gts-cxx fixes GTS, with the check for each
  run-cases.sh           runs all cases through one program
  run-ref.sh             runs them through gts-ref
  compare.sh, compare.py runs both and compares; output in test/compare/out/
test/data/gts-boolean/   boolean test surfaces from GTS, see its README.md
test/data/delaunay/      Delaunay test points and constraints, see its README.md
test/data/graph/         a graph in the Jostle format, see its README.md
spike/                   (branch spike/core) the three designs compared
```
