# gts-cxx

A rewrite in C++23 of GTS, the GNU Triangulated Surface Library by Stéphane
Popinet (<https://gts.sourceforge.net/>), without glib and without any other
dependency.

**Status: every GTS function that
[Gled](https://github.com/gled-project/gled) calls is ported,** and some
more. The library has its core data structures, surface statistics, reading and
writing of `.gts` files, the edge heap, tessellation, refinement by edge
splitting, generated spheres, coarsening by edge collapse with the
volume-optimized vertex and cost, iso-surfaces by marching cubes and by
tetrahedra, Shewchuk's robust predicates, bounding-box trees, booleans, the
test for self-intersection and the cluster grid. It also has the distances
between points, faces, edges and surfaces, the test for a point inside a
surface, discrete curvature, the kd-tree, the merging of vertices and
edges, and the removal of duplicate faces.

The port keeps GTS's algorithms and its floating-point arithmetic, so that
its results can be compared with the C library's; redesigns come later.

- `docs/design.md`: the data structures and the primitives.
- `QUIRKS.md`: where the port keeps odd behaviour of GTS, and where it
  departs from GTS.
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

About 9.8k of GTS's 28.2k lines in `src/*.c` have no counterpart, as
nothing needs them yet. They can be ported the same way when something
does:

| GTS files | lines | what |
|---|---|---|
| `graph.c`, `pgraph.c`, `partition.c` | 3.6k | graphs, progressive graphs, graph partitioning |
| `split.c`, `psurface.c`, `hsurface.c` | 2.7k | progressive and hierarchical surfaces: vertex splits, levels of detail |
| `cdt.c`, `refine.c` | 1.6k | constrained Delaunay triangulation in 2D, and its refinement |
| `stripe.c` | 0.8k | triangle strips for rendering |
| `named.c`, `fifo.c` | 0.4k | named objects, a queue |

Left out by design, as the C++ structures replace them: the object system
(`object.c`), containers (`container.c`), the generic heap (`heap.c`), the
`GtsFile` tokenizer and the general matrix type. Within ported files, some
smaller functions are left out: the in-circle and in-sphere tests,
circumcircle and enclosing triangle (for Delaunay), vertex normals, fans
and contact tests, connected components (`gts_surface_split`), and the
Geomview (OOGL) and VTK writers.

## Layout

```
COPYING                  LGPL-2.1
QUIRKS.md                behaviour kept from GTS to improve later, and departures
docs/design.md           the data structures and the primitives
Makefile
include/gts/
  Geometry.h             value types and measures: Vec3, Segment, Triangle, BBox
  Elements.h             Vertex, Edge, Face, their typed indices and flags
  Attribute.h            attributes and their hooks
  Geometrium.h           Geometrium, Surface, rings, collection, renumbering
  EHeap.h                the edge heap, as GtsEHeap
  Refine.h               tessellation, refinement, spheres
  Coarsen.h              coarsening, the volume-optimized vertex and cost
  Iso.h                  iso-surfaces: marching cubes and tetrahedra
  Predicates.h           Shewchuk's orient2d and orient3d, Simulation of Simplicity
  BBTree.h               bounding-box trees of faces or edges
  Boolean.h              intersection of surfaces, booleans, self-intersection
  Cluster.h              the cluster grid: simplification by vertex clustering
  Distance.h             distances, closest points, inside test, surface distance
  Curvature.h            discrete mean and Gaussian curvature, principal directions
  KdTree.h               the kd-tree of vertices, for merging them
  Stats.h                area, volume, closed/manifold/orientable, boundary, centres
  io.h                   .gts reading and writing
  all.h                  all of the above
src/                     the implementation
test/unit/test-core.cxx  unit tests of the core
test/README.md           testing against the C library
test/compare/
  gts-ref.cxx            the reference program, linked with the C GTS library
  gts-cxx.cxx            its twin, linked with this library
  cases.txt              the cases
  ties.txt               cases decided by GTS's memory order, compared by distance
  run-cases.sh           runs all cases through one program
  run-ref.sh             runs them through gts-ref
  compare.sh, compare.py runs both and compares; output in test/compare/out/
test/data/gts-boolean/   boolean test surfaces from GTS, see its README.md
spike/                   (branch spike/core) the three designs compared
```
