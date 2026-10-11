# Quirks

gts-cxx 0.1.0 reproduced GTS exactly, so that its results could be compared
with the C library's case by case. Where GTS behaves oddly, the port kept
the behaviour and listed it here, with what could be done instead. Where the
port had to depart from GTS, or chose to, the departure is listed too.

The entries are in four groups. Fixed: errors of GTS that gts-cxx no longer
has. Errors: errors of GTS that gts-cxx still has. Choices: behaviour that
is consistent but may be changed, to be decided. Kept: behaviour that stays.
Nothing under Choices is decided until it is taken up.

## Fixed

`test/compare/fixed.txt` lists the comparison cases whose results these
fixes change, with the check each must pass instead of equality with GTS.

- Reading `.gts` files: a comment at the end of a line swallowed the
  newline and the line after it, and the C library fails on such files. The
  comment now ends the line.
- `isosurface_tetra_bounded()` is meant to close the surface where it leaves
  the grid by adding a layer of points with f - iso = -1 round it. GTS's
  `copy_to_bounded()` built that layer with two indexing errors. Its loop
  over columns stopped two columns early, so the last two columns of the
  function's values were dropped. It wrote the last border column at index
  ny instead of nx + 1, so two or three columns of the enlarged slice were
  never written, and when ny exceeded nx + 1 it wrote past the slice. The
  unwritten columns kept what malloc returned. In GledGTS a 24 x 20 x 18
  grid gave 13255 and 12892 vertices in two runs, the second with a NaN
  area, and valgrind reported GTS reading uninitialised memory. gts-cxx now
  copies every value into the enlarged slice inside a complete border, and
  the surfaces close.
- `isosurface_tetra_bcl()` subtracted iso from the first slice twice and
  never from the second, so the first layer of cubes was wrong and the
  surfaces of the sphere and the blobs were not closed. It now subtracts iso
  from each slice once.
- `SurfaceInter::check()` stopped looking at the neighbours of a segment
  once a check had failed, so all later segments counted as not closed. It
  now counts them, and `closed` says whether the curve is closed.
- `refine()` computed the cost of each new edge before the edge had faces,
  so with a cost that depends on the faces, such as GledGTS's angle cost,
  every new edge got the default cost. The new edges are now costed after
  the split.
- `quality_stats()` took the edge angles between every pair of faces on an
  edge, also faces of other surfaces. It now takes only the faces of the
  surface.
- `boundary_distance()` and `distance()` counted as boundary an edge with
  one face in any surface, as GTS passed no surface to
  `gts_edge_is_boundary()`. They now count the faces of the surface
  measured.
- `is_inside()` ended its ray at x_max + |x_max| / 10 of the tree's box, so
  when x_max was 0 the ray ended on the box. It now ends a tenth of the
  larger of |x_max| and the box's width beyond the box.
- `gaussian_curvature()` divided by the area of the vertex's region without
  checking it. It now gives none where the area is not positive, as
  `mean_curvature_normal()` does. The C library's bounded torus has such
  vertices, where GTS gave infinity.
- `isosurface_tetra_bcl()` built, for each cube, the tetrahedra towards
  the cubes above it and beside it, but skipped the last layer of cubes in
  z and the last column in x and in y altogether. Their tetrahedra towards
  the neighbours they do have were never built, and a surface that crossed
  them was left open. Every cube is now visited, and each group of
  tetrahedra is built where its neighbour exists.
- `tessellate()` left a face of another surface on a split edge on the old
  edge, beside the two halves, so the two surfaces no longer shared the
  seam. `Geometrium::split_remaining()` now splits such faces in two at the
  same vertex, with the same halves.
- When a collapse left an edge without faces (only with duplicate faces in
  the input), GTS gave up the new vertex if that edge was the last in the
  vertex's list, even when the vertex still had other edges. The new
  vertex is now given up only when it has no edges left.

## Errors

None known.

## Choices

- Refinement and coarsening fill their edge heap in the order of GTS's face
  hash table, which follows memory addresses. Where costs tie, the result
  then changes from run to run (see `test/compare/ties.txt`). gts-cxx fills
  the heap in face index order, so it is repeatable, but it chooses
  differently from GTS. A tie-break by geometry, such as by the edge's
  midpoint, would make the result independent of the order of the input
  too.
- `isosurface_tetra()` and its two variants orient the faces towards
  f < iso. `isosurface_cartesian()`, and GTS's documentation of all four,
  orient them towards f > iso. A sphere from the tetrahedra has negative
  volume. GledGTS's `GTSIsoMaker` inverts the tetra surfaces
  (`bInvertTetra`).
- The layer of `isosurface_tetra_bounded()` marks as outside the side where
  f - iso < 0, which is the inside for a function that grows outwards, such
  as the squared distance from a centre. The surface of such a function
  then gets a second sheet round the grid. The layer lies a grid step
  beyond the grid.
- `tessellate()` without a splitting function projects every new vertex
  onto the unit sphere, as GTS's `gts_surface_tessellate` does. A cube
  tessellated this way becomes a rough sphere. The default suits only
  `generate_sphere()`.
- `gts_matrix_compatible_row()` drops a constraint row only when it is
  exactly zero. melowntech's GTS changed the test to `|row|^2 < 1e-12`
  against instabilities when the row is nearly zero. That threshold is
  absolute, and for a mesh a few centimetres across measured in metres,
  valid rows are of that size. A threshold relative to the scale of the
  faces around the edge would serve both.
- `sample()` steps through a plane by adding dx and dy, as GledGTS's
  `GTSIsoMaker` did, and the cartesian method steps z by adding and
  subtracting dz. The tetrahedra methods compute positions as
  x0 + dx * (i + t). The two kinds of methods put the same grid points at
  positions that differ in the last bits. Computing every coordinate as
  x0 + i * dx would be exact to one rounding.
- `principal_curvatures()` is given half the length of the mean curvature
  normal, so the mean curvature has no sign, and a concave region gets the
  same curvatures as the convex one.

## Kept

- GtsEHeap's randomized mode is not ported.
- The fold test (`collapse_folds()`) finds faces that fold back onto each
  other across an edge, and faces 1e4 times larger than a neighbour. It does
  not find faces that turn over without folding onto a neighbour, nor faces
  that pass through other parts of the surface.
- After a collapse, the costs of the edges within two edges of the new
  vertex are recomputed. An edge there that was not in the heap, such as an
  edge of another surface, is put in. GTS did so after a CRITICAL message.
- The neighbourhood lists of the volume-optimized vertex are built as GTS
  built them, with a linear search for each element. Vertices with
  hundreds of edges, as at the caps of `1.gts`, make this quadratic.
  gts-cxx is still about 4 times faster than GTS on that case.
- `isosurface_tetra_bcl()` covers the grid only between the centres of the
  outermost cubes, half a grid step inside its boundary, as its tetrahedra
  join cube centres. A surface that comes closer to the boundary is left
  open there, where the other methods close it. GledGTS's `GTSTorus` pads
  its grid by 5% of the radii, which is enough unless the step exceeds a
  tenth of the minor radius. With a minor radius of 0.8 and a step of 0.1,
  the outermost centres in z are at 0.7906 and the torus reaches 0.8: the
  bcl torus has 2016 boundary edges, all at |z| > 0.79, where the
  cartesian and tetra tori are closed. The `iso_torus.C` demo's torus
  (1.5 and 0.05) closes with every method.
- GTS's `get_vertex()` prints to stderr when it finds a crossing in an
  unexpected slice. That cannot happen with its own tetrahedra, and the port
  leaves the print out.
- The faces that the intersection curve cuts are triangulated by ear
  clipping, starting from the first segment of lists whose order follows
  the bounding-box tree. GTS builds that tree from its faces in the order of
  its hash table, so the same polygon can get different diagonals in GTS
  and in gts-cxx. In `horse5 horse6`, 5 faces of a few thousand differ;
  every vertex agrees. A triangulation that does not depend on the order,
  or a constrained Delaunay one, would remove the difference and give
  better triangles.
- `new_ear()` refuses an ear when any edge joins its two outer vertices,
  through `gts_vertices_are_connected()`, not only edges of the face being
  cut. The outer vertices are corners of the face or points of the curve,
  so such an edge is an edge of the same surface, a segment of the curve,
  or, where the two surfaces have coplanar faces, a diagonal of the other
  surface's pieces; refusing is right for all of them. In the 261
  comparison cases, every edge that refused an ear was of the surface
  being cut or had no faces yet.
- `gts_surface_inter_new()` took whether each surface is open, and did not
  use it. `SurfaceInter` does not take it.
- `tribox3.c`, the box-triangle overlap test, is not ported: neither the
  booleans nor the self-intersection test use it.
- `Geometrium::normal()` uses `gts_triangle_normal`'s own choice of
  vertices, which gives the same orientation as `vertices_of()` but not
  always the same first vertex, so that areas agree with GTS to the bit.
- `BBTree::point_closest()` gives the square of the distance; GTS
  documented it as the distance.
- Where two leaves are equally close, `point_distance()` and
  `point_closest()` take the first in the tree's order, which follows the
  order of the faces. GTS's order of faces was its hash order.
- `principal_directions()` is arbitrary where the two curvatures are equal,
  as everywhere on a sphere, and the last bit of a sum decides it.

## Departures from GTS

- `add_edge()` refuses an edge from a vertex to itself, and `add_face()`
  three edges that do not form a triangle. GTS accepted both and failed
  later. An edge that has a duplicate is not collapsible, as in GTS, so
  coarsening never makes such an edge.
- `generate_sphere()` throws for order 0, where GTS printed a CRITICAL
  message and left the surface empty.
- `volume_optimized_vertex()` throws when the constraints do not fix the
  position, where GTS aborted on an assertion.
- New faces from `split_edge()` and `split_face4()` join every surface of
  the face they come from. GTS added them only to the surface being
  refined, also for faces of other surfaces on the split edge.
- `isosurface_cartesian()` makes a vertex where a face first uses a
  crossing. GTS made a vertex at every crossing of the grid when it filled a
  slice, and left those no face used without edges.
- `isosurface_tetra_bcl()` keys its vertices by the grid indices of the
  corners. GTS keyed them by strings printed from the indices, in hash
  tables.
- Simulation of Simplicity, which breaks the ties of exactly degenerate
  orientation tests, ranks points by vertex index. GTS ranked them by
  memory address, truncated to 32 bits. Where two surfaces touch, or have
  coplanar faces, GTS's results therefore depend on everything the program
  allocated before. The same library on the same `cube` and `cube1` found
  14 segments of the intersection curve in a small test program and none
  in `gts-ref`. gts-cxx is repeatable. `SurfaceInter::sos_ties()` counts
  the ties between the two surfaces, and `test/compare/compare.py` checks
  the cases with any such tie for validity instead of equality with GTS.
  In all 73 of them that differ, the result of gts-cxx is closed, manifold
  and orientable. Where GTS gave an error for touching cubes (`cube cube1`,
  `cube6 cube7`, `cube9 cube10`), and where it aborted (`a23`, `a24`,
  `a25`), gts-cxx gives valid results.
- Where GTS aborted on an assertion in the booleans (`connect_interior_loop`,
  `new_ear`, `intersect_edges`, `triangle_triangle_orientation`), gts-cxx
  throws `std::runtime_error` or `std::logic_error`, and the surfaces are
  left as they were.
- A vertex where an edge crosses a face gets its attributes interpolated
  along the edge. The pieces of a cut edge get the edge's attributes by the
  split hooks, and the pieces of a cut face the face's. GTS called the
  vertex class's `intersection_attributes` and the face class's
  `attributes`, and gave the pieces of an edge nothing.
- `fast_expansion_sum_zeroelim()` in Shewchuk's predicates read one element
  past the end of an expansion it had used up, and never used the value. It
  no longer reads it. The predicates agree with GTS's on 2 million cases.
- `swap_edge()` gives the two new faces the attributes of the faces they
  replace, by the split hooks. GTS's `gts_edge_swap()` gave them none. It
  throws `std::invalid_argument` where GTS printed a CRITICAL message or
  aborted on an assertion: when the edge does not have exactly two faces in
  the surface, or they are not oriented consistently.
- `is_inside()` breaks ties by Simulation of Simplicity with the point and
  the ray's end ranked after all vertices. GTS ranked them by address, as
  everything else.
- `distance()` gives zero for a range with nothing to measure: no faces in
  the second surface, or no boundary edges. GTS reset the count, extremes
  and sums, and left the mean and standard deviation as they were.
- `principal_directions()` throws `std::invalid_argument` where GTS aborted
  on an assertion, at an edge without exactly two faces in the surface, and
  where it returned without a result, when no normal can be found.
- `merge_vertices()` removes an edge between two vertices it merges,
  together with its faces. GTS kept the edge, joining the vertex to itself,
  and left removing it to the caller, as its example `cleanup.c` does.
- `remove_duplicate_faces()` removes a face from the surface where GTS's
  example `cleanup.c` destroyed it, and counts as duplicates only faces of
  the same surface. With one surface the result is the same.
- `KdTree` keeps vertices with equal coordinates in the order it was given
  them. GTS's comparison never called two coordinates equal, which made
  their order depend on qsort.
- `ClusterGrid` puts a point outside its box in cluster (0, 0, 0) without a
  word. GTS did the same after printing a CRITICAL message.
- `ClusterGrid` sums the points of a cluster apart from its vertex, so
  `update()` can be called again after more triangles are added. GTS summed
  into the vertex's coordinates and divided them in place. The points are
  summed in the order the triangles come, which for a surface follows
  GTS's hash order there, so positions can differ in the last bits.
- `ClusterGrid` takes positions, as GTS did, so the cluster vertices get
  default attribute values. Merging the attributes of the vertices that
  fall into a cluster, by the merge hooks, would carry colours and the like
  through the simplification.
- Elements are not destroyed but marked dead, and `collect()` reclaims them.
  Attributes of dead elements stay readable until then, so that callouts
  can read the attributes of what a split or a collapse removed.
