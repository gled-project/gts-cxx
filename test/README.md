# Testing against the C library

`make check` runs the unit tests, `test/unit/test-core.cxx`, under
AddressSanitizer. `make compare` runs the comparison with the C library
described here; it needs the C GTS library and glib.

`test/compare/gts-ref` runs one GTS operation through the C library and prints
a summary of the result that does not depend on the order of vertices, edges
and faces: counts, area, volume, closedness, manifoldness, orientability,
self-intersection, bounding box, and ranges of face quality, face area, edge
length and edge angle. With `-o` it also writes the resulting surface. The
operations and their parameters follow how Gled calls GTS.

`test/compare/gts-cxx` does the same through the C++ library and prints the
same summary; operations not ported yet are reported as unsupported.
`test/compare/run-cases.sh` runs every case in `test/compare/cases.txt`,
every surface in `test/data/gts-boolean` through `stats`, and every pair in
`test/data/gts-boolean/tests` through union, intersection and difference,
through either program. `compare.sh` runs both and compares the summaries
with `compare.py`: integers exactly, floating-point values to a relative
tolerance (default 1e-9). It then reads every surface the C++ library wrote
back through the C library, to check the writer.

All 307 cases run through both libraries: the 240 of the port of what
Gled calls, 21 for the distances, the inside test, curvature, merging and
the centres, and 46 for what was ported after: Delaunay triangulation, the
iso-surface slices, triangle strips, progressive and hierarchical surfaces,
and graphs. The 21 print a line per vertex or sample point, keyed by its
coordinates, and read surfaces the C library wrote (`@REF@` in
`cases.txt`), so that both programs measure the same input. The summaries
include whether the surface intersects itself. On 2026-10-11:

- 185 cases agree with the C library in every value to 1e-9, and so do
  the surfaces written by the C++ library when read by the C one.
- 15 are in `test/compare/fixed.txt`: gts-cxx fixes an error of GTS there,
  or departs from it on purpose (see `QUIRKS.md`), so the results differ,
  or the C library crashes. `compare.py` runs the check given for each
  instead: the iso-surfaces of the bounded and the body-centred tetrahedra
  must be closed, manifold and orientable, also when the C library reads
  them back; those of the cartesian method and of the slices may differ
  only in the sign of the volume, as their faces are oriented the other
  way; in the curvature of the C library's bounded torus, only vertices
  where gts-cxx gives no Gaussian curvature may differ; the hierarchical
  surfaces, which crash GTS, must succeed.
- 16 are in `test/compare/ties.txt`: GTS chooses between equal costs or
  candidates in the order of its hash tables, which follows memory
  addresses. `compare.sh` checks instead that the two surfaces lie within
  a tenth of the mean edge length of each other, and `compare.py` that the
  lines listed with the case agree: for strips and partitions, that every
  face or node is taken once.
- 1 is self-tied: a flat surface, whose self-intersection Simulation of
  Simplicity decides; gts-cxx counts those ties (`self_sos_ties`).
- 82 are degenerate booleans, and 8 crash in the C library: six
  booleans; the refinement of `crossing.gts`, which does not end (see
  `QUIRKS.md`); and reading back the progressive surface of an open
  surface, an error gts-cxx keeps. gts-cxx reports the errors. See below.

82 boolean cases are degenerate: the two surfaces touch or have coplanar
faces, and Simulation of Simplicity decides the exactly degenerate
orientation tests. GTS breaks those ties by memory address, so its own
results change with the program's allocations; gts-cxx breaks them by
vertex index. `compare.py` treats a boolean case with any such tie between
the two surfaces as degenerate, and checks that the result of gts-cxx is
closed, manifold and orientable. All 82 are.

The C library aborts in some degenerate boolean cases, and which ones
depends on the build of `gts-ref`: three builds on 2026-10-09 aborted in
15, 3 and 6 cases. `horse4 horse5` aborted in all three, and gts-cxx
reports an error there; so it does for `t19 t20`. Where `a23`, `a24` and
`a25` abort, gts-cxx gives valid results. Of the cases where the C library
reports an error, gts-cxx succeeds in 9, the touching cubes.

GTS keeps faces in hash tables keyed by pointer, so its results depend on
memory addresses. In two ordinary runs, 26 of 240 summaries differed:

- means and sums in the last digit, from the order of summation;
- coarsening of a symmetric sphere, where many costs tie, by 3e-5 in area;
- self-intersection of one boolean result;
- whether six degenerate booleans abort.

Run under `setarch -R`, which turns off address space randomisation, the
reference repeats exactly: two runs gave identical summaries for all 240 cases
and identical surfaces for all 148 that wrote one. `run-cases.sh` therefore
uses it. gts-cxx orders things differently from the C library, and the
comparison allows for that.

Results of the reference for the 240 cases of the port, on 2026-10-09
(Fedora's GTS 0.7.6, snapshot 2012-11-30, glib 2.84.4):

- 150 cases succeed.
- 75 fail with an error. All of them are booleans whose intersection curve is
  not closed; GTS's own `set` example reports the same. They are the
  degenerate single-triangle pairs (`t1 t2` and so on), which GTS's test suite
  uses only to compute the intersection curve.
- 15 abort on a `g_assert` in `boolean.c`, in `new_ear` or
  `connect_interior_loop`. All are on pairs that GTS marks `XFAIL`
  (`horse4 horse5`) or that Debian disabled (`t19 t20`, `a23 b23`,
  `a24 b24`, `a25 b25`).
- One case, `sphere-0`, trips GTS's argument check
  (`geodesation_order != 0`) and prints a glib CRITICAL message.

