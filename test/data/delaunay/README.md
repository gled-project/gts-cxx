# Delaunay test inputs

The format is the one GTS's example `delaunay.c` reads: a line `nv ne nf`,
nv vertices `x y z`, then ne constraints as pairs of 1-based vertex numbers.
Faces are ignored.

`squarehole.gts`, `too_close.gts` and `two_segments.gts` come unchanged from
`test/delaunay/` of GTS, the GNU Triangulated Surface Library, as of the head
of Stéphane Popinet's gts-stable darcs repository (tag of 2012-11-30). GTS is
distributed under the GNU Library General Public License, version 2 or
later.

The others were made on 2026-10-11 with Python's `random.Random(20261011)`:

- `random-100.gts`, `random-2000.gts`: points in the unit square, with a
  random height.
- `cartesian-10.gts`: a 10 x 10 grid turned by 0.1 radians, in the order of
  GTS's test `cartesian.c`; its points are cocircular in fours.
- `ring.gts`: a disc of 24 sides, counterclockwise, round a square hole,
  clockwise, as constraints, and 150 points between them.
- `crossing.gts`: 300 points in the unit square and two constraints that
  cross each other.
