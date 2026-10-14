# Height-field test inputs

Made on 2026-10-13 with Python's `random.Random(20261013)`, from a terrain
of four Gaussian hills and valleys on a slope, with a ripple and noise of
1% of the hills' height:

- `terrain-120x90.pgm`: the terrain on a 120 x 90 grid, as an ASCII PGM
  (P2) with heights from 0 to 65535, so that equal errors are rare.
- `terrain8-64x48.pgm`: the same on a 64 x 48 grid, as a binary PGM (P5)
  with heights from 0 to 255, where many errors are equal.
- `points-2000.xyz`: 2000 random points in the unit square with the
  terrain's height, one `x y z` per line, as GTS's `happrox -f` reads them.
