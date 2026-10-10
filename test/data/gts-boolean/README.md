# Boolean test surfaces from GTS

These surfaces and the `tests` list come unchanged from `test/boolean/` of
GTS, the GNU Triangulated Surface Library, as of the head of Stéphane
Popinet's gts-stable darcs repository (tag of 2012-11-30). GTS is
distributed under the GNU Library General Public License, version 2 or
later.

The comments in `tests` name the contributors of some surfaces: T. Binder,
B. Petit, Darie Aron, Wayne Gramlich and Romain Behar.

Each line of `tests` names two surfaces. GTS's test script runs the union,
intersection and difference of each pair. `XFAIL` marks a pair that GTS
expected to fail. Debian and the melowntech fork disabled further pairs that
fail on some platforms: `t19 t20`, `t31 t32`, `cube cube1`, `cube cube2`,
`cube cube4`, `cube6 cube7`, `sponge cutter`, `sphere.gts sphere.gts`,
`a23`–`a25`, `gts-crash.gts gts-crash1.gts`, the `buildings` pairs and the
`gerris` pairs.
