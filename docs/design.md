# Design of gts-cxx

- `gts::Geometrium` owns the vertices, edges and faces, and their attributes.
  Surfaces are sets of faces inside it, and a face can belong to several
  surfaces. One thread writes to a Geometrium at a time, or any number read
  it; separate Geometria are independent.
- `Vertex`, `Edge` and `Face` are final, trivially copyable structs of 32, 24
  and 32 bytes, linked by typed 32-bit indices. A vertex stores the last edge
  of the ring of edges around it, and each edge the next edge around each of
  its vertices; faces around an edge are linked the same way. Rings list the
  newest element first, as GTS's lists do.
- A face that leaves its last surface is unlinked, with any edge left without
  faces and any vertex left without edges, as GTS destroys them.
  `KeepOrphans` keeps them, as GTS's `gts_allow_floating_*` did. Unlinked
  elements are marked dead; `collect()` reclaims the storage, renumbers
  everything, and tells subscribers through `on_renumber()`.
- Attributes are arrays parallel to the element arrays, with hooks given as
  lambdas for making a value when an element is made from others: between
  two (`interpolate`), as a piece of one (`split`), or from two that become
  one (`merge`). They take the place of GTS's subclassing of vertices, edges
  and faces. A hook that is not given comes from `DefaultHooks<T>`, which
  interpolates floating-point and other affine values linearly, copies
  anything else, and can be specialized for a type.
- `split_edge()` and `split_face4()` are the primitives of refinement and
  tessellation. They link the new elements in GTS's order, set their
  attributes with the hooks, put new faces in every surface of the face
  they come from, and report each split to `on_split()` subscribers. The
  report includes the cuts, edges made inside a face, which have no parent
  of their own kind and get their attributes from a subscriber.
- `collapse_edge()` is the primitive of coarsening, with `collapsible()`
  and `collapse_folds()` for GTS's tests before a collapse. It merges the
  attributes of the two vertices and of each pair of edges that become
  one, and reports to `on_collapse()` subscribers.
- `swap_edge()` replaces the edge between two faces by the edge between
  their opposite vertices, as GTS's `gts_edge_swap()`.
- Geometry for computation uses the value types `Vec3`, `Segment`,
  `Triangle` and `BBox`.
- `add_edge()` refuses an edge from a vertex to itself, and `add_face()`
  three edges that do not form a triangle.

