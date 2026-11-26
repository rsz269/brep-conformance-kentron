# brep-conformance-kentron

The [brep-conformance](https://github.com/springcraft-ai/brep-conformance) adapter for kentron, a
boundary-representation CAD kernel in Rust: it puts kentron behind the suite's protocol, so that the
suite can score it.

## Who maintains this

brep-conformance scores a B-rep kernel through an adapter: a small program that answers the suite's
requests with the kernel's own calls.

kentron's maintainers maintain this adapter, beside kentron itself at
[rsz269/kentron](https://github.com/rsz269/kentron).

- **It adds nothing to the kernel.** The adapter is pass-through, as the suite's
  [contract](https://github.com/springcraft-ai/brep-conformance/blob/main/doc/contract.md#what-an-adapter-may-do)
  asks: one call of the kernel for each function, and anything the kernel has no one call for is
  listed as unsupported. What is scored is kentron, not this code.
- **The suite decides what is scored.** An official run uses only a commit of this adapter that the
  suite's maintainers have reviewed as pass-through.

## Overview

The adapter answers the requests of the suite's
[api.md](https://github.com/springcraft-ai/brep-conformance/blob/main/doc/api.md#calling-a-kernel)
on stdin and stdout, through `serve` in the `brep-protocol` crate. [`src/main.rs`](src/main.rs)
holds every solid in one kentron `Model` and answers each function with kentron's one call for it:

| Function                                       | kentron                                                               |
| ---------------------------------------------- | --------------------------------------------------------------------- |
| `box`, `cylinder`, `sphere`, `torus`, `cone`   | `primitives::Cuboid`, `Cylinder`, `Sphere`, `Torus`, `Cone` (`build`) |
| `slab`                                         | `primitives::Slab` (`build`)                                          |
| `extrude`, `revolve`                           | `sweep::extrude`, `sweep::revolve`                                    |
| `screw`                                        | `sweep::screw_extrude`, lead `2π·rise`, height `rise·angle`           |
| `loft`, `loft_through`                         | `loft::loft`, ruled and smooth, closed into a solid                   |
| `sweep`                                        | `sweep::pipe` along the path's lines and arcs                         |
| `sweep_path`                                   | `sweep::pipe_with`, `Guide::Discrete` (least rotations), its `twist`  |
| `union`, `intersection`, `difference`          | `boolean::union`, `intersection`, `difference`                        |
| `union_all`, `difference_all`                  | `boolean::combine` (`Combination::Union`, `Difference`)               |
| `fillet`, `fillet_radii`                       | `blend::fillet`, `blend::fillet_edges`                                |
| `fillet_variable`                              | `blend::fillet_linear`, the radii ordered by the edge's vertices      |
| `fillet_faces`, `full_round`                   | `blend::fillet_faces`, `blend::full_round`                            |
| `fillet_setback`                               | `blend::fillet_setback`                                               |
| `fillet_setback_g2`                            | `blend::fillet_setback_g2`                                            |
| `fillet_conic`, `fillet_g2`                    | `blend::fillet_conic`, `blend::fillet_g2`                             |
| `chamfer`, `chamfer_two`, `chamfer_distances`  | `blend::chamfer`, `blend::chamfer_edges`                              |
| `offset`                                       | `offset::shape`, joined by arcs                                       |
| `offset_faces`                                 | `offset::shape_with_faces`, the other faces by 0, joined by meeting   |
| `delete_faces`                                 | `defeature::remove_faces`                                             |
| `move_faces`                                   | `draft::move_faces`                                                   |
| `thicken`                                      | `offset::thicken`                                                     |
| `shell`                                        | `offset::thick_solid` inward, joined by meeting                       |
| `draft`                                        | `draft::draft`                                                        |
| `write_step`, `read_step`                      | `io::step::write`, `io::step::read`                                   |
| `sew`                                          | `sew::sew` of the faces `io::step::read` reads                        |
| `volume`, `area`, `centroid`, `inertia`        | `properties::volume`, `properties::area`                              |
| `normal`, `curvature`                          | `properties::normals`, `properties::curvatures`                       |
| `interference`                                 | `boolean::interference`                                               |
| `raycast`                                      | `raycast::raycast`, `raycast::raycast_all` for a batch                |
| `bounding_box`                                 | `properties::tight_bounding_box`                                      |
| `check`                                        | `validity::check`, its errors                                         |
| `classify`                                     | `boolean::PointClassifier`                                            |
| `closest`, `distance`                          | `distance::distance`                                                  |
| `section`                                      | `boolean::section` with the unbounded plane                           |
| `split`                                        | `boolean::cells` with the unbounded plane, pieces by side             |
| `project`                                      | `hidden_lines::hidden_lines` in the protocol's view axes              |
| `mesh`                                         | `meshing::triangulate`                                                |
| `display`                                      | `meshing::display`                                                    |
| `transform`, `affine`                          | `modify::transform`, `modify::transform_affine`                       |

What the adapter adds is translation:

- The suite's profiles become kentron's edges, wires and planar faces (`make::face_bounded_on`):
  lines, circular arcs, a conic as the rational quadratic B-spline on its three points, a periodic
  spline as kentron's periodic B-spline of one period. A frame becomes kentron's `Frame` of the
  same axes, and a sweep's path its spine of lines and arcs.
- A point becomes the edge or face it lies on through kentron's distance query
  (`distance::distance` from a vertex at the point) at `10·ε`. A point that names no edge or face is
  refused as invalid; one within `10·ε` of several is refused as unsupported, not settled here.
  A seam, where a face meets itself, names no edge: the protocol lists seams apart from edges.
- A plane for `section` is kentron's unbounded planar face: a face on the plane with no boundary.
- `sew`'s one closed shell becomes the shell of a solid.
- The empty solid is kentron's empty compound, a shape with no face, edge or vertex (the test
  kentron's booleans make). On it, `classify` gives every point outside, `centroid` none and
  `inertia` zeros; any other shape goes to kentron.
- Where a kentron call takes a choice the protocol's definition makes, the adapter names it: a
  whole offset rounds by arcs (`Join::Arc`), the offsets of `offset_faces` and the walls of `shell`
  meet by extending the faces (`Join::Intersection`), and the ruled `loft` joins each point of the
  profile to its image in the scaled copy, its two sections' edges matched as given
  (`Matching::AsGiven`).

Each kentron error is the protocol's refusal of the same kind: `InvalidInput` as invalid, `Failed`
as failed.

### Unsupported

`hello` leaves these out; each has no one kentron call:

| Function                                                                  | Why                                                                                           |
| ------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------- |
| `betti`, `section_area`                                                   | Each only through other calls: counts of shapes, a section's curves                           |

Also refused as unsupported: a blend section (`gn`) segment in a profile, a loft section with holes
(kentron lofts wires), and `thicken` of more than one face (kentron thickens one sheet).

[`Cargo.toml`](Cargo.toml) pins one kentron release from crates.io, 0.21.1, and `VERSION` in
[`src/main.rs`](src/main.rs) reports it; moving to another release is a change made here.

## Setup

The Rust toolchain [`rust-toolchain.toml`](rust-toolchain.toml) pins.

## Usage

The suite runs the adapter as a separate process. Build the adapter here, then run the suite from a
checkout of [brep-conformance](https://github.com/springcraft-ai/brep-conformance) beside this one:

```bash
cargo build --release
cd ../brep-conformance
cargo build --release -p brep-conformance
RAYON_NUM_THREADS=1 target/release/brep-conformance \
    --adapter ../brep-conformance-kentron/target/release/kentron-adapter --out kentron.json
```

kentron runs the independent steps of an operation on a thread pool. The suite runs one worker per
thread, so `RAYON_NUM_THREADS=1` keeps each worker to one; it changes no result.

## Licences

The adapter is under the MIT licence: see [LICENSE](LICENSE). That licence covers the adapter and
nothing else.

kentron is a separate work under its own licence, `AGPL-3.0-or-later AND BSD-3-Clause AND MIT` as
its `Cargo.toml` declares (its SLSQP optimizer is translated from NLopt's). The adapter links kentron
statically, so a built adapter binary is subject to kentron's terms as well as this repository's.
kentron's own repository is the authority on its licence, not this page.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).
