# Reader fixtures

The six original CGNS assets are retained baseline regressions. Their Git LFS
pointers record their content hashes; this work does not establish new
redistribution provenance for those pre-existing assets.

## two-quads.cgns

New synthetic fixture generated locally using ONERA Maia 1.10, revision
`54560ee1ba0f61a65f698ca07f8b5d25d6050f95`
([upstream](https://github.com/onera/Maia/tree/54560ee1ba0f61a65f698ca07f8b5d25d6050f95)).
Original mesh/data values were defined for this project and are distributed
under the repository MIT license; no upstream mesh or solution data was copied.
Maia is optional for regeneration and is not a test dependency. No Maia
Python/YAML sources are committed.

- Size: 16,744 bytes.
- SHA-256: `4a4cced20d00e3642b6083c709f9e703cd259faa69b25e8df84546adfbb3ba73`.
- One base with cell/physical dimensions `(2,2)`, one unstructured zone
  with six vertices and two cells, no boundary vertices in the zone size.
- Cartesian coordinates, in node order:
  `(0,0), (1,0), (3,0), (0,2), (1,2), (3,2)`.
- `Quads`: QUAD_4, ElementRange `[1,2]`, flattened connectivity
  `[1,2,5,4, 2,3,6,5]` (counterclockwise).
- `Left`: BCWall, Vertex PointList `[1,4]`, the single left exterior edge.
- `Nodes`: Vertex FlowSolution, `node_id = [1,2,3,4,5,6]` (Float64).
- `Cells`: CellCenter FlowSolution, `cell_id = [10,20]` (Float64).

To regenerate, construct this tree using `maia.pytree.new_*` constructors and
`maia.io.write_tree` in a single Python process with the Maia revision linked
above. Use Int32 zone size, element range, connectivity and point list arrays;
use Float64 coordinates and fields.
Maia/HDF5 versions may change serialization and thus the hash: verify the tree
and semantic expectations before updating the recorded digest.

Expected topology and data in `test_reader_semantics.jl` are hand-computed,
independently of the BcubeCGNS writer. Cell centers are `(0.5,1)` and `(2,1)`;
signed areas are `2` and `4`. The unequal widths and distinct field markers
make cell permutations and field association errors visible.

## structured-1d.cgns, structured-2d.cgns, structured-3d.cgns

Synthetic structured fixtures created with the same Maia 1.10 revision and
MIT provenance as `two-quads.cgns`; no upstream data was copied. Normal tests
require only the committed CGNS files, not Maia.

| File | Bytes | SHA-256 |
| --- | ---: | --- |
| structured-1d.cgns | 13,649 | `d1f90fe381809b687abb99c194da4439746b4fe56ecb627f6093f149eb379b78` |
| structured-2d.cgns | 14,388 | `4aea3774a3831ebd8f9b4e2934dc70e248abe84e1cc21cdc253267bbc85d4aed` |
| structured-3d.cgns | 15,149 | `cd426b23346f8fbe5f6fa3657ee5fd190882f22bf12c9375ba9973b874c94000` |

Each file contains `Base/Grid`, with matching cell/physical dimensions 1, 2
or 3. Vertex axes are `x=(0,1,3)`, `y=(0,2,5)`, `z=(0,4)`, truncated to the
file dimension. Zone size rows are `(n,n-1,0)`; coordinate arrays have the
full vertex shape, in Fortran order. `Left` is a Vertex BCWall PointRange
from all-ones to `(1,ny,nz)`, truncated to the dimension. `Nodes/id` contains
one-based linear vertex IDs; `Cells/id` contains ten times the one-based
linear cell IDs. Both fields have their location's multidimensional shape.

To regenerate serially, use Maia `new_CGNSTree`, `new_CGNSBase`, `new_Zone`
with `type='Structured'`, `new_GridCoordinates`, `new_DataArray`, `new_ZoneBC`,
`new_BC` and `new_FlowSolution`, then `maia.io.write_tree`. Use Int32 sizes
and point ranges, Float64 coordinates/fields, and Fortran-contiguous arrays.
No element sections are required. Verify semantic expectations before accepting
serialization/hash differences. Test variants reuse these trees for directional
BC lists/ranges and malformed selections; no package writer supplies expectations.

The hand-checked signed lengths are `(1,2)`, areas `(2,4,3,6)` and volumes
`(8,16,12,24)`, in i-fastest cell order. Expected corner connectivity and
boundary face node sets are written explicitly in `test_reader_structured.jl`.

## mixed-quads-tris.cgns

Synthetic fixture created with the same Maia 1.10 revision and MIT provenance
as `two-quads.cgns`. It reuses those six coordinates and node field markers;
no upstream mesh/data was copied. Maia is not required by normal tests.

- Size: 20,545 bytes.
- SHA-256: `0bd8f0a9ebbb6806e2c3c2be0994adf6ad8282d0ed496dbc23385991556f83ab`.
- `Base` dimensions `(2,2)`; unstructured `MixedGrid` size `(6,3,0)`.
- The left quad is unchanged; the right quad is split along edge `2->6`.
- `Later`, MIXED, IDs `3:5`, appears first in child order. Its encoded
  connectivity is `[5,2,6,5, 3,1,4, 3,3,6]`, with offsets `[0,4,7,10]`.
- `First`, MIXED, IDs `1:2`, appears second. Encoded connectivity is
  `[7,1,2,5,4, 5,2,3,6]`, with offsets `[0,5,9]`.
- Embedded codes `7`, `5`, `3` mean QUAD_4, TRI_3, BAR_2.
- `Ends`: BCWall, FaceCenter PointList `[5,4]`, selecting the right and left
  edges, not the triangle sharing their section. This retains the historical
  2D FaceCenter compatibility convention documented in the reader matrix.
- `Cells/cell_id = [10,20,30]`; `Nodes/node_id = [1,2,3,4,5,6]`.

To regenerate serially, read the independent two-quad tree with Maia, replace
its element sections by the above `new_Elements(type='MIXED')` records, add
ElementStartOffset DataArrays, set the zone cell count to three, replace the
BC and cell field, and use `maia.io.write_tree`. Sizes, ranges, connectivity,
offsets and lists are Int32; coordinates/fields are Float64. Verify a Maia
readback is the same tree and recheck expectations before updating hashes.

Hand-computed cell centers in CGNS ID order are `(1/2,1)`, `(7/3,2/3)` and
`(5/3,4/3)`; signed areas are `(2,2,2)`. Original IDs and cell field markers
expose section-order errors. The BC edges have node sets `{1,4}` and `{3,6}`.

## timesteps.cgns

Synthetic two-step fixture created serially with Maia 1.10, using the same
revision and MIT provenance as `two-quads.cgns`; no upstream data was copied.
Size: 24,996 bytes.
SHA-256: `d3550d49592a1cc4540d0099e3c740e98de0927699f7744a9454115b6a2acb8c`.
Normal tests require only the committed binary, not Maia.

To regenerate, read `two-quads.cgns` with `maia.io.read_tree`, retain its
topology and BCs, and replace its grids/solutions with two copies:
`GridAlpha` has the original coordinates; `GridZulu` adds 10 to CoordinateX.
`SolutionAlpha` is the original CellCenter `Cells` solution (`cell_id=[10,20]`);
`SolutionZulu` adds 100 to that field. Remove the original Nodes solution.
Add SimulationType `TimeAccurate`, BaseIterativeData `StoredSteps` with Int32
NumberOfSteps 2, TimeValues `[9,2]` and IterationValues `[90,20]`.
ZoneIterativeData `StepMap` contains standard space-padded C1 `(32,2)`
GridCoordinatesPointers `[GridZulu,GridAlpha]` and FlowSolutionPointers
`[SolutionZulu,SolutionAlpha]`. Write with `maia.io.write_tree`.

The grids and solutions are inserted Alpha before Zulu, opposite pointer order.
Hand-computed expectations select x+10 and `[110,120]` at stored step 1,
original coordinates and `[10,20]` at step 2. Nonmonotonic times prevent
accidental time sorting. Connectivity, original IDs and BCs stay constant.
Tree variants exercise invalid metadata and complete MIXED validation.

## poly-cubes.cgns

Synthetic fixture created serially with Maia 1.10, revision
`54560ee1ba0f61a65f698ca07f8b5d25d6050f95`. Original synthetic data,
redistributed under this repository's MIT license; no upstream mesh copied.
Size: 19,395 bytes. SHA-256:
`bb089f72d5278bb7928ebee3784cd51fbec3c193d1db3e5bd5326deece737a0b`.
Normal tests use the committed binary without Maia.

Regenerate with `maia.pytree` constructors and `maia.io.write_tree` in one
serial process. Base dimensions are `(3,3)`, unstructured zone size `(12,2,0)`.
Vertices are `(x,y,z)` in nested loops z=0:1, y=0:1, x=0:2 (x fastest).
`Faces` is NGON_n IDs 1:11 with loops:
`[1,4,5,2], [1,2,8,7], [2,5,11,8], [5,4,10,11], [1,7,10,4],
[7,8,11,10], [2,5,6,3], [2,3,9,8], [3,6,12,9], [6,5,11,12], [8,9,12,11]`.
Its ElementStartOffset is `0:4:44`. `Cells` is NFACE_n IDs 12:13 with
records `[1,2,3,4,5,6]` and `[7,8,9,10,-3,11]`, offsets `[0,6,12]`.
All connectivity, ranges and offsets are Int32; coordinates/fields Float64.
`CellValues/marker` is CellCenter `[12,13]`; `Ends` is a FaceCenter BC of type Null with
PointList `[5,9]`. Use ElementStartOffset DataArrays under Elements_t.

The two cubes have centers `(0.5,0.5,0.5)`, `(1.5,0.5,0.5)` and signed
volumes +1. Face 3 is shared with opposite signs. End boundary vertex sets
are `{1,4,7,10}` and `{3,6,9,12}`. Tests also construct literal tetrahedron,
pyramid, prism and planar triangle/quad trees with hand-checked outward loops;
these do not use the package writer as an oracle. Malformed variants reuse
the same data to exercise offsets, references, orientation and geometry.

Serial Maia readback verified zone size, X coordinates, signed NFACE references
and cell markers against the literal expectations above.

## Literal Eclipps fixtures

`test/test_reader_eclipps.jl` supplies original synthetic trees under this
repository's MIT license; no external mesh or solver output was copied.
They require neither Maia nor committed binary fixtures. The source literals
are the regeneration recipe and provenance; binary hashes are inapplicable to
the temporary HDF5 files created and deleted by the tests.

The 2D fixture uses the two-quad coordinates `(0,1,3) × (0,2)`, embedded in
3D with zero Z. Original quad IDs are 1–2; bottom Bar2 IDs 3–4 have
connectivity `[1,2]`, `[2,3]`. Surface order `[4,3]` reverses those IDs.
Vertex order `[3,2,1,6]` includes unused node 6; nodal temperatures are
`[30,20,10,60]`. Face BC order `[3,4]` carries temperatures `[300,400]`,
so extracted cells must receive `[400,300]`. Structured equivalents use
IndexDimension 2 corner tuples, preserving the same geometry and ordering.
Ranges, explicit dataset selections and malformed variants reuse these trees.

The 3D fixture is one `1 × 2 × 3` hexahedron, vertices i fastest on
`(0,1) × (0,2) × (0,3)`. The z=0 quad follows `[1,2,4,3]` with signed
XY area +2 and center `(0.5,1,0)`, deliberately preserving its +z stored
orientation. Structured and unstructured versions select the same reversed
vertex order `[3,4,2,1]`. A separate triangle `(1,2,4)` and its three CCW
boundary edges have equal surface cell/node counts to detect incorrect inference
of data location from array length. A two-grid
variant uses reversed pointers and decreasing times to verify stored-step
selection during surface extraction.

All expected coordinates, connectivity, signed measures and data permutations
are hand-computed. Temporary public-API fixtures use explicit CGNS/HDF5
attributes and literal payloads, independently of the package mesh writer.
Existing NACA/RG Eclipps assets remain broad compatibility regressions.
The Montreuil et al. paper is a convention reference only and is not a test
dependency. See the
[paper audit and convention limits](../../docs/eclipps-conventions.md).

## Adding or regenerating fixtures

Record the original data author/license, upstream URL and exact tool revision,
byte size and SHA-256 (`sha256sum <file>`) for each committed binary. Check
redistribution terms before downloading public meshes; the six pre-existing
assets above remain accepted regressions with provenance limitations. Use the
individual sizes/digests above as the verification authority. Keep the recipe
in prose or Julia literals; do not
commit Maia Python/YAML sources, environment artifacts or the supplied PDF.

Regenerate one file at a time using the Maia revision linked above.
Read it back independently and compare topology, coordinates, ordering and
fields against the recorded literal recipe and hand calculations before
accepting a new serialization digest. A byte-identical file is optional across
HDF5 versions; equivalent semantics are required. Add direct public-reader
assertions for the new behavior, reuse the small geometry for malformed variants,
and explain expected numbering/orientation and failure modes in the tests.
The package writer may supply regression output, but must not be the sole oracle.
