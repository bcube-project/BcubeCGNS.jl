# Reader support and limitations

The HDF5 `Bcube.read_file` reader supports one CGNS base and one zone.
Eclipps surface extraction is described separately in
[Eclipps conventions](eclipps-conventions.md). Implementation presence does not
establish complete CGNS compliance or extensive validation of every mapped type.

| Area | Current implementation | Tested evidence | Limits / qualifications |
| --- | --- | --- | --- |
| Base/zone selection | Exactly one base and zone; no `domains` selection | Single-zone assets | No domain selection or multizone assembly |
| Coordinates | Exact Cartesian axis names, consecutive axes, validated lengths and space dimension | Public two-quad assertions; malformed/missing axes | Cartesian axes only; automatic detection retains leading axes |
| Fixed-arity unstructured sections | Mapped Bcube types, ordered by element ID; validated arity, node bounds, ranges and cell count | Line mesh; multisection Eclipps; rectangle; descending-order two-quad tree | Mapping is broader than independently tested geometry; see below |
| Full-zone fields | Validated Vertex/CellCenter arrays; Vertex default; variable filtering; sorted cell association, including boundary-first numbering | Direct field placement/default/length, section-order, MIXED association, lower-dimensional equal-count rejection and selected-step mesh/field consistency | Partial/rind and general face fields excluded; no volume-to-surface projection |
| Unstructured BCs | Vertex lists/ranges with Vertex default; exact exterior face matching for 2D EdgeCenter and FaceCenter/legacy CellCenter element selections; shared-family aggregation | Two-quad tree: descending/noncontiguous selections, three sections, invalid IDs/gaps/dimensions, exact face membership and consistent tags; disjoint family patches | Per-patch Vertex node-set semantics; 2D FaceCenter/CellCenter compatibility aliases |
| Structured zones | Bar2/Quad4/Hexa8 connectivity; Vertex and directional BC lists/ranges | Independent 1D/2D/3D public-reader fixtures; exact boundary membership, fields and signed geometry | Directional 2D aliases are compatibility extensions; structured external global-numbering metadata ignored |
| MIXED sections | All embedded types must map to fixed-arity Bcube entities; modern offsets and legacy streams | Independent quad/tri/edge fixture; exact connectivity, fields and BCs; complete-section rejection | Independent MIXED geometry covers BAR_2/TRI_3/QUAD_4; other mapped codes lack equivalent fixtures |
| Time selection | Explicit one-based `timestep`, iterative grid/solution pointers; bare/qualified zone pointers and validated activation counts; constant unambiguous data | Independent two-step fixture; ordering, missing/malformed metadata, inactive steps, qualified-name limits, mesh/field consistency | Motion transforms and changing topology excluded |
| NGON/NFACE | Modern offsets; convex planar TRI/QUAD loops and oriented TETRA/PYRA/PENTA/HEXA shells | Independent two-cube file; literal tetrahedron, pyramid, prism and 2D cells; malformed representations | Arbitrary polygons/polyhedra, warped faces and legacy count-prefixed streams excluded |
| Eclipps extraction | One selected BC; structured/unstructured surface lists; Vertex/FaceCenter datasets with explicit association, ordering and remapping | Four accepted VTK regressions; literal 2D/3D trees and temporary public-file checks | Bounded compatibility conventions in [Eclipps conventions](eclipps-conventions.md); partial datasets excluded |
| Unsupported representations | Explicit errors for unsupported element codes, solution locations, partial solutions and rind | Focused connectivity and unsupported field-location tests | Validation is bounded by the supported encodings, not complete SIDS conformance |

Mapped element codes are NODE, BAR_2/3, TRI_3/6/9/10, QUAD_4/8/9,
TETRA_4/10, PYRA_5, PENTA_6 and HEXA_8 (`src/common.jl`). Independent fixed/MIXED
geometric evidence covers BAR_2, TRI_3 and QUAD_4;
structured connectivity covers BAR_2/QUAD_4/HEXA_8, and polyhedral conversion
covers TRI_3/QUAD_4/TETRA_4/PYRA_5/PENTA_6/HEXA_8. A mapped higher-order
code alone does not establish independently verified geometry or orientation.

Existing `test/test_cgns_mesh.jl` additionally compares cell counts and integrated
volume between the Bcube mesh and CGNSMesh paths for the rectangle and
multisection NACA assets. Those checks share input and implementation and are
not independent fixture oracles. Existing writer tests remain regressions.

## Field and boundary association

Full-zone Vertex and CellCenter arrays must match the zone's vertex/cell count.
Cells follow increasing original element IDs; structured arrays flatten with i
fastest. Compact cell-field ordering for boundary-first or interspersed element
numbering preserves the existing reader's compatibility convention.

Lower-dimensional unstructured extraction supports mesh-only reads and full-zone
Vertex arrays. Requested CellCenter arrays are rejected because they belong to
the original zone cells; no volume-to-surface projection is performed. Partial
solutions, rind and general face-centered FlowSolutions are unsupported.

Vertex BCs use per-patch node-set semantics. Explicit element selections must
match exterior mesh faces exactly. Patches sharing a family produce one public
boundary name; each patch resolves independently before its faces are unioned.
Reported nodes retain first occurrence in the decoded patches. Overlapping
boundary names retain Bcube's first-boundary precedence. Empty ZoneBC containers
are accepted.

## Stored timesteps

`timestep=k` selects the kth stored step, one-based, independently of physical
time values and child order. `timestep=:latest` selects the last stored step.
Files with iterative metadata require a selection, including mesh-only reads.
Without iterative metadata, `timestep=nothing` reads all present FlowSolutions;
`:latest` does the same with a warning, and integer indices are rejected.
Malformed metadata never uses this fallback.

Grid and solution pointers resolve exact names. Without a pointer, at most one
grid/solution can be treated as constant. Null solution pointers yield no fields;
Null grids cannot yield a mesh. Bare and base-qualified zone names are supported;
activation counts and pointers must agree. Inactive zones, motion transforms,
selected grid-connectivity pointers and changing topology are unsupported.

## NGON/NFACE limits

Modern zero-based ElementStartOffset arrays are required. NGON loops must form
convex planar TRI/QUAD faces. Signed NFACE references must form oriented,
convex, nonzero-volume TETRA/PYRA/PENTA/HEXA shells. Matching uses complete
oriented face cycles. Geometry tolerances scale with local diameter (1e-10).
Arbitrary polygons/polyhedra, concave cells, warped faces, open or inconsistent
shells, nonmanifold reuse and legacy count-prefixed streams are rejected.
MIXED sections reject embedded NGON/NFACE and unsupported codes, including in
unselected records. Ordinary fixed-arity cells do not receive these polyhedral
geometry checks.

## Validation and tests

Malformed or unsupported file content raises contextual errors. Invalid public
argument choices use `ArgumentError`; internal invariants use assertions.
Some direct iterative-metadata checks also use `ArgumentError`.

Tests combine original file regressions with small independent fixtures and
literal trees. See [fixture provenance and regeneration](../test/assets/README.md).
The suite was validated on Julia 1.13.3 with Bcube 0.5.0: 447 assertions pass,
including the accepted VTK checksums and writer regressions. This does not claim
validation of every declared older Bcube version. Normal tests require neither
an external fixture generator nor the Eclipps paper.
