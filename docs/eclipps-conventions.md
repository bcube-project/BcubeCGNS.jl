# CGNS-ECLIPPS surface reader conventions

These conventions refer to Montreuil et al., *Enhancement of Prediction
Capability in Icing Accretion and related Phenomena, Part 1*, AIAA 2009-3969
(2009), pages 4–5, Tables 1–3 and Figures 3–5. The paper is not redistributed
and is not a test dependency.

## What the paper documents

The neutral layer restricts CGNS for interoperability between icing modules.
Page 4 requires dimensional exchange, with SI dimension exponents in the
tables. These are module input requirements, not a requirement that a general
reader reject an incomplete set of quantities.

| Source | Quantities and tree organization |
| --- | --- |
| Table 1 / Figure 3 | `ReferenceState`: `VelocityMagnitudeReference`, `AngleOfAttackReference`, `TemperatureReference`, `PressureReference`, `LengthReference`, `WaterDensityReference`, `LiquidWaterContentGlobal`, `MedianVolumicDiameterGlobal`, `DurationAccretion`, `NumberDropletClass`; per-class `LiquidWaterContentClass`, `MedianVolumicDiameterClass` beneath `DropletClass#01` in the figure |
| Table 2 / Figure 4 | Wall gas `Temperature`, `Pressure` in `IcingWallData:Gas/DirichletData`; `SkinFrictionX/Y/Z`, `CoefConvectiveExchange` in `NeumannData`; droplet-class `DropletVelocityX/Y/Z`, `DropletTemperature` in Dirichlet and `MassFlux` in Neumann |
| Table 3 / Figure 5 | Volume `VelocityX/Y/Z`, `Temperature`, `Pressure` under `Predictor:Gas`, with `GridLocation`; turbulence variables shown but not mandatory |

Figure 4 also shows `ZoneBC/Wall#01`, `PointList`, `GridLocation` and
`SurfacicElementList`. It does **not** define the latter's payload, element
numbering, corner order, or dataset remapping algorithm. Nor does it establish
that every dataset is node-centered. The reader preserves values without
renaming, dimensional conversion, beta/mass-flux conversion, or enforcement of
these solver-required quantity sets. The printed MassFlux dimension row is not
used to infer or alter units. Volume fields remain the separate full-zone
Vertex/CellCenter reader's responsibility; general face fields are excluded.

## CGNS rules and implementation compatibility conventions

Dataset location inheritance, optional dataset selections, and scalar versus
local BC arrays follow [CGNS SIDS boundary conditions, sections 9.3–9.5](https://cgns.org/standard/SIDS/boundary.html).
The following bounded extraction rules supplement that standard; they are not
claimed as a complete ECLIPPS specification.

| Area | Implemented contract / provenance |
| --- | --- |
| Selection | One BC node name, independently of ZoneBC child order; string and singleton collection APIs preserved. Exactly one base/zone. |
| Unstructured surface list | `SurfacicElementList/List` contains original zone element IDs in desired output cell order. Observed in accepted assets; uses validated fixed-arity/MIXED/representable polyhedral decoding. Optional global-numbering metadata does not replace zone element IDs. |
| Structured surface list | Flattened IndexDimension-tuples for each Bar2 (2D) or Quad4 (3D) corner, in element order. Observed 3D encoding, extended consistently to 2D, including embedding in 3D. SpaceDimension is not IndexDimension. Stored orientation is preserved. |
| Vertex BCs | Missing location defaults to Vertex. PointList order and inclusive PointRange order (i fastest, descending supported) determine output nodes. Unused nodes remain. Repeated entries in the accepted multisection asset remain separate nodes/data slots; connectivity uses the last occurrence as historically. All surface vertices must be present. |
| FaceCenter BCs | Cells follow surface-list order, values follow BC selection order and are remapped. Unstructured nodes are sorted; structured nodes follow first occurrence, preserving historical compaction behavior. |
| Structured face indices | Minimum connecting vertex tuple identifies a face already defined by its corner list. This is an explicit compatibility convention, not specified by the paper. Ambiguous keys are rejected when face data is requested; directional/general face datasets remain excluded. |
| Datasets | Omitted location inherits BC location; explicit Vertex/FaceCenter supported. A changed location requires an explicit selection. Explicit selections must cover the entire extracted node/cell set and are reordered to mesh order. Partial/foreign/ambiguous selections fail. |
| Values | Dirichlet-only, Neumann-only, empty groups/datasets and arbitrary named numeric quantities are accepted. One-entry arrays broadcast; other arrays must match selection length. Location determines MeshPointData/MeshCellData even when node and cell counts coincide. |
| Errors | Missing/malformed lists, invalid indices/lengths, duplicate surface IDs, unsupported locations and variable-name collisions fail clearly. Unsupported embedded MIXED records anywhere still abort, even if unreferenced by the BC. |

## Public API decisions

`extract_surf_from_eclipps(path, name; verbose=false, timestep=nothing)` retains
`(; mesh, data, refState)` and `data[dataset][variable]` wrappers. Reference
scalars retain the existing scalar return representation. Artificial contiguous
mesh global IDs and valid existing fixture output are preserved. `timestep`
is an additive keyword with the general reader's explicit one-based stored-step
semantics; iterative grids cannot be read accidentally without selection.

The flat data dictionary cannot encode both a Dirichlet and a Neumann variable
with the same name. Such a collision now raises an error instead of silently
overwriting one. Representing both would require a public return-shape change;
none is introduced. Partial datasets would likewise need an
association/mask contract and remain explicitly unsupported. No solver-input
completeness option or multizone assembly is added.

## Independent verification

`test/test_reader_eclipps.jl` contains literal synthetic 2D/3D trees with
hand-checked coordinates, connectivity, signed orientation, centers and data
association. It reuses small geometries for ordering, location, range,
remapping and malformed-input checks. A temporary HDF5 file uses explicit
node/attribute serialization to exercise the public API; it does not use the
package mesh writer as an oracle. The four accepted real-file VTK regressions
remain compatibility checks. See the [fixture provenance](../test/assets/README.md).

## Structured surface topology

Structured corner lists must describe elementary exterior edges/faces in logical
index space: one exterior fixed axis, unit tangential spans, distinct complete
corners and a cyclic quad order. Diagonals, skipped vertices, interior surfaces
and crossed/non-elementary quads fail before mesh construction. Both valid
orientations and supplied element order are preserved. IndexDimension remains
independent of coordinate embedding; repeated/unused Vertex selections and
full-selection dataset permutations retain their compatibility behavior.
