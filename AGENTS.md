# Contributor instructions

## Project scope

BcubeCGNS implements Bcube's IO interface as a basic single-zone CGNS
reader/writer. Keep changes focused on that interface; distributed mesh
processing and multizone assembly are outside the project scope.

## Implementation and validation

Keep code readable in explicit stages, reuse existing helpers, and document
nonobvious choices. Keep tests compact and well commented: explain fixture
layouts, expected results, reference provenance and the failure each test detects.
Reuse small fixtures for related checks. Normal tests must be self-contained
apart from the committed Git LFS assets.

Run scripts sequentially with a single process. Use
`julia --project=<directory-containing-Project.toml> <script.jl>` and Julia's
package manager for dependency changes. Format changed Julia files with
`jlfmt --inplace <path-to-file>` before validation and commit coherent changes
with relevant validation.

Preserve accepted Julia-dependent VTK checksum references. Investigate output
differences before updating hashes; do not replace them merely to make tests pass.
See `test/assets/README.md` for fixture provenance and regeneration recipes,
`docs/reader-coverage.md` for reader support, and
`docs/eclipps-conventions.md` for ECLIPPS conventions.
