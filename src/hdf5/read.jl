"""
    Bcube.read_file(
        ::CGNSHDF5IoHandler,
        filepath::String;
        domains = String[],
        varnames = nothing,
        topodim = 0,
        spacedim = 0,
        verbose = false,
        timestep = nothing,
        kwargs...,
    )

Read a single-zone CGNS/HDF5 file as `(; mesh, data, refState)`.

`timestep=k` selects the kth stored step (one-based), through iterative grid
and solution pointers, independently of time values or child order.
`:latest` selects the last stored step. Iterative files require selection,
including mesh-only reads. Without iterative metadata, `nothing` reads all
present FlowSolutions; `:latest` does the same with a warning, and integer
indices are rejected. Invalid metadata always raises an error.

For `varnames`, `domains`, `spacedim` and `verbose`, see `Bcube.read_file`.
Domain selection is not supported yet. `topodim=0` uses base topology; positive
overrides select stored elements of that dimension. Extra keywords are accepted
and currently ignored.

Only Vertex/CellCenter FlowSolutions are supported; grid-motion transforms,
multizone assembly and changing topology are excluded.

Cell order follows increasing original element
IDs (structured i fastest), not section order. Invalid encodings and field
associations raise errors as described by the decoding helpers.

# Return
See `Bcube.read_file` for more details. `refState` is a nested dictionary or `nothing`.
"""
function Bcube.read_file(
    ::CGNSHDF5IoHandler,
    filepath::String;
    domains = String[],
    varnames = nothing,
    topodim = 0,
    spacedim = 0,
    verbose = false,
    timestep = nothing,
    kwargs...,
)
    return h5open(filepath, "r") do file
        tree = CGNS.hdf5_to_node(file)
        return read_tree(
            tree;
            domains,
            varnames,
            topodim,
            spacedim,
            verbose,
            timestep,
            kwargs...,
        )
    end
end

"""
    read_tree(
        tree;
        domains = String[],
        varnames = nothing,
        topodim = 0,
        spacedim = 0,
        verbose = false,
        timestep = nothing,
        kwargs...,
    )

Read one CGNS base/zone into `(; mesh, data, refState)`.

Arguments, timestep selection and return values follow `Bcube.read_file` for
CGNS/HDF5. The input tree is not mutated. Multiple/missing bases or zones and nonempty domain
selection are errors. Full-dimensional reads validate the zone cell count;
lower-dimensional unstructured reads select stored boundary elements instead.
Requested full-zone CellCenter arrays cannot be associated with that selection;
mesh-only and Vertex-only reads remain valid.
"""
function read_tree(
    tree;
    domains = String[],
    varnames = nothing,
    topodim = 0,
    spacedim = 0,
    verbose = false,
    timestep = nothing,
    kwargs...,
)
    isempty(domains) || throw(ArgumentError("Reading selected domains is not supported"))

    # Read (unique) CGNS base
    cgnsBase = get_cgns_base(tree)

    # Reference state
    refState = read_ref_state(cgnsBase)

    # Read base dimensions (topo and space)
    dims = get_value(cgnsBase)
    topodim = topodim > 0 ? topodim : dims[1]
    zone_space_dim = dims[2]
    verbose && println("topodim = $topodim, spacedim = $(zone_space_dim)")

    # Find the list of Zone_t and select the one matching the `timestep` argument
    zones = get_children(cgnsBase; label = CGNS.Zone_t)
    if length(zones) == 0
        error("Could not find any Zone_t node in the file")
    elseif length(zones) > 1
        error("The file contains several Zone_t nodes, only one zone is supported for now")
    end
    zone = select_timestep(cgnsBase, first(zones), timestep)

    # Read zone
    zoneCGNS = read_zone(
        zone,
        varnames,
        topodim,
        spacedim,
        verbose;
        validate_cell_count = topodim == dims[1],
        zone_dim = dims[1],
    )

    # Build Bcube Mesh
    mesh = cgns_mesh_to_bcube_mesh(zoneCGNS)

    # Build Bcube MeshCellData & MeshPointData
    if !isnothing(varnames)
        data = flow_solutions_to_bcube_data(zoneCGNS.fSols)
    else
        data = nothing
    end

    return (; mesh, data, refState)
end

"""
    iterative_pointer(iterative, name, nsteps, step)

Return the selected, whitespace/NUL-trimmed name from a standard iterative
pointer array, or `nothing` if the array is absent.

The HDF5 tree decoder represents C1 `(32, nsteps)` arrays as a vector of strings. Reject malformed
lengths, empty names and names exceeding CGNS's 32-byte limit. `"Null"` is
returned literally so the caller can distinguish absent data from no pointer.
`iterative` may be `nothing`; `name` is the exact DataArray name, `nsteps` the
positive stored-step count and `step` an already validated index in `1:nsteps`.
Duplicate arrays also raise errors; bounds are the caller's responsibility.
"""
function iterative_pointer(iterative, name, nsteps, step)
    isnothing(iterative) && return nothing

    # Get the iterative DataArray
    dataArray = get_child(iterative; name, label = CGNS.DataArray_t, unique = true)
    isnothing(dataArray) && return nothing

    # Read the values
    values = get_value(dataArray)

    # Values should be a vector of strings and of length NumberOfSteps (nsteps)
    (values isa AbstractVector{<:AbstractString} && length(values) == nsteps) ||
        error("$name must contain $nsteps stored-step names")

    # Read the names and select the `step`-th one.
    names = [strip(rstrip(v, '\0')) for v in values]
    all(n -> !isempty(n) && ncodeunits(n) <= 32, names) ||
        error("$name contains an empty or invalid CGNS name")

    return names[step]
end

"""
    timestep_nodes(zone, iterative, label, pointer, nsteps, step)

Resolve a grid or solution pointer by exact name and label.

Without pointers,
a single node is constant across steps; multiple candidates are ambiguous and
rejected. `Null` selects no solution; a missing grid is rejected by coordinate
validation. Return a node vector, preserving names and array ordering.
"""
function timestep_nodes(zone, iterative, label, pointer, nsteps, step)
    nodes = get_children(zone; label)
    name = iterative_pointer(iterative, pointer, nsteps, step)
    if isnothing(name)
        length(nodes) <= 1 || error("Missing $pointer for multiple $label nodes")
        return nodes
    end
    # Null explicitly selects no data; preserve the node vector's element type.
    name == "Null" && return nodes[1:0]
    selected = filter(n -> get_name(n) == name, nodes)
    length(selected) == 1 || error("$pointer references missing or ambiguous '$name'")
    return selected
end

"""
    validate_zone_activation(base, zone, iterative, nsteps, step)

Validate optional single-zone activation metadata and return nothing.

NumberOfZones must be a vector of length `nsteps`, with counts zero or one.
ZonePointers is a decoded string matrix `(MaxNumberOfZones, nsteps)`; with
counts present its row count equals their maximum and each column's non-Null
entries equal that step's count. Without counts, one pointer slot is supported.
Names are bare zone names or `Base/Zone` (32 bytes per token, 65 total), as in
CGNS SIDS 11.1.1. Reject malformed paths, inconsistent counts, unsupported zone
counts, inactive selected steps and selected targets other than this base/zone.
Without pointers, the existing single zone is constant unless counts mark the
selected step inactive. Inputs are not mutated; no zone assembly is performed.
"""
function validate_zone_activation(base, zone, iterative, nsteps, step)
    # Read number of zones
    node_NumberOfZones = get_child(
        iterative;
        name = "NumberOfZones",
        label = CGNS.DataArray_t,
        unique = true,
    )
    numberOfZones = isnothing(node_NumberOfZones) ? nothing : get_value(node_NumberOfZones)

    # Check the size of the array `numberOfZones`, it must match `nsteps`
    if !isnothing(node_NumberOfZones)
        (
            numberOfZones isa AbstractVector{<:Integer} &&
            !(eltype(numberOfZones) <: Bool) &&
            length(numberOfZones) == nsteps &&
            all(n -> n in (0, 1), numberOfZones)
        ) || error("NumberOfZones must contain $nsteps integer counts (zero or one)")
    end

    # Read ZonePointers, and make should only one Zone is targeted for step `step`
    node_ZonePointers =
        get_child(iterative; name = "ZonePointers", label = CGNS.DataArray_t, unique = true)
    if isnothing(node_ZonePointers)
        (isnothing(numberOfZones) || numberOfZones[step] == 1) ||
            error("NumberOfZones: selected timestep has no active zone")
        return nothing
    end
    zonePointers = get_value(node_ZonePointers)

    # Several checks about the ZonePointers node
    slots = isnothing(numberOfZones) ? 1 : maximum(numberOfZones)
    (
        zonePointers isa AbstractMatrix{<:AbstractString} &&
        size(zonePointers) == (slots, nsteps)
    ) || error("ZonePointers dimensions must match NumberOfZones and NumberOfSteps")

    # All zone names should smaller than 32 chars.
    names = map(v -> strip(rstrip(v, '\0')), zonePointers)
    for name in names
        tokens = split(name, '/'; keepempty = true)
        (length(tokens) in (1, 2) && all(t -> 1 <= ncodeunits(t) <= 32, tokens)) ||
            error("ZonePointers contains malformed zone path '$name'")
    end

    if !isnothing(numberOfZones)
        all(j -> sum(!=("Null"), view(names, :, j)) == numberOfZones[j], 1:nsteps) ||
            error("NumberOfZones contradicts active ZonePointers")
    end
    selected = filter(!=("Null"), names[:, step])
    length(selected) == 1 || error("ZonePointers: selected timestep has no active zone")
    name = only(selected)
    name in (get_name(zone), get_name(base) * "/" * get_name(zone)) ||
        error("ZonePointers: selected timestep does not reference the supported base/zone")
end

"""
    select_timestep(base, zone, timestep)

Return a non-mutating zone view containing the selected grid and solution.

BaseIterativeData's positive integer NumberOfSteps defines one-based indices;
TimeValues/IterationValues, when present, must match its length. TimeAccurate
SimulationType or ZoneIterativeData without base metadata is invalid. Iterative
files require an explicit integer index (Bool is rejected) or `:latest`.
Without iterative metadata, integer selection is rejected; `:latest` warns and
retains all solutions. Missing pointers permit only unambiguous constant data.
ZonePointers may use bare or base-qualified names and must select this single
zone. Optional NumberOfZones is validated against all stored-step pointers;
see `validate_zone_activation`.

Reject active motion/connectivity pointers, which cannot be applied by this
reader. All Elements_t and BC nodes remain present for complete validation.
"""
function select_timestep(base, zone, timestep)
    # Fetch BaseIterativeData and ZoneIterativeData
    all_BaseIterativeData = get_children(base; label = CGNS.BaseIterativeData_t)
    all_ZoneIterativeData = get_children(zone; label = CGNS.ZoneIterativeData_t)

    # At most one BaseIterativeData and/or ZoneIterativeData can exist
    (length(all_BaseIterativeData) <= 1 && length(all_ZoneIterativeData) <= 1) ||
        throw(ArgumentError("Expected at most one base/zone iterative metadata node"))

    # Check if a "timed" solution is expected.
    simulation = get_child(base; label = CGNS.SimulationType_t)
    timed = !isnothing(simulation) && get_value(simulation) == "TimeAccurate"

    # If no BaseIterativeData is present
    if isempty(all_BaseIterativeData)
        # No BaseIterativeData and no ZoneIterativeData but timed : error
        (isempty(all_ZoneIterativeData) && !timed) ||
            throw(ArgumentError("Missing BaseIterativeData_t"))

        if timestep === :latest
            @warn "No iterative metadata found; timestep=:latest reads all present FlowSolutions."
        elseif !isnothing(timestep)
            throw(ArgumentError("timestep requires stored-step metadata"))
        end

        return zone
    end

    # Read BaseIterativeData
    value = get_value(only(all_BaseIterativeData))
    numberOfSteps = if value isa Number
        value
    elseif value isa AbstractArray && length(value) == 1
        only(value)
    else
        nothing
    end
    (numberOfSteps isa Integer && !(numberOfSteps isa Bool) && numberOfSteps > 0) ||
        throw(ArgumentError("BaseIterativeData NumberOfSteps must be a positive integer"))

    # Safe check : number of entries of "TimeValues" or "IterationValues", if exists,
    # must match "NumberOfSteps"
    for name in ("TimeValues", "IterationValues")
        node = get_child(
            only(all_BaseIterativeData);
            name,
            label = CGNS.DataArray_t,
            unique = true,
        )
        (isnothing(node) || length(get_value(node)) == numberOfSteps) ||
            throw(ArgumentError("$name length must equal NumberOfSteps"))
    end

    # If we reach this point, we should have a value for `timestep`
    isnothing(timestep) &&
        throw(ArgumentError("Time-dependent files require explicit timestep=k (one-based)"))

    # Convert `:latest` into an integer
    timestep = timestep === :latest ? numberOfSteps : timestep
    (timestep isa Integer && !(timestep isa Bool) && 1 <= timestep <= numberOfSteps) ||
        throw(ArgumentError("timestep must be an integer in 1:$numberOfSteps"))

    # Check that the selected zone and timestep is correct
    validate_zone_activation(
        base,
        zone,
        only(all_BaseIterativeData),
        numberOfSteps,
        timestep,
    )

    # Select the ZoneIterativeData node
    zoneIterativeData =
        isempty(all_ZoneIterativeData) ? nothing : only(all_ZoneIterativeData)

    # Raise an error if some specific (unsupported) nodes are detected
    for name in (
        "RigidGridMotionPointers",
        "ArbitraryGridMotionPointers",
        "ZoneGridConnectivityPointers",
    )
        pointer = iterative_pointer(zoneIterativeData, name, numberOfSteps, timestep)
        (isnothing(pointer) || pointer == "Null") ||
            throw(ArgumentError("Selected $name is not supported"))
    end

    # Select the GridCoordinate and FlowSolution(s) nodes corresponding to the timestep
    grids = timestep_nodes(
        zone,
        zoneIterativeData,
        CGNS.GridCoordinates_t,
        "GridCoordinatesPointers",
        numberOfSteps,
        timestep,
    )
    solutions = timestep_nodes(
        zone,
        zoneIterativeData,
        CGNS.FlowSolution_t,
        "FlowSolutionPointers",
        numberOfSteps,
        timestep,
    )

    # Copy of the zone content without any GridCoordinate and FlowSolution(s) nodes, and then "recreate" the zone
    # with the timestep-selected GridCoordinate and FlowSolution(s).
    # In other words, we eliminate the flow solutions that don't correspond to the selected timestep.
    children = filter(
        n -> !(get_label(n) in (CGNS.GridCoordinates_t, CGNS.FlowSolution_t)),
        get_children(zone),
    )
    return CGNS.CGNSNode(
        get_name(zone),
        get_value(zone),
        vcat(children, grids, solutions),
        "Zone_t",
    )
end

"""
    read_zone(zone, args...; validate_cell_count = false, zone_dim = nothing)

Dispatch one Zone_t to the structured or unstructured decoder.

Positional arguments follow `Bcube.read_file`; `usr_space_dim` is its `spacedim`
argument.
`validate_cell_count` checks full-dimensional unstructured cell counts and
explicit BC dimensions. `zone_dim` is the original base cell dimension,
independent of selected `topo_dim`, for field and BC-location validation.
Return `(; coords, c2t, c2n, bcs, fSols, vtx_ln_to_gn, cell_ln_to_gn)`:
coordinates are `(nvertices, spacedim)`, types have one CGNS code per cell,
and connectivity concatenates one-based vertex IDs in cell/corner order.
Absolute IDs follow local vertex/cell order. BCs are records or `nothing`;
solutions are a name-keyed dictionary or `nothing`. See the specialized
readers for supported encodings and validation. Invalid ZoneType asserts.
"""
function read_zone(zone, args...; validate_cell_count = false, zone_dim = nothing)
    # Zone type
    zoneType = read_zone_type(zone)

    if zoneType == "Unstructured"
        read_zone_unstr(zone, args...; validate_cell_count, zone_dim)
    elseif zoneType == "Structured"
        read_zone_struct(zone, args...)
    else
        error("Unhandled ZoneType ($zoneType)")
    end
end

"""
    read_zone_unstr(
        zone,
        varnames,
        topo_dim,
        usr_space_dim,
        verbose;
        validate_cell_count = false,
        zone_dim = nothing,
    )

Read a single unstructured zone with fixed-arity, MIXED or representable
NGON/NFACE sections.

Polyhedral conversion requires simple convex planar
triangles/quads and closed outward-oriented tetrahedron/pyramid/prism/hexahedron
shells; arbitrary polygons/polyhedra and warped faces are rejected.

Validate all section records and vertex IDs before selecting entities of `topo_dim`,
including records outside that dimension.

Unsupported embedded types always raise errors; no partial mesh is returned.

Cells and their original IDs are sorted by CGNS element ID, independently of section child order.

Full-zone CellCenter arrays follow the selected cells in increasing element
ID order, preserving historical support for boundary-first section numbering.

For full-dimensional reads (`validate_cell_count=true`), explicit BC dimensions
are checked before filtering, and the selected count
must equal ZoneCellCount. Lower-dimensional surface extraction deliberately
selects boundary elements instead.

Requested CellCenter arrays are rejected
when `zone_dim` exceeds `topo_dim`; unrequested arrays and Vertex fields remain valid.

Return coordinates `(nvertices, spacedim)`, flattened connectivity, CGNS type
codes, BC records, optional solutions, and absolute vertex/cell IDs.

For `varnames` and `verbose`, see `Bcube.read_file`; `usr_space_dim` is its
`spacedim` argument. The NamedTuple fields are specified by `read_zone`.
Absent Maia metadata, vertex IDs are `1:nvertices` and cell IDs are the selected
original element IDs. `:CGNS#GlobalNumbering/Vertex` and `/Cell`, when present,
are used as supplied, without validating their lengths or permuting them.

Overlapping section ranges, out-of-bounds vertices and missing selected cells
raise errors. Geometric validity of ordinary fixed-arity cells is not
checked by the NGON/NFACE geometry validator.
"""
function read_zone_unstr(
    zone,
    varnames,
    topo_dim,
    usr_space_dim,
    verbose;
    validate_cell_count = false,
    zone_dim = nothing,
)

    # Number of elements
    nvertices, ncells, nbnd = get_value(zone)
    verbose && println("nvertices = $nvertices, ncells = $ncells")

    # Read GridCoordinates
    coords = read_grid_coordinates(zone, topo_dim, usr_space_dim, verbose)

    # Read all elements
    elts = map(read_connectivity, get_children(zone; label = CGNS.Elements_t))
    elts = convert_poly_sections(elts, coords)

    # CGNS cell fields use element ID order, independently of tree child order.
    sort!(elts; by = elt -> first(elt.erange))

    # Check the connectivity
    for (i, elt) in enumerate(elts)
        all(v -> 1 <= v <= nvertices, elt.c2n) ||
            error("Section '$(elt.name)': vertex ID outside 1:$nvertices")
        if i > 1 && first(elt.erange) <= last(elts[i - 1].erange)
            error("Section '$(elt.name)': overlapping ElementRange")
        end
    end

    # Select cells matching the desired topo_dim. A MIXED section can contain both
    # cells and boundary entities. Every record has already been validated.
    c2t, c2n, cell_ids = Int[], Int[], Int[]
    for elt in elts
        filter_volumic_cells!(c2t, c2n, cell_ids, elt, topo_dim)
    end
    isempty(c2t) && error("No elements of topological dimension $topo_dim")
    # Preserve boundary-first numbering used by accepted files: full-zone cell
    # fields follow the selected cells in increasing element ID order.
    (!validate_cell_count || length(cell_ids) == ncells) ||
        error("Volume element count does not match Zone cell count $ncells")

    # Read optional global numbering, retaining original selected IDs by default.
    vtx_ln_to_gn, cell_ln_to_gn =
        read_maia_unstr_numbering(zone, nvertices, cell_ids, verbose)

    # Read all BCs and then keep only the ones whose topo dim is equal to the base topo dim minus 1
    bcs = read_zoneBC(zone, elts, verbose; zone_dim)
    if !isnothing(bcs)
        if validate_cell_count
            for bc in bcs
                (bc.bcdim == -1 || bc.bcdim == topo_dim - 1) || error(
                    "BC '$(bc.bcname)': selected element dimension $(bc.bcdim) must equal boundary dimension $(topo_dim - 1)",
                )
            end
        end
        filter!(bc -> (bc.bcdim == topo_dim - 1) || (bc.bcdim == -1), bcs)
    end

    # Read FlowSolutions
    if !isnothing(varnames)
        fSols = read_solutions(zone, varnames, verbose)
        if !isnothing(zone_dim) && topo_dim < zone_dim
            for solution in values(fSols)
                if solution.gridLocation == "CellCenter" && !isempty(solution.data)
                    error(
                        "FlowSolution '$(solution.name)': CellCenter arrays belong to original zone cells, not dimension $topo_dim elements",
                    )
                end
            end
        end
    else
        fSols = nothing
    end

    return (; coords, c2t, c2n, bcs, fSols, vtx_ln_to_gn, cell_ln_to_gn)
end

"""
    read_maia_unstr_numbering(zone, nvertices, cell_ids, verbose)

Return vertex and cell global numbering from the optional Maia metadata.

Without metadata, use `1:nvertices` and the selected original `cell_ids`.
Required Vertex/Cell arrays are looked up exactly and uniquely.
"""
function read_maia_unstr_numbering(zone, nvertices, cell_ids, verbose)
    # Read absolute indices : if a maia-node is found, we select it. Otherwise, we use the "erange"
    maia_gnum = get_child(
        zone;
        name = ":CGNS#GlobalNumbering",
        label = CGNS.UserDefinedData_t,
        unique = true,
    )
    if isnothing(maia_gnum)
        vtx_ln_to_gn = 1:nvertices
        cell_ln_to_gn = cell_ids
        verbose && println("Default global numbering set")
    else
        vtx_ln_to_gn = get_value(
            get_child(
                maia_gnum;
                name = "Vertex",
                label = CGNS.DataArray_t,
                required = true,
                unique = true,
            ),
        )
        cell_ln_to_gn = get_value(
            get_child(
                maia_gnum;
                name = "Cell",
                label = CGNS.DataArray_t,
                required = true,
                unique = true,
            ),
        )
        verbose && println("Maia global numbering set")
    end

    return vtx_ln_to_gn, cell_ln_to_gn
end

"""
    filter_volumic_cells!(c2t, c2n, cell_ids, elt, topo_dim)

From the `elt` info, only select cells matching the given `topo_dim`.

Specialize on the decoded section type outside the per-record loop. Homogeneous
sections can copy their arrays in bulk; heterogeneous MIXED/polyhedral sections
resolve each distinct element dimension once. Validation of every record and
section range remains the caller's responsibility, before this selection.
"""
function filter_volumic_cells!(c2t, c2n, cell_ids, elt, topo_dim)
    code = first(elt.c2t)
    if all(==(code), elt.c2t)
        if is_volumic_entity(code, topo_dim)
            append!(c2t, elt.c2t)
            append!(c2n, elt.c2n)
            append!(cell_ids, elt.erange)
        end
    else
        selected = Dict{Int, Bool}(
            code => is_volumic_entity(code, topo_dim) for code in unique(elt.c2t)
        )
        for j in eachindex(elt.c2t)
            selected[elt.c2t[j]] || continue
            push!(c2t, elt.c2t[j])
            append!(c2n, view(elt.c2n, (elt.offsets[j] + 1):elt.offsets[j + 1]))
            push!(cell_ids, first(elt.erange) + j - 1)
        end
    end
end

"""
    read_zone_struct(zone, varnames, topo_dim, usr_space_dim, verbose)

Convert a full structured zone of dimension 1, 2 or 3 to Bar2/Quad4/Hexa8
connectivity.

Coordinates, cells, fields and default absolute IDs use CGNS
Fortran ordering (i fastest).

Arguments follow `read_zone`.

Return the same arrays/BC records as the
unstructured reader. Reject inconsistent zone sizes and dimension overrides;
structured Maia global numbering is currently ignored with a warning.
"""
function read_zone_struct(zone, varnames, topo_dim, usr_space_dim, verbose)
    # Number of elements
    zonedims = get_value(zone)
    (topo_dim in 1:3 && size(zonedims) == (topo_dim, 3)) ||
        error("Structured zone dimensions must match topodim in 1:3")
    nvertices_struct = view(zonedims, 1:topo_dim, 1)
    nvertices = prod(nvertices_struct)
    ncells_struct = view(zonedims, 1:topo_dim, 2)
    ncells = prod(ncells_struct)
    verbose && println(
        "nvertices = $(nvertices_struct) ($nvertices), ncells = $(ncells_struct) ($ncells)",
    )
    (all(nvertices_struct .>= 2) && all(ncells_struct .== nvertices_struct .- 1)) ||
        error("Structured cell sizes must equal vertex sizes minus one")

    # Numbering conversion struct->unstr
    N = length(ncells_struct)
    tdims_cell = ntuple(l -> 1:ncells_struct[l], N)
    ijk2I = LinearIndices(Dims(nvertices_struct))

    # Read GridCoordinates
    coords = read_grid_coordinates(zone, topo_dim, usr_space_dim, verbose)

    # Read all elements
    if topo_dim == 3
        c2n = [
            (
                ijk2I[i, j, k],
                ijk2I[i + 1, j, k],
                ijk2I[i + 1, j + 1, k],
                ijk2I[i, j + 1, k],
                ijk2I[i, j, k + 1],
                ijk2I[i + 1, j, k + 1],
                ijk2I[i + 1, j + 1, k + 1],
                ijk2I[i, j + 1, k + 1],
            ) for (i, j, k) in Iterators.product(tdims_cell...)
        ]
        c2t = fill(BCUBE_ENTITY_TO_CGNS_ENTITY[Bcube.Hexa8_t], ncells)
    elseif topo_dim == 2
        c2n = [
            (ijk2I[i, j], ijk2I[i + 1, j], ijk2I[i + 1, j + 1], ijk2I[i, j + 1]) for
            (i, j) in Iterators.product(tdims_cell...)
        ]
        c2t = fill(BCUBE_ENTITY_TO_CGNS_ENTITY[Bcube.Quad4_t], ncells)
    elseif topo_dim == 1
        c2n = [(i, i + 1) for i in tdims_cell[1]]
        c2t = fill(BCUBE_ENTITY_TO_CGNS_ENTITY[Bcube.Bar2_t], ncells)
    else
        error("Structured file with topodim = $(topo_dim) not implemented yet")
    end
    c2n = collect(Iterators.flatten(c2n))

    # Read absolute indices
    vtx_ln_to_gn = 1:nvertices
    cell_ln_to_gn = 1:ncells
    verbose && println("Default global numbering set")
    maia_gnum = get_child(
        zone;
        name = ":CGNS#GlobalNumbering",
        label = CGNS.UserDefinedData_t,
        unique = true,
    )
    if !isnothing(maia_gnum)
        @warn "Reading maia GlobalNumbering is not implemented for structured zone yet"
    end

    # Read all BCs and then keep only the ones whose topo dim is equal to the base topo dim minus 1
    bcs = read_zoneBC(zone, nothing, verbose)
    if !isnothing(bcs)
        filter!(bc -> (bc.bcdim == topo_dim - 1) || (bc.bcdim == -1), bcs)
    end

    # Read FlowSolutions
    fSols = isnothing(varnames) ? nothing : read_solutions(zone, varnames, verbose)

    return (; coords, c2t, c2n, bcs, fSols, vtx_ln_to_gn, cell_ln_to_gn)
end

"""
    read_grid_coordinates(zone, topo_dim, usr_space_dim, verbose)

Return Cartesian coordinates as an `(nvertices, spacedim)` matrix in original
vertex order (i fastest for structured arrays).

Only exact CoordinateX/Y/Z
DataArray names in the first GridCoordinates_t are axes; consecutive axes
starting with X and at least `topo_dim` axes are required. Each must contain
ZoneVertexCount values. `usr_space_dim=0` detects embedding with
`compute_space_dim`; a positive override retains that many leading axes.

Missing grids/axes, inconsistent lengths or dimensions raise errors.

For `verbose`, see `Bcube.read_file`.

Non-Cartesian encodings are excluded.
"""
function read_grid_coordinates(zone, topo_dim, usr_space_dim, verbose)
    gridCoordinates =
        get_child(zone; label = CGNS.GridCoordinates_t, required = true, unique = true)
    axes = [
        get_child(
            gridCoordinates;
            name = "Coordinate" * suffix,
            label = CGNS.DataArray_t,
            unique = true,
        ) for suffix in ("X", "Y", "Z")
    ]
    ncoords = findlast(!isnothing, axes)
    isnothing(ncoords) && error("Missing Cartesian coordinates")
    (ncoords >= topo_dim && all(!isnothing, axes[1:ncoords])) || error(
        "Cartesian coordinates must contain consecutive X/Y/Z axes for dimension $topo_dim",
    )
    nvertices, _, _ = get_zone_dims(zone)
    arrays = [flatten_struct_to_unstr(get_value(n)) for n in axes[1:ncoords]]
    all(a -> length(a) == nvertices, arrays) ||
        error("Coordinate array lengths must match Zone vertex count $nvertices")
    coords = hcat(arrays...)
    space_dim =
        usr_space_dim > 0 ? usr_space_dim : compute_space_dim(topo_dim, coords; verbose)
    topo_dim <= space_dim <= ncoords || throw(
        ArgumentError(
            "Requested space dimension $space_dim is incompatible with $ncoords coordinate axes",
        ),
    )
    return coords[:, 1:space_dim]
end

"""
    read_connectivity(node, reshape = false)

Decode one Elements_t section with fixed-arity, MIXED (20), NGON (22) or
NFACE (23) records.

NGON/NFACE require modern zero-based ElementStartOffset;
legacy count-prefixed polygon streams are excluded. NFACE entries remain signed
face IDs until `convert_poly_sections` resolves them against NGON sections.

MIXED connectivity prefixes each record with its CGNS type; every embedded
type must map to a Bcube entity, even for non-cell or unselected records.
Reject unsupported types, truncated/trailing records and malformed ranges.
Modern ElementStartOffset arrays must equal the zero-based record starts,
including type prefixes; legacy MIXED streams without offsets are accepted.

Return `(; erange, c2t, c2n, offsets, nelts, name)`: one type per element,
flattened node IDs (signed face IDs for unresolved NFACE) without MIXED type
prefixes, and zero-based offsets into that record array (length `nelts+1`). Original CGNS IDs follow `erange` order.
`reshape=true` returns an `(nelts, arity)` matrix only for uniform arity;
heterogeneous arity raises ArgumentError. Node bounds are checked by the
owning zone reader. No input connectivity or offsets are mutated.
"""
function read_connectivity(node, reshape = false)
    # Read ElementRange and ElementConnectivity nodes
    @assert get_label(node) == CGNS.Elements_t
    name = get_name(node)
    code, _ = get_value(node)
    range_node = get_child(
        node;
        name = "ElementRange",
        label = CGNS.IndexRange_t,
        required = true,
        unique = true,
    )
    bounds = get_value(range_node)
    (length(bounds) == 2 && eltype(bounds) <: Integer && 1 <= bounds[1] <= bounds[2]) ||
        error("Section '$name': invalid ElementRange")
    erange = bounds[1]:bounds[2]
    nelts = length(erange)
    conn_node = get_child(
        node;
        name = "ElementConnectivity",
        label = CGNS.DataArray_t,
        required = true,
        unique = true,
    )
    raw = vec(get_value(conn_node))
    eltype(raw) <: Integer || error("Section '$name': connectivity must contain integers")

    # Parse connectivity
    c2t, c2n, offsets = if code in (CGNS_NGON_n, CGNS_NFACE_n)
        read_ngon_nface_connectivity(node, code, raw, nelts)
    elseif code == CGNS_MIXED
        read_mixed_connectivity(node, raw, erange)
    else
        read_fixed_connectivity(node, code, raw, nelts)
    end

    # Reshape connectivity (nelts, n_nodes_by_elt) if asked
    if reshape
        arities = diff(offsets)
        all(==(first(arities)), arities) || throw(
            ArgumentError("Section '$name': cannot reshape heterogeneous element arities"),
        )
        c2n = permutedims(Base.reshape(c2n, first(arities), nelts))
    end

    return (; erange, c2t, c2n, offsets, nelts, name)
end

"""
    read_ngon_nface_connectivity(node, code, raw, nelts)

Decode NGON vertex loops or signed NFACE references using modern offsets.

Both encodings share offset validation; face geometry and signed references
are resolved later by `convert_poly_sections`.

Input `node` must be an `Elements_t` node, `code` is the node value (element code),
and `raw` is the `ElementConnectivity` value.
"""
function read_ngon_nface_connectivity(node, code, raw, nelts)
    name = get_name(node)

    # Read offset array
    offset_node = get_child(
        node;
        name = "ElementStartOffset",
        label = CGNS.DataArray_t,
        required = true,
        unique = true,
    )
    offsets = vec(get_value(offset_node))

    # Check this offset array is correct
    (
        eltype(offsets) <: Integer &&
        length(offsets) == nelts+1 &&
        first(offsets) == 0 &&
        last(offsets) == length(raw) &&
        all(>(0), diff(offsets))
    ) || error("Section '$name': invalid NGON/NFACE ElementStartOffset")

    c2t = fill(code, nelts)
    c2n = raw

    return c2t, c2n, offsets
end

"""
    read_mixed_connectivity(node, raw, erange)

Decode every type-prefixed MIXED record, including unselected elements.

Return element types, flattened nodes and decoded offsets; optional stored
offsets must match the starts of the original type-prefixed records.

Input `node` must be an `Elements_t` node and `raw` is the `ElementConnectivity` value.
"""
function read_mixed_connectivity(node, raw, erange)
    name = get_name(node)
    c2t, c2n, offsets, starts = Int[], Int[], [0], [0]
    pos = 1
    for id in erange
        pos <= length(raw) || error("Section '$name', element $id: missing MIXED type")
        embedded = raw[pos]
        haskey(CGNS_ENTITY_TO_BCUBE_ENTITY, embedded) || error(
            "Section '$name', element $id: unsupported embedded CGNS element code $embedded",
        )
        n = nnodes(cgns_entity_to_bcube_entity(embedded))
        pos + n <= length(raw) ||
            error("Section '$name', element $id: truncated MIXED connectivity")
        push!(c2t, embedded)
        append!(c2n, view(raw, (pos + 1):(pos + n)))
        push!(offsets, length(c2n))
        pos += n + 1
        push!(starts, pos - 1)
    end
    pos == length(raw) + 1 || error("Section '$name': trailing MIXED connectivity")
    offset_node = get_child(
        node;
        name = "ElementStartOffset",
        label = CGNS.DataArray_t,
        unique = true,
    )
    if !isnothing(offset_node)
        supplied = vec(get_value(offset_node))
        (eltype(supplied) <: Integer && supplied == starts) ||
            error("Section '$name': ElementStartOffset does not match MIXED records")
    end
    return c2t, c2n, offsets
end

"""
    read_fixed_connectivity(node, code, raw, nelts)

Decode a fixed-arity section, validating its complete connectivity length.

Input `node` must be an `Elements_t` node, `code` is the node value (element code),
and `raw` is the `ElementConnectivity` value.
"""
function read_fixed_connectivity(node, code, raw, nelts)
    name = get_name(node)

    # Check element type is supported by Bcube
    haskey(CGNS_ENTITY_TO_BCUBE_ENTITY, code) ||
        error("Section '$name': unsupported CGNS element code $code")

    n_nodes_by_elt = nnodes(cgns_entity_to_bcube_entity(code))

    # Check that the total number of nodes indicated by the connectivity matches the number of elements
    # times the number of nodes by elt.
    length(raw) == nelts * n_nodes_by_elt || error(
        "Section '$name': connectivity length does not match ElementRange and element arity",
    )

    c2t = fill(code, nelts)
    c2n = raw
    offsets = collect(0:n_nodes_by_elt:length(raw))

    return c2t, c2n, offsets
end

"""
    read_index(node)

Read a one-dimensional IndexRange_t or IndexArray_t as CGNS IDs.

Ranges
include both endpoints, ascending or descending; lists return a flattened view
without sorting. Invalid range lengths and unsupported labels are rejected.
"""
function read_index(node)
    type = get_label(node)
    if type == CGNS.IndexRange_t
        erange = get_value(node)
        length(erange) == 2 || error("Unstructured index range requires two endpoints")
        return erange[1]:(erange[1] <= erange[2] ? 1 : -1):erange[2]
    elseif type == CGNS.IndexArray_t
        return vec(get_value(node))
    else
        error("Reading $type not implemented yet")
    end
end

"""
    read_zoneBC(zone, elts, verbose; zone_dim = nothing)

Read ZoneBC_t children in stored BC order; `elts` contains decoded sections
for unstructured lookup (use `nothing` for structured zones).

For `verbose`, see `Bcube.read_file`.

Return BC records from `read_bc_unstr`, or `nothing` when absent or empty. Validate
Vertex IDs against the zone vertex count. Structured selections are converted
by `read_bc_struct` using multidimensional zone indices.
"""
function read_zoneBC(zone, elts, verbose; zone_dim = nothing)
    zoneBC = get_child(zone; label = CGNS.ZoneBC_t)

    # Premature exit if no ZoneBC is present
    (isnothing(zoneBC) || isempty(get_children(zoneBC; label = CGNS.BC_t))) &&
        (return nothing)

    # Read each BC
    bcs = if read_zone_type(zone) == "Unstructured"
        map(
            bc -> read_bc_unstr(bc, elts, verbose; zone_dim),
            get_children(zoneBC; label = CGNS.BC_t),
        )
    else
        dims = Int.(get_value(zone)[:, 1])
        map(
            bc -> read_bc_struct(bc, dims, verbose),
            get_children(zoneBC; label = CGNS.BC_t),
        )
    end

    # Check that the BC nodes id are not greater that the total number of nodes in the Zone
    if !isnothing(bcs)
        nvertices, _, _ = get_zone_dims(zone)
        for bc in bcs
            all(v -> 1 <= v <= nvertices, bc.bcnodes) ||
                error("BC '$(bc.bcname)': vertex ID outside 1:$nvertices")
        end
    end

    return bcs
end

"""
    structured_bc_indices(selection, dims, bcname)

Read structured PointList columns or inclusive PointRange endpoints of shape
`(IndexDimension, n)` into tuples, i fastest for ranges.

Descending endpoints are allowed; list ordering is preserved without changing the source array.
Reject empty, noninteger, malformed or out-of-bounds selections. `dims` gives
location-specific upper bounds (vertices or directional faces).
"""
function structured_bc_indices(selection, dims, bcname)
    ids = get_value(selection)

    d = length(dims)
    (ndims(ids) == 2 && size(ids, 1) == d && size(ids, 2) > 0 && eltype(ids) <: Integer) ||
        error("BC '$bcname': expected a nonempty $d-row structured index array")

    isrange = get_label(selection) == CGNS.IndexRange_t
    (!isrange || size(ids, 2) == 2) ||
        error("BC '$bcname': PointRange requires two index columns")
    all(1 <= ids[k, j] <= dims[k] for k in 1:d for j in axes(ids, 2)) ||
        error("BC '$bcname': structured index outside location bounds")

    if isrange
        ranges = ntuple(k -> ids[k, 1]:(ids[k, 1] <= ids[k, 2] ? 1 : -1):ids[k, 2], d)
        return vec(collect(Iterators.product(ranges...)))
    else
        return [Tuple(ids[:, j]) for j in axes(ids, 2)]
    end
end

"""
    read_bc_name(bc)

Return the FamilyName when present, otherwise the BC node name.
"""
function read_bc_name(bc)
    family = get_child(bc; label = CGNS.FamilyName_t)
    return isnothing(family) ? get_name(bc) : get_value(family)
end

"""
    read_bc_location(bc)

Return the BC GridLocation string, defaulting to Vertex.
"""
function read_bc_location(bc)
    location = get_child(bc; label = CGNS.GridLocation_t)
    return isnothing(location) ? "Vertex" : parse_to_string(get_value(location))
end

"""
    read_bc_selection(bc)

Return the unique required IndexRange_t or IndexArray_t child.

The owning structured/unstructured reader decodes its indices.
"""
function read_bc_selection(bc)
    return get_child_from_predicate(
        bc,
        n -> get_label(n) in (CGNS.IndexRange_t, CGNS.IndexArray_t);
        required = true,
        unique = true,
    )
end

"""
    read_bc_struct(bc, dims, verbose)

Resolve structured Vertex BC lists/ranges (default location) or directional
IFaceCenter/JFaceCenter/KFaceCenter selections to linear vertex IDs and exact
face node sets.

CGNS face indices are the minimum connecting vertex index:
the normal axis has vertex bounds, tangential axes have cell bounds. Only
exterior faces are valid. Directional locations also accept 2D edges as a
compatibility extension. Standard 2D EdgeCenter selections infer direction
only when exactly one axis is constant at an exterior plane; ambiguous patches
are rejected. Generic FaceCenter and other locations are unsupported.

Return `(; bcname, bcnodes, bcdim, face_nodes)` for the mesh adapter; Vertex
BCs retain node-set semantics. FamilyName overrides the BC node name.
"""
function read_bc_struct(bc, dims, verbose)
    bcname = read_bc_name(bc)
    verbose && println("Reading BC '$bcname'")
    location = read_bc_location(bc)
    index = read_bc_selection(bc)
    d = length(dims)
    linear = LinearIndices(Tuple(dims))

    if location == "Vertex"
        points = structured_bc_indices(index, dims, bcname)
        return (;
            bcname,
            bcnodes = unique([linear[p...] for p in points]),
            bcdim = -1,
            face_nodes = nothing,
        )
    end

    normal = findfirst(==(location), ("IFaceCenter", "JFaceCenter", "KFaceCenter"))
    if location == "EdgeCenter" && d == 2
        points = structured_bc_indices(index, dims, bcname)
        candidates = [
            k for k in 1:d if all(p -> p[k] == first(points)[k], points) &&
                first(points)[k] in (1, dims[k])
        ]

        length(candidates) == 1 || error("BC '$bcname': ambiguous EdgeCenter direction")

        normal = only(candidates)
    end

    (!isnothing(normal) && d >= 2 && normal <= d) ||
        error("BC '$bcname': unsupported structured GridLocation '$location'")

    bounds = [k == normal ? dims[k] : dims[k] - 1 for k in 1:d]
    points = structured_bc_indices(index, bounds, bcname)
    all(p -> p[normal] in (1, dims[normal]), points) ||
        error("BC '$bcname': selected face is not exterior")
    offsets = Iterators.product((k == normal ? (0:0) : (0:1) for k in 1:d)...)
    face_nodes =
        [[linear[ntuple(k -> p[k] + o[k], d)...] for o in offsets] |> vec for p in points]

    return (; bcname, bcnodes = unique(vcat(face_nodes...)), bcdim = d - 1, face_nodes)
end

"""
    read_bc_unstr(bc, elts, verbose; zone_dim = nothing)

Resolve an unstructured BC against decoded fixed-arity, MIXED or converted
NGON/NFACE sections `elts`, using each selected record's type and node offsets.

Omitted GridLocation defaults to Vertex.

Lists preserve input order ; ranges include both endpoints in either direction.

Return `(; bcname, bcnodes, bcdim, face_nodes)`, with unique vertex IDs and, for
2D EdgeCenter and FaceCenter/legacy CellCenter selections, individual element node vectors.
`face_nodes` is a list of Vector{Int}.

`zone_dim` is the original zone dimension; EdgeCenter support is limited to 2D.

Reject missing/ambiguous selections, unknown element IDs and mixed dimensions.

Vertex bounds are validated by `read_zoneBC` against the owning zone.
"""
function read_bc_unstr(bc, elts, verbose; zone_dim = nothing)
    bcname = read_bc_name(bc)
    verbose && println("Reading BC '$bcname'")
    location = read_bc_location(bc)
    index = read_bc_selection(bc)
    ids = collect(read_index(index))
    isempty(ids) && error("BC '$bcname': empty selection")

    # If GridLocation is vertex, early exit
    if location == "Vertex"
        return (; bcname, bcnodes = unique(ids), bcdim = -1, face_nodes = nothing)
    end

    # If not, assert that grid location is FaceCenter or CellCenter
    (
        location in ("FaceCenter", "CellCenter") ||
        (location == "EdgeCenter" && zone_dim == 2)
    ) || error("BC '$bcname': unsupported GridLocation '$location'")

    # For FaceCenter elements, for each face index we must find the corresponding elements in "Elements_t"
    # Once found, we add to `face_nodes` with the face2node corresponding to this face.
    sorted_elts = sort(elts; by = elt -> first(elt.erange))
    ends = [last(elt.erange) for elt in sorted_elts]
    face_nodes = Vector{Int}[]
    dimensions = Int[]
    for id in ids
        k = searchsortedfirst(ends, id)
        (k <= length(sorted_elts) && id >= first(sorted_elts[k].erange)) ||
            error("BC '$bcname': element $id does not belong to any Elements_t node")
        elt = sorted_elts[k]
        j = id - first(elt.erange) + 1
        entity = cgns_entity_to_bcube_entity(elt.c2t[j])
        push!(dimensions, Bcube.topodim(entity))
        push!(face_nodes, Int.(elt.c2n[(elt.offsets[j] + 1):elt.offsets[j + 1]]))
    end
    all(==(first(dimensions)), dimensions) ||
        error("BC '$bcname': selected elements have mixed dimensions")

    return (;
        bcname,
        bcnodes = sort!(unique(vcat(face_nodes...))),
        bcdim = first(dimensions),
        face_nodes,
    )
end

"""
    read_solutions(zone, varnames, verbose)

Return `Dict(solution_name => read_solution(zone, fs, varnames))` for every
FlowSolution_t child in the selected zone.
"""
function read_solutions(zone, varnames, verbose)
    fSols = Dict(
        get_name(fs) => read_solution(zone, fs, varnames; verbose) for
        fs in get_children(zone; label = CGNS.FlowSolution_t)
    )

    return fSols
end

"""
    read_solution(zone, fs, varnames; verbose)

Return `(; name, gridLocation, data)` for one full-zone FlowSolution_t, where
`data` maps variable names to flat vectors; for `varnames`, see `Bcube.read_file`.

Missing GridLocation defaults to Vertex; only Vertex and CellCenter are supported.

Selected arrays must have ZoneVertexCount or ZoneCellCount entries respectively. Structured arrays flatten
i fastest; unstructured cell values retain stored order, associated with selected
cells in increasing original element ID order by the zone reader.

Partial lists/ranges, rind and unsupported locations raise errors even
when no names are requested.

Values are not converted.
"""
function read_solution(zone, fs, varnames; verbose)
    name = get_name(fs)

    # Raise an error if an unsupported type is encountered
    for label in (CGNS.IndexArray_t, CGNS.IndexRange_t, CGNS.Rind_t)
        isnothing(get_child(fs; label)) ||
            error("FlowSolution '$name': partial arrays and rind are not supported")
    end

    # Read grid location and defaults to "Vertex" it not found
    node = get_child(fs; label = CGNS.GridLocation_t)
    gridLocation = isnothing(node) ? "Vertex" : get_value(node)
    gridLocation in ("Vertex", "CellCenter") ||
        error("FlowSolution '$name': unsupported GridLocation '$gridLocation'")

    # Read zone dims
    nvertices, ncells, _ = get_zone_dims(zone)
    expected = gridLocation == "Vertex" ? nvertices : ncells

    # Read variables matching asked "varnames"
    dArrays = get_children(fs; label = CGNS.DataArray_t)
    if varnames != "*"
        # filter to obtain only the desired variables names
        filter!(dArray -> get_name(dArray) in varnames, dArrays)
    end
    for array in dArrays
        length(get_value(array)) == expected || error(
            "FlowSolution '$name', array '$(get_name(array))': expected $expected values at $gridLocation",
        )
    end
    data = Dict(
        get_name(dArray) => flatten_struct_to_unstr(get_value(dArray)) for
        dArray in dArrays
    )

    # Flow solution name
    name = get_name(fs)

    return (; name, gridLocation, data)
end

"""
    cgns_mesh_to_bcube_mesh(zoneCGNS)

Convert the NamedTuple returned by `read_zone` to a Bcube Mesh.

Coordinate rows become nodes, CGNS codes map to Bcube entities, and flat one-based connectivity
is partitioned by entity arity without changing cell/corner order. Preserve
absolute node/cell IDs. BC records may be `nothing` or an empty vector.
Vertex BCs use exterior faces whose nodes are all selected. Element BCs use
exact exterior face node sets; invalid or interior selections raise errors.
Shared family names aggregate nodes and independently resolved patch faces.
"""
function cgns_mesh_to_bcube_mesh(zoneCGNS)
    nodes = [Bcube.Node(zoneCGNS.coords[i, :]) for i in 1:size(zoneCGNS.coords, 1)]
    # nodes = map(row -> Bcube.Node(row), eachrow(zoneCGNS.coords)) # problem with Node + Slice

    c2n = Int.(zoneCGNS.c2n)
    c2t = map(cgns_entity_to_bcube_entity, zoneCGNS.c2t)
    c2nnodes = map(nnodes, c2t)

    if isnothing(zoneCGNS.bcs)
        return Bcube.Mesh(
            nodes,
            c2t,
            Bcube.Connectivity(c2nnodes, c2n);
            absoluteNodeIndices = zoneCGNS.vtx_ln_to_gn,
            absoluteCellIndices = zoneCGNS.cell_ln_to_gn,
        )
    else
        names = unique([bc.bcname for bc in zoneCGNS.bcs])
        bc_names = Dict{Int, String}(i => name for (i, name) in enumerate(names))
        bc_nodes = Dict{Int, Vector{Int}}(
            i => unique(
                Int[v for bc in zoneCGNS.bcs if bc.bcname == name for v in bc.bcnodes],
            ) for (i, name) in enumerate(names)
        )
        mesh = Bcube.Mesh(
            nodes,
            c2t,
            Bcube.Connectivity(c2nnodes, c2n);
            bc_names,
            bc_nodes,
            absoluteNodeIndices = zoneCGNS.vtx_ln_to_gn,
            absoluteCellIndices = zoneCGNS.cell_ln_to_gn,
        )
    end
    # Bcube resolves one Vertex patch per name exactly. Shared-family patches
    # must resolve separately before union to avoid inventing connecting faces.
    if all(bc -> isnothing(bc.face_nodes), zoneCGNS.bcs) &&
       allunique(bc.bcname for bc in zoneCGNS.bcs)
        return mesh
    else
        return exact_bc_faces(mesh, zoneCGNS.bcs)
    end
end

"""
    exact_bc_faces(mesh, bcs)

Rebuild boundary associations after node-based Mesh construction.

Match each explicit BC element to one exterior face by its complete sorted node IDs;
reject elements that are not mesh boundary faces. Vertex BCs retain per-patch
node-set semantics. Shared-family patches union resolved face IDs and reported
nodes, never vertex selections before resolution. First output BC wins for
overlapping selections, as in Bcube.

Return a mesh with consistent `bc_faces` and `f2bc` connectivity and unchanged geometry.
"""
function exact_bc_faces(mesh, bcs)
    isempty(bcs) && return mesh

    # Bcube builds named tuples in its dictionary iteration order, which need
    # not match CGNS child order. Tags must follow that same named-tuple order.
    names = collect(keys(mesh.bc_nodes))
    groups = [[bc for bc in bcs if Symbol(bc.bcname) == name] for name in names]
    f2n = Bcube.connectivities_indices(mesh, :f2n)
    f2c = Bcube.connectivities_indices(mesh, :f2c)
    exterior = Dict(
        Tuple(sort(collect(f2n[f]))) => f for f in 1:length(f2n) if length(f2c[f]) == 1
    )
    tags = zeros(Int, length(f2n))
    faces = [Int[] for _ in groups]
    for (tag, patches) in enumerate(groups)
        selected = Int[]
        for bc in patches
            patch_faces = if isnothing(bc.face_nodes)
                nodes = Set(bc.bcnodes)
                [f for f in values(exterior) if all(in(nodes), f2n[f])]
            else
                map(bc.face_nodes) do nodes
                    key = Tuple(sort(nodes))
                    haskey(exterior, key) || error(
                        "BC '$(bc.bcname)': selected element is not an exterior mesh face",
                    )
                    exterior[key]
                end
            end
            append!(selected, patch_faces)
        end
        for f in sort!(unique(selected))
            if tags[f] == 0
                tags[f] = tag
                push!(faces[tag], f)
            end
        end
    end
    f2bc = Bcube.MeshConnectivity(
        :face,
        :face,
        Bcube.Connectivity(Int.(tags .> 0), tags[tags .> 0]),
    )
    bc_faces = (; (name => faces[i] for (i, name) in enumerate(names))...)
    return typeof(mesh)(
        mesh.nodes,
        mesh.entities,
        merge(mesh.connectivities, (; f2bc)),
        mesh.bc_nodes,
        bc_faces,
        mesh.absoluteNodeIndices,
        mesh.absoluteCellIndices,
        mesh.metadata,
        mesh.backend,
    )
end

"""
    flow_solutions_to_bcube_data(fSols)

Wrap the validated name-keyed records from `read_solutions` as
`Dict(solution_name => Dict(variable_name => MeshData))`.

Vertex vectors become
MeshPointData and CellCenter vectors become MeshCellData, retaining vector
length, values and ordering. Assumes locations have already been validated;
empty solutions remain empty dictionaries. No interpolation is performed.
"""
function flow_solutions_to_bcube_data(fSols)
    return Dict(
        fname => Dict(
            varname =>
                t.gridLocation == "Vertex" ? MeshPointData(array) : MeshCellData(array)
            for (varname, array) in t.data
        ) for (fname, t) in fSols
    )
end

"""
    compute_space_dim(topodim, coords, tol = 1e-15; verbose = true)

Infer the retained leading Cartesian dimension from coordinate extents using
Bcube's `_compute_space_dim`.

`coords` is a nonempty `(nvertices, naxes)` matrix
with 1–3 consecutive X/Y/Z columns; absent Y/Z extents are zero. `tol` is
relative to the largest axis extent and `verbose` is forwarded to Bcube.
Return an integer at least `topodim`; this is axis-based detection, not an
arbitrary-plane rotation. Bcube rejects a 1D mesh with varying Z but constant Y;
degenerate extents are not validated here.
"""
function compute_space_dim(topodim, coords, tol = 1e-15; verbose = true)
    spacedim = size(coords, 2)

    xmin, xmax = extrema(view(coords, :, 1))
    lx = xmax - xmin

    ly = lz = 0.0

    if spacedim > 1
        ymin, ymax = extrema(view(coords, :, 2))
        ly = ymax - ymin
    end

    if spacedim > 2
        zmin, zmax = extrema(view(coords, :, 3))
        lz = zmax - zmin
    end

    return Bcube._compute_space_dim(topodim, lx, ly, lz, tol, verbose)
end

"""
    get_zone_dims(zone)

Return `(nvertices, ncells, nbnd)` from Zone_t size data.

Structured size columns
are multiplied over index directions; unstructured size values are flattened
in stored order. Assumes valid CGNS size shape; `read_zone_type` validates the
zone kind. The boundary-size value is metadata, not a count of decoded BC faces.
"""
function get_zone_dims(zone)
    dims = get_value(zone)
    if read_zone_type(zone) == "Structured"
        return Int(prod(dims[:, 1])), Int(prod(dims[:, 2])), Int(prod(dims[:, 3]))
    end
    return Int(dims[1]), Int(dims[2]), Int(dims[3])
end

"""
    read_ref_state(base)

Read the ReferenceState node in the base.

Return a Dict of (String) keys and associated values. If the ReferenceState node doesn't exist,
`nothing` is returned.
The child must be named `ReferenceState` and labelled ReferenceState_t.
Parse to depth 2 with `_recursive_parse`: singleton DataArrays become scalars,
other arrays retain shape/order, and deeper values warn and remain raw.
No dimensional conversion or solver-quantity completeness validation is applied.
"""
function read_ref_state(base)
    refState = get_child(
        base;
        name = "ReferenceState",
        label = CGNS.ReferenceState_t,
        unique = true,
    )
    return isnothing(refState) ? refState : _recursive_parse(refState, 0, 2)
end

"""
    read_grid_location_child(node)

Return the raw value of the GridLocation_t child. This low-level accessor
requires the child to exist; it neither supplies the Vertex default nor validates
supported locations. Use the BC/solution readers for those policies.
"""
function read_grid_location_child(node)
    get_value(get_child(node; label = CGNS.GridLocation_t, required = true))
end

"""
    read_zone_type(node)

Return `"Structured"` or `"Unstructured"` from the ZoneType_t child of a zone.

Decode string/byte-iterable values with `parse_to_string`. Missing or unsupported
zone kinds raise errors; no default zone type is inferred.
"""
function read_zone_type(node)
    zoneTypeNode = get_child(node; label = CGNS.ZoneType_t)
    isnothing(zoneTypeNode) && error("Zone '$(get_name(node))': missing ZoneType_t node")
    zoneType = parse_to_string(get_value(zoneTypeNode))
    zoneType in ("Structured", "Unstructured") ||
        error("Zone '$(get_name(node))': unsupported ZoneType '$zoneType'")
    return zoneType
end

"""
    parse_to_string(a::String)
    parse_to_string(iterable)

Return a String unchanged, or build a String from the first byte of each
entry in an iterable (including Python byte-array iterables).

Preserves stored
order; does not trim padding or validate CGNS names. Iterable entries must be
nonempty and their first values must form a valid String constructor input.
"""
parse_to_string(a::String) = a
parse_to_string(iterable) = String([first(c) for c in iterable])

"""
    _recursive_parse(node, depth, max_depth)

Parse ReferenceState descendants into name-keyed dictionaries.

DataArrays with
one entry become scalars; longer arrays retain shape/order. Descriptors and
nodes carrying DimensionalExponents return their raw value. At `max_depth`,
warn and return the current raw value instead of descending. Callers supply
zero-based `depth`; `read_ref_state` starts at 0 with limit 2. Assumes nonempty
DataArrays and does not convert units or enforce named solver quantities.
"""
function _recursive_parse(node, depth, max_depth)
    # If it's a DataArray, return the value (scalar if necessary)
    if get_label(node) == CGNS.DataArray_t
        x = get_value(node)
        return length(x) > 1 ? x : first(x)
    end

    if (get_label(node) == CGNS.Descriptor_t) ||
       has_child(node; label = CGNS.DimensionalExponents_t)
        return get_value(node)
    end

    if depth == max_depth
        @warn "Reached max depth when parsing ReferenceState, some entries might be missing"
        return get_value(node)
    end

    Dict(
        get_name(child) => _recursive_parse(child, depth + 1, max_depth) for
        child in get_children(node)
    )
end

"""
    flatten_struct_to_unstr(array)

Return `vec(array)` in Julia/CGNS column-major order (i fastest), preserving
values and sharing storage when `vec` permits it.

No sorting or copying is
promised; the input must support `vec`.
"""
flatten_struct_to_unstr(array) = vec(array)

"""
    validate_poly_face(nodes, coords, context)

Validate a distinct-vertex TRI/QUAD loop against zone coordinates.

Require a nondegenerate convex planar polygon; preserve its supplied cyclic orientation.
Return its unnormalized three-component normal from the first three corners.
Two-dimensional coordinates are padded with zero Z. `nodes` contains one-based
row IDs into the vertex-by-axis matrix `coords`. Errors include
`context`; tolerance is 1e-10 scaled by the face diameter.
"""
function validate_poly_face(nodes, coords, context)
    length(nodes) in (3, 4) || error("$context: arbitrary polygons are unsupported")
    length(unique(nodes)) == length(nodes) || error("$context: repeated vertex")
    all(v -> 1 <= v <= size(coords, 1), nodes) || error("$context: invalid vertex ID")

    p = [vcat(collect(coords[v, :]), zeros(3-size(coords, 2))) for v in nodes]

    all(x -> all(isfinite, x), p) || error("$context: nonfinite coordinates")

    scale = maximum(sqrt(sum(abs2, x-p[1])) for x in p)
    normal = cross(p[2]-p[1], p[3]-p[1])
    tol = 1e-10*scale^2
    sum(abs2, normal) > tol^2 || error("$context: degenerate face")

    for i in eachindex(p)
        abs(sum(normal .* (p[i]-p[1]))) <= tol*scale || error("$context: nonplanar face")

        turn = cross(
            p[mod1(i+1, length(p))]-p[i],
            p[mod1(i+2, length(p))]-p[mod1(i+1, length(p))],
        )

        sum(turn .* normal) > tol^2 || error("$context: nonconvex or crossed face")
    end
    return normal
end

"""
    poly_cycle_equal(a, b)

Compare oriented node loops up to cyclic rotation, never reflection.

Return Bool; `a` and `b` are one-based indexed node-ID vectors. Unequal
lengths do not match; geometry and duplicate IDs are validated by the caller.
"""
function poly_cycle_equal(a, b)
    length(a)==length(b) &&
        any(k -> all(i -> a[i]==b[mod1(i+k, length(b))], eachindex(a)), 0:(length(b) - 1))
end

"""
    poly_assignment(template, faces, vertices, assignment)

Extend a partial canonical node assignment by backtracking over sorted vertex
IDs.

Prune every completed face against oriented supplied loops. Return a full
node vector, or `nothing` when no extension exists; restore the input vector.
`template` contains one-based canonical corner indices, `faces` contains oriented
zone-node loops, `vertices` is the candidate ID vector and `assignment` is a
mutable prefix of the canonical node vector. This helper checks topology only.
"""
function poly_assignment(template, faces, vertices, assignment)
    for f in template
        maximum(f) <= length(assignment) || continue
        any(g -> poly_cycle_equal(assignment[collect(f)], g), faces) || return nothing
    end
    length(assignment) == length(vertices) && return copy(assignment)
    for v in vertices
        v in assignment && continue
        push!(assignment, v)
        result = poly_assignment(template, faces, vertices, assignment)
        pop!(assignment)
        isnothing(result) || return result
    end
    return nothing
end

"""
    match_poly_cell(faces, context)

Find a deterministic orientation-preserving node assignment to Bcube's linear
TETRA_4, PYRA_5, PENTA_6 or HEXA_8 face templates.

Return `(code, nodes)`.
Backtracking prunes completed faces; all supplied faces must match cyclically.
Reject shells which share counts but not the representable incidence topology.
"""
function match_poly_cell(faces, context)
    vertices = sort!(unique(vcat(faces...)))
    for code in (10, 12, 14, 17)
        entity = cgns_entity_to_bcube_entity(code)
        nnodes(entity) == length(vertices) || continue
        template = Bcube.faces2nodes(entity)
        sort(collect(length.(template))) == sort(length.(faces)) || continue
        result = poly_assignment(template, faces, vertices, Int[])
        isnothing(result) || return code, result
    end
    error("$context: unsupported polyhedron topology or inconsistent face orientation")
end

"""
    convert_poly_sections(elts, coords)

Resolve modern NGON vertex loops and signed NFACE face references across all
sections before dimensional selection.

Positive face IDs use the stored outward
loop, negative IDs reverse it. Faces must be simple convex planar triangles or
quads. Cell shells must match linear tetrahedra, pyramids, prisms or hexahedra,
have positive volume and be convex. Repeated/missing faces, open or inconsistent
shells, inward normals, nonmanifold reuse and same-sign shared faces raise errors.
NGON-only 2D cells retain cyclic order (positive XY orientation when embedded in
2D). Preserve original element ranges and compact sorted cell-field association;
converted NGON faces remain available for exact BC selection.
    
Return a vector of section records with the same fields as `read_connectivity`; converted types,
node vectors and zero-based offsets replace raw polygon/face references.
`coords` is the vertex-by-axis matrix; errors include section/ID
context.
"""
function convert_poly_sections(elts, coords)
    # Check that elements contains at least one NGON, ie at least one no-canonical element.
    # Otherwise, no need for conversion.
    any(e -> any(c -> c in (CGNS_NGON_n, CGNS_NFACE_n), e.c2t), elts) || return elts

    faces = Dict{Int, Vector{Int}}()
    usage = Dict{Int, Vector{Int}}()
    has_nface = any(e -> CGNS_NFACE_n in e.c2t, elts)

    for e in elts
        CGNS_NGON_n in e.c2t || continue
        for j in 1:e.nelts
            id = first(e.erange)+j-1
            nodes = collect(e.c2n[(e.offsets[j] + 1):e.offsets[j + 1]])

            normal = validate_poly_face(nodes, coords, "Section '$(e.name)', element $id")
            !has_nface &&
                size(coords, 2)==2 &&
                normal[3] <= 0 &&
                error("NGON element $id: nonpositive orientation")
            haskey(faces, id) && error("Overlapping NGON element ID $id")

            faces[id] = nodes
        end
    end

    result = map(elts) do e
        any(c -> c in (CGNS_NGON_n, CGNS_NFACE_n), e.c2t) || return e

        types, nodes, offsets = Int[], Int[], [0]
        for j in 1:e.nelts
            id = first(e.erange)+j-1
            context = "Section '$(e.name)', element $id"
            if e.c2t[j]==CGNS_NGON_n
                loop = faces[id]
                code = length(loop)==3 ? 5 : 7
            else
                refs = collect(e.c2n[(e.offsets[j] + 1):e.offsets[j + 1]])
                length(unique(abs.(refs)))==length(refs) ||
                    error("$context: repeated face")
                loops = Vector{Int}[]
                for ref in refs
                    ref != 0 && haskey(faces, abs(ref)) ||
                        error("$context: missing NGON face $ref")
                    push!(get!(usage, abs(ref), Int[]), sign(ref))
                    push!(loops, ref>0 ? faces[ref] : reverse(faces[-ref]))
                end
                code, loop = match_poly_cell(loops, context)
                size(coords, 2)==3 || error("$context: NFACE requires 3D coordinates")
                points = [collect(coords[v, :]) for v in loop]
                center = sum(points)/length(points)
                scale = maximum(sqrt(sum(abs2, p-center)) for p in points)
                for f in loops
                    normal = validate_poly_face(f, coords, context)
                    origin = collect(coords[first(f), :])
                    sum(normal .* (center-origin)) < -1e-10*scale^3 ||
                        error("$context: inward or degenerate shell")
                    all(p -> sum(normal .* (p-origin)) <= 1e-10*scale^3, points) ||
                        error("$context: nonconvex shell")
                end
            end
            push!(types, code)
            append!(nodes, loop)
            push!(offsets, length(nodes))
        end
        return merge(e, (c2t = types, c2n = nodes, offsets = offsets))
    end
    for (id, signs) in usage
        (length(signs)<=2 && (length(signs)<2 || sum(signs)==0)) ||
            error("NGON face $id: nonmanifold or inconsistent shared orientation")
    end
    return result
end
