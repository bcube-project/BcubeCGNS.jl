"""
    extract_surf_from_eclipps(
        filepath::String,
        bcnames;
        verbose = false,
        timestep = nothing,
    )

Extract one named BC from a single-zone CGNS-ECLIPPS file as
`(; mesh, data, refState)`.

A string or singleton collection selects the BC by
node name. `timestep` follows the CGNS/HDF5 `Bcube.read_file` contract, including
`:latest` selection and its warning when iterative metadata is absent.
No complete solver-input set is required:
reference quantities and all named Dirichlet/Neumann arrays retain their values.

Surface cells follow `SurfacicElementList/List` order and stored orientation.
Unstructured lists contain original element IDs; structured lists contain
IndexDimension-tuples for each Bar2/Quad4 corner (2D/3D respectively).
Structured corners must form an elementary exterior logical edge/face with a
valid cyclic order; both orientations are retained independently of embedding.
These encodings are compatibility assumptions from existing files, not a
complete specification in Montreuil et al. (2009), Figures 3–5.
Vertex BC nodes retain PointList/expanded PointRange order, including unused
and repeated nodes (connectivity uses the last occurrence, for compatibility); all surface vertices must be present. FaceCenter BCs compact referenced
vertices (sorted for unstructured, first occurrence for structured).

`data[dataset][variable]` contains MeshPointData or MeshCellData according to
the dataset location (inherited from the BC when absent), never guessed from
array length. Explicit dataset selections must cover the complete extracted
node/cell set and are remapped to mesh order; partial datasets are unsupported.
Scalar BC arrays broadcast over that set. Structured FaceCenter indices identify
the minimum connecting vertex tuple of each surface-list face; ambiguous keys
are rejected. Different locations require an explicit selection. Duplicate names across Dirichlet
and Neumann groups, invalid lengths/indices and unsupported locations raise
errors. General face-centered FlowSolutions and multizone assembly are outside
this API.
"""
function extract_surf_from_eclipps(
    filepath::String,
    bcnames;
    verbose = false,
    timestep = nothing,
)
    tree = h5open(CGNS.hdf5_to_node, filepath, "r")
    return extract_eclipps_tree(tree, bcnames; verbose, timestep)
end

"""
    extract_surf_from_eclipps(filepath::String, bcname::String; kwargs...)

Extract one BC selected by a string; see `extract_surf_from_eclipps`.
"""
function extract_surf_from_eclipps(filepath::String, bcname::String; kwargs...)
    extract_surf_from_eclipps(filepath, (bcname,); kwargs...)
end

"""
    eclipps_location(node, default = "Vertex")

Return the GridLocation string, inheriting `default` when absent.
"""
function eclipps_location(node, default = "Vertex")
    loc = get_child(node; label = CGNS.GridLocation_t)
    return isnothing(loc) ? default : parse_to_string(get_value(loc))
end

"""
    eclipps_selection(node; required = true)

Return the sole PointList/PointRange child, or nothing for an optional missing
selection.

Reject ambiguous selections; nested user-defined lists are ignored.
"""
function eclipps_selection(node; required = true)
    return get_child_from_predicate(
        node,
        n -> get_name(n) in ("PointList", "PointRange");
        unique = true,
        required = required,
    )
end

"""
    eclipps_indices(selection, dims; allow_repeated = false)

Resolve ordered positive integer IDs, expanding inclusive descending
ranges.

`dims=nothing` means unstructured IDs; otherwise validate structured
IndexDimension-tuples and convert to i-fastest linear vertex IDs. Repeated IDs
are rejected unless `allow_repeated=true` for historical Vertex BC lists.
"""
function eclipps_indices(selection, dims; allow_repeated = false)
    a = get_value(selection)

    (a isa AbstractArray && eltype(a) <: Integer) ||
        error("Eclipps indices must be integer arrays")

    ids = if isnothing(dims)
        collect(read_index(selection))
    else
        linear = LinearIndices(Tuple(dims))
        [linear[p...] for p in structured_bc_indices(selection, dims, get_name(selection))]
    end

    (!isempty(ids) && all(>(0), ids) && (allow_repeated || allunique(ids))) ||
        error("Eclipps selection must contain distinct positive indices")

    return Int.(ids)
end

"""
    extract_eclipps_tree(tree, bcnames; verbose = false, timestep = nothing)

Entry point for one named BC.

Validate the entire unstructured connectivity before selecting surface cells, retaining whole-MIXED
rejection. Local mesh IDs are artificial contiguous one-based IDs, as in the
historical public API. BC point order and surface-list order are independent.

Return `(; mesh, data, refState)` with the same wrappers, association and
limitations as `extract_surf_from_eclipps`. For `verbose` and `timestep`, see
`Bcube.read_file`. 

Missing/ambiguous base, zone, BC or surface lists, unsupported dimensions/locations and invalid
surface connectivity raise errors; no writer or solver completeness is required.
"""
function extract_eclipps_tree(tree, bcnames; verbose = false, timestep = nothing)
    # The argument `bcnames` can be either an AbstractString or a collection,
    # but for now only one bcname is supported
    names = bcnames isa AbstractString ? [bcnames] : collect(bcnames)
    length(names) == 1 ||
        throw(ArgumentError("Eclipps extraction requires exactly one BC name"))

    # Read the base and check the dimension
    base = get_cgns_base(tree)
    topodim, spacedim = get_value(base)
    topodim in (2, 3) || error("Eclipps extraction supports 2D/3D zones")

    # Select the zone : only one zone supported for now
    zone =
        select_timestep(base, get_child(base; label = CGNS.Zone_t, unique = true), timestep)

    # Read bcs matching the name, only one supported for now
    zonebc = get_child(zone; label = CGNS.ZoneBC_t, required = true, unique = true)
    bcs = filter(n -> get_name(n) == only(names), get_children(zonebc; label = CGNS.BC_t))
    length(bcs) == 1 ||
        throw(ArgumentError("Selected BC '$(only(names))' not found uniquely"))
    bc = only(bcs)

    # Grid location
    location = eclipps_location(bc)
    location in ("Vertex", "FaceCenter") ||
        error("Unsupported Eclipps BC GridLocation '$location'")

    # Read SurfacicElementList
    surfacicEltList = get_child(
        bc;
        name = "SurfacicElementList",
        label = CGNS.UserDefinedData_t,
        required = true,
        unique = true,
    )
    list = get_child(
        surfacicEltList;
        name = "List",
        label = CGNS.DataArray_t,
        required = true,
        unique = true,
    )
    raw = get_value(list)
    (raw isa AbstractArray && eltype(raw) <: Integer && !isempty(raw)) ||
        error("SurfacicElementList must be a nonempty integer array")

    is_structured = read_zone_type(zone) == "Structured"
    dims = is_structured ? Int.(get_value(zone)[:, 1]) : nothing
    geometry = if is_structured
        read_eclipps_struct_geometry(
            zone,
            raw,
            dims,
            topodim,
            spacedim,
            location,
            get_name(bc),
            verbose,
        )
    else
        read_eclipps_unstr_geometry(zone, raw, topodim, spacedim, verbose)
    end
    (; coords, c2t, cells, cellids) = geometry

    # Checks
    all(c -> allunique(c), cells) || error("Surface element repeats a vertex")
    allunique(Tuple.(sort.(cells))) || error("SurfacicElementList repeats an element")

    selection = eclipps_selection(bc)
    pointids = if location == "Vertex"
        eclipps_indices(selection, dims; allow_repeated = true)
    elseif is_structured
        # Preserve the historical first-occurrence structured node compaction.
        unique(vcat(cells...))
    else
        sort!(unique(vcat(cells...)))
    end
    all(i -> i <= size(coords, 1), pointids) || error("BC vertex index outside zone")

    node_map = Dict(id => i for (i, id) in enumerate(pointids))
    all(id -> haskey(node_map, id), vcat(cells...)) ||
        error("BC vertex selection omits surface vertices")

    if location == "FaceCenter"
        ids = eclipps_indices(selection, dims)
        Set(ids) == Set(cellids) ||
            error("BC face selection differs from SurfacicElementList")
    end

    zoneCGNS = (;
        coords = coords[pointids, :],
        c2t,
        c2n = [node_map[id] for id in vcat(cells...)],
        bcs = nothing,
        fSols = nothing,
        vtx_ln_to_gn = 1:length(pointids),
        cell_ln_to_gn = 1:length(cells),
    )

    # Build the Bcube mesh
    mesh = cgns_mesh_to_bcube_mesh(zoneCGNS)

    # Read the dataset
    data = eclipps_datasets(bc, location, pointids, cellids, dims)

    return (; mesh, data, refState = read_ref_state(base))
end

"""
    read_eclipps_struct_geometry(
        zone,
        raw,
        dims,
        topodim,
        spacedim,
        location,
        bcname,
        verbose,
    )

Decode structured surface corner tuples from `SurfacicElementList`.

Return `(; coords, c2t, cells, cellids)` in supplied surface order and orientation.
`dims` contains vertex dimensions; face IDs use minimum connecting vertex tuples.
Validate elementary exterior topology before constructing the mesh.
"""
function read_eclipps_struct_geometry(
    zone,
    raw,
    dims,
    topodim,
    spacedim,
    location,
    bcname,
    verbose,
)
    length(dims) == topodim || error("Structured IndexDimension does not match base")

    coords = read_grid_coordinates(zone, topodim, spacedim, verbose)

    # Number of nodes by elements : quadrangles (4 nodes) in 3D, lines (2 nodes) in 2D
    n_nodes_by_elt = topodim == 3 ? 4 : 2
    length(raw) % (n_nodes_by_elt * topodim) == 0 ||
        error("Malformed structured SurfacicElementList length")

    tuples = reshape(vec(raw), topodim, :)
    indexnode = CGNS.CGNSNode("PointList", tuples, CGNS.CGNSNode[], "IndexArray_t")

    # Corners shared by adjacent elements may repeat, unlike point selections.
    points = structured_bc_indices(indexnode, dims, bcname)
    validate_structured_surface(points, dims, bcname)
    linear = LinearIndices(Tuple(dims))
    conn = reshape([linear[p...] for p in points], n_nodes_by_elt, :)
    cells = [collect(c) for c in eachcol(conn)]
    ctype = bcube_entity_to_cgns_entity(topodim == 3 ? Bcube.Quad4_t() : Bcube.Bar2_t())
    c2t = fill(ctype, length(cells))
    # The corner tuples identify each structured face by its minimum
    # connecting vertex index; reject patches where that key is ambiguous.
    cartesian = CartesianIndices(Tuple(dims))
    cellids = [
        linear[ntuple(k -> minimum(cartesian[id][k] for id in c), topodim)...] for
        c in cells
    ]

    (location != "FaceCenter" || allunique(cellids)) ||
        error("Ambiguous structured surface face indices")
    return (; coords, c2t, cells, cellids)
end

"""
    read_eclipps_unstr_geometry(zone, raw, topodim, spacedim, verbose)

Decode surface element IDs from `SurfacicElementList` after full zone validation.

Return `(; coords, c2t, cells, cellids)` in requested order, using original zone
IDs rather than optional Maia global numbering. The input tree is not mutated.
"""
function read_eclipps_unstr_geometry(zone, raw, topodim, spacedim, verbose)
    # Generic FaceCenter structured BC conventions must not leak into the
    # full-zone BC reader. Unstructured sections are all decoded here.
    # Surface lists refer to zone element IDs, never Maia global IDs.
    # A shallow zone view omits that optional metadata without mutating input.
    local_zone = CGNS.CGNSNode(
        get_name(zone),
        get_value(zone),
        get_children_from_predicate(zone, n -> get_name(n) != ":CGNS#GlobalNumbering"),
        "Zone_t",
    )
    z = read_zone_unstr(
        local_zone,
        nothing,
        topodim-1,
        spacedim,
        verbose;
        zone_dim = topodim,
    )
    coords = z.coords
    cellids = Int.(vec(raw))
    c2t, cells = select_eclipps_cells(z, cellids)
    return (; coords, c2t, cells, cellids)
end

"""
    select_eclipps_cells(zone, cellids)

Return types and corner lists for the requested original surface element IDs.

`zone` comes from complete `read_zone_unstr` validation with optional global
numbering omitted, so its original cell IDs are sorted. Binary search preserves
the requested order without materializing corner lists for unrequested faces.
Duplicate or absent selections raise errors.
"""
function select_eclipps_cells(zone, cellids)
    positions = map(id -> searchsortedfirst(zone.cell_ln_to_gn, id), cellids)
    (
        allunique(cellids) && all(
            j <= length(zone.cell_ln_to_gn) && zone.cell_ln_to_gn[j] == id for
            (j, id) in zip(positions, cellids)
        )
    ) || error("SurfacicElementList contains duplicate or unknown surface element IDs")

    offsets = Vector{Int}(undef, length(zone.c2t) + 1)
    offsets[1] = 0
    for j in eachindex(zone.c2t)
        code = zone.c2t[j]
        n_nodes_by_elt = nnodes(cgns_entity_to_bcube_entity(code))
        offsets[j + 1] = offsets[j] + n_nodes_by_elt
    end
    cells = [collect(view(zone.c2n, (offsets[j] + 1):offsets[j + 1])) for j in positions]

    return zone.c2t[positions], cells
end

"""
    validate_structured_surface(points, dims, bcname)

Validate ordered structured Bar2/Quad4 corners in logical index space.

Each
element must have distinct complete corners of one elementary exterior edge
or face: exactly one axis fixed at an exterior plane and unit tangential spans.
Consecutive quad corners must share an edge, rejecting crossed cycles. Both
orientations are accepted unchanged; embedding coordinates are irrelevant.
`points` are already bounds-checked tuples, `dims` are vertex dimensions.
Return nothing or raise a contextual error before building a mesh.
"""
function validate_structured_surface(points, dims, bcname)
    d = length(dims)
    arity = d == 2 ? 2 : 4
    for (i, corners) in enumerate(Iterators.partition(points, arity))
        lo = [minimum(p[k] for p in corners) for k in 1:d]
        hi = [maximum(p[k] for p in corners) for k in 1:d]
        fixed = findall(==(0), hi .- lo)
        valid = allunique(corners) && length(fixed) == 1
        if valid
            normal = only(fixed)
            valid =
                lo[normal] in (1, dims[normal]) &&
                all(k -> k == normal || hi[k] - lo[k] == 1, 1:d) &&
                all(
                    j ->
                        sum(
                            abs(corners[j][k] - corners[mod1(j+1, arity)][k]) for k in 1:d
                        ) == 1,
                    1:arity,
                )
        end
        valid || error(
            "BC '$bcname': structured surface element $i must be an elementary exterior edge/face with cyclic corners",
        )
    end
    return nothing
end

"""
    eclipps_datasets(bc, location, pointids, cellids, dims)

Read named Dirichlet/Neumann arrays into the historical flat dataset dictionaries.

Resolve inherited or explicit locations; explicit selections must be permutations
of the complete target set. Scalar arrays broadcast, local arrays must match
selection length. Preserve numeric values and reject name collisions rather than
silently losing one equation group. `location` is the parent BC location;
`pointids` and `cellids` give original IDs in extracted mesh order. `dims` is
`nothing` for unstructured IDs or the structured vertex dimensions. Return
`Dict{String,Dict{String,Bcube.MeshData}}`, with flat MeshPointData/MeshCellData
vectors at Vertex/FaceCenter respectively. Changed locations require explicit
selections. Invalid/non-numeric arrays, incomplete selections, unsupported
locations and duplicate dataset/variable names raise errors. Empty groups
are allowed. No solver completeness checks are imposed.
"""
function eclipps_datasets(bc, location, pointids, cellids, dims)
    data = Dict{String, Dict{String, Bcube.MeshData}}()
    for dataset in get_children(bc; label = CGNS.BCDataSet_t)
        name = get_name(dataset)
        haskey(data, name) && error("Duplicate BC dataset '$name'")
        loc = eclipps_location(dataset, location)
        loc in ("Vertex", "FaceCenter") ||
            error("Dataset '$name': unsupported GridLocation '$loc'")
        target = loc == "Vertex" ? pointids : cellids
        (loc != "FaceCenter" || allunique(cellids)) ||
            error("Dataset '$name': ambiguous structured surface face indices")
        selection = eclipps_selection(dataset; required = false)
        source = if !isnothing(selection)
            eclipps_indices(selection, dims)
        elseif loc != location
            error("Dataset '$name': changed location requires explicit point selection")
        elseif loc == "Vertex"
            target
        else
            eclipps_indices(eclipps_selection(bc), dims)
        end
        (length(source) == length(target) && Set(source) == Set(target)) ||
            error("Dataset '$name': selection must cover the complete extracted $(loc) set")
        positions = Dict(id => i for (i, id) in enumerate(source))
        permutation =
            source == target ? collect(eachindex(target)) : [positions[id] for id in target]
        # Bcube/VTK dispatches on dictionary value types, not just contents.
        # Keep a MeshData bound while allowing different numeric array types.
        variables = Dict{String, Bcube.MeshData}()
        for groupname in ("DirichletData", "NeumannData")
            group =
                get_child(dataset; name = groupname, label = CGNS.BCData_t, unique = true)
            isnothing(group) && continue
            for array in get_children(group; label = CGNS.DataArray_t)
                variable = get_name(array)
                haskey(variables, variable) && error(
                    "Dataset '$name': duplicate variable '$variable' across BC data groups",
                )
                values = get_value(array)
                (values isa AbstractArray && eltype(values) <: Number) ||
                    error("Dataset '$name/$variable': expected numeric array")
                n = length(values)
                n in (1, length(source)) || error(
                    "Dataset '$name/$variable': array length $n does not match selection length $(length(source))",
                )
                ordered =
                    n == 1 ? fill(only(values), length(target)) : vec(values)[permutation]
                variables[variable] =
                    loc == "Vertex" ? MeshPointData(ordered) : MeshCellData(ordered)
            end
        end
        data[name] = variables
    end
    return data
end
