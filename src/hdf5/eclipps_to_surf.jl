"""
    extract_surf_from_eclipps(filepath::String, bcnames; verbose = false)

Extract Bcube compatible mesh and data from a list of BC of a CGNS-ECLIPPS file.

Warning : for now it is limited to only one BC because Bcube doesn't support multi zones.

# Dev notes:
    * todo : add an option to retrieve boundary conditions of the extracted surface(s)
    * find a way to factorize the beginning of the code with `read_file`
"""
function extract_surf_from_eclipps(filepath::String, bcnames; verbose = false)
    # Open the file
    file = h5open(filepath, "r")
    root = file

    # Read (unique) CGNS base
    cgnsBase = get_cgns_base(root)

    # Reference state
    refState = read_ref_state(cgnsBase)

    # Read base dimensions (topo and space)
    topodim, spacedim = get_value(cgnsBase)
    verbose && println("topodim = $topodim, spacedim = $(spacedim)")

    # Find the list of Zone_t
    zones = get_children(cgnsBase; type = "Zone_t")
    if length(zones) == 0
        error("Could not find any Zone_t node in the file")
    elseif length(zones) > 1
        error("The file contains several Zone_t nodes, only one zone is supported for now")
    end
    zone = first(zones)

    # Read zone
    coords, c2t, c2n, bcs, _ = read_zone(zone, nothing, topodim - 1, spacedim, 0, verbose)

    # Reshape c2n into a vector of vectors
    offsets = cumsum(nnodes.(cgns_entity_to_bcube_entity.(c2t))) .+ 1
    prepend!(offsets, 1)
    c2n = [view(c2n, offsets[i]:(offsets[i + 1] - 1)) for i in 1:length(c2t)]

    # Read ZoneBC and filter the ones we are interested in
    zoneBC = get_child(zone; type = "ZoneBC_t")
    bcs = get_children(zoneBC; type = "BC_t")
    filter!(z -> get_name(z) in bcnames, bcs)
    bcs = map(_read_eclipps_bc, bcs)
    @assert length(bcs) == 1
    bc = first(bcs)

    # Identify the subset of required nodes
    node2selected = zeros(Bool, size(coords, 1))
    for ielt in bc.ielts
        node2selected[c2n[ielt]] .= true
    end
    new2old_nodes = findall(node2selected)
    n_new = length(new2old_nodes)
    verbose && println("Number of nodes in extracted surface : $(length(new2old_nodes))")
    old2new_nodes = zeros(Int, size(coords, 1))
    old2new_nodes[new2old_nodes] .= 1:n_new

    # Build the new c2n
    _c2n = map(bc.ielts) do ielt
        old2new_nodes[c2n[ielt]]
    end
    _c2n = reduce(vcat, _c2n)

    # Build the Bcube mesh
    new_ZoneBC = (;
        coords = view(coords, new2old_nodes, :),
        c2t = view(c2t, bc.ielts),
        c2n = _c2n,
        bcs = nothing,
        fSols = nothing,
    )
    mesh = cgns_mesh_to_bcube_mesh(new_ZoneBC)

    # Shape data for Bcube : in the CGNS-ECLIPPS convention, data are
    # node centered
    #- merge dirichlet and neumann
    _data = map(bc.data) do data
        d = vcat(data.dirichlet, data.neumann)
        x = [v for v in values(d) if !isnothing(v)]
        (data.name, x)
    end
    data = Dict(
        fname => Dict(name => if (length(array) == n_new)
            MeshPointData(array)
        else
            MeshCellData(array)
        end for (name, array) in d) for (fname, d) in _data
    )

    # Close the file
    close(file)

    return (; mesh, data, refState)
end
function extract_surf_from_eclipps(filepath::String, bcname::String; kwargs...)
    extract_surf_from_eclipps(filepath, (bcname,); kwargs...)
end

"""
Specificity of ECLIPPS BC is that is contains a list of elements (to
be picked in the zone connectivity) defining the BC.
"""
function _read_eclipps_bc(bc)
    # List of surfacic elements (to be picked in zone connectivity)
    # Rq : it is often stored as a (1,n) array in CGNS so we vectorize
    ielts = vec(
        get_value(
            get_child(
                get_child(bc; name = "SurfacicElementList");
                name = "List",
                type = "DataArray_t",
            ),
        ),
    )

    bcDataSets = get_children(bc; type = "BCDataSet_t")
    filter!(
        bcs ->
            has_child(bcs; name = "DirichletData", type = "BCData_t") ||
            has_child(bcs; name = "NeumannData", type = "BCData_t"),
        bcDataSets,
    )
    data = map(_read_eclipps_bcdataset, bcDataSets)

    return (; bcname = get_name(bc), ielts, data)
end

"""
Return a NamedTuple with two entries : "dirichlet" and "neumann". For
each entry, a vector of (name,value) for each "data" (ex: temperature)
"""
function _read_eclipps_bcdataset(bcs)
    bcData = get_child(bcs; name = "DirichletData", type = "BCData_t")
    dirichlet = if isnothing(bcData) || !has_child(bcData; type = "DataArray_t")
        nothing
    else
        map(get_children(bcData; type = "DataArray_t")) do data
            (name = get_name(data), value = get_value(data))
        end
    end

    bcData = get_child(bcs; name = "NeumannData", type = "BCData_t")
    neumann = if isnothing(bcData) || !has_child(bcData; type = "DataArray_t")
        nothing
    else
        map(get_children(bcData; type = "DataArray_t")) do data
            (name = get_name(data), value = get_value(data))
        end
    end

    return (; name = get_name(bcs), dirichlet, neumann)
end
