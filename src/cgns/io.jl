function node_to_hdf5(parent, node::CGNSNode{L}) where {L}
    g = create_group(parent, get_name(node))

    attrs = attributes(g)
    attrs["label"] = enum_to_string(L)

    if node.data !== nothing
        g[" data"] = get_value(node)
    end

    for child in get_children(node)
        node_to_hdf5(g, child)
    end
end

function hdf5_to_node(g)
    attrs = attributes(g)

    name = read(attrs["name"])
    label = read(attrs["label"])

    data = read_data(g)

    cgnsChildrenNames = findall(is_cgns_node, g)
    children = if length(cgnsChildrenNames) > 0
        map(k -> hdf5_to_node(g[k]), cgnsChildrenNames)
    else
        CGNSNode{Empty, Nothing, Nothing}[]
    end

    # Root node get a special treatment
    if HDF5.name(g) == "/"
        name = "CGNSTree"
        label = CGNS.enum_to_string(CGNS.CGNSTree_t)
    end
    return CGNSNode(name, data, children, label)
end

function read_data(obj)
    haskey(obj, " data") || return nothing
    data_type = read(attributes(obj)["type"])
    data = read(obj[" data"])
    if (data_type == "C1") && (data isa AbstractVector)
        return String(UInt8.(data))
    elseif (data_type == "C1") && (data isa AbstractArray)
        x = extract_strings(data)
        return x
    elseif data_type in ("I4", "I8", "R4", "R8")
        return data
    else
        error("Datatype '$(data_type)' with shape '$(typeof(data))' not handled")
    end
end

function extract_strings(array)
    dims = size(array)
    tail_dims = dims[2:end]
    result = Array{String}(undef, tail_dims...)

    for I in CartesianIndices(tail_dims)
        full_index = (Colon(), Tuple(I)...)
        bytes = UInt8.(array[full_index...])

        # Find first 0x00
        pos = findfirst(==(0x00), bytes)

        if pos === nothing
            result[I] = String(bytes)
        else
            result[I] = String(@view bytes[1:(pos - 1)])
        end
    end

    return result
end

""" WARNING : this is not a CGNS compliant file yet, it's only for debug """
function write_cgns_file(filename::String, root::CGNSNode)
    h5open(filename, "w") do file
        write_cgns_node(file, root)
    end
end

function read_cgns_file(filename::String)
    file = h5open(filename, "r")

    # Find CGNSBase
    baseKeys = findall(
        n -> is_cgns_node(n) && (read(attributes(n)["label"]) == "CGNSBase_t"),
        file,
    )
    bases = map(key -> hdf5_to_node(file[key]), baseKeys)
    if length(bases) == 1
        close(file)
        return first(bases)
    elseif length(baseKeys) > 1
        close(file)
        return bases
    end

    # Find Zone
    zoneKeys =
        findall(n -> is_cgns_node(n) && (read(attributes(n)["label"]) == "Zone_t"), file)
    zones = map(hdf5_to_node, zoneKeys)
    if length(zones) == 1
        close(file)
        return first(zones)
    elseif length(zones) > 1
        close(file)
        return zones
    end

    close(file)

    error("Could not find any CGNSBase or Zone")
end

function is_cgns_node(n)
    attrs = attributes(n)
    haskey(attrs, "label") || (return false)
    return haskey(STRING_TO_CGNS_LABEL, read(attrs["label"]))
end

"""
    parse_tree_from_list(tree_as_list)

Construct a CGNS tree from a recursive list of nodes where each node is defined as
- the node name
- the node value
- a list of child nodes
- the node label
"""
function parse_tree_as_list(tree_as_list)
    @assert length(tree_as_list) == 4
    name = tree_as_list[1]
    value = tree_as_list[2]
    children = tree_as_list[3]
    label = tree_as_list[4]
    _children = if length(children) > 0
        map(parse_tree_as_list, children)
    else
        CGNSNode{Empty, Nothing, Nothing}[]
    end
    return CGNSNode(name, value, _children, label)
end

print_tree(tree) = print_tree(stdout, tree)

function print_tree(io, node::CGNSNode, indent::Int = 0)
    pad = "  "^indent
    colorSupported = get(io, :color, false)

    if colorSupported
        print(io, pad)
        printstyled(io, get_name(node); color = :blue)
        print(io, " ")
        printstyled(io, enum_to_string(get_label(node)); color = 0)
    else
        println(io, pad, get_name(node), " ", enum_to_string(get_label(node)))
    end
    if !isnothing(get_value(node))
        print(io, " ", get_value(node))
    end
    print(io, "\n")

    # if node.data !== nothing
    #     print(io, "  data=", summary(node.data))
    # end
    # println(io)

    for child in get_children(node)
        print_tree(io, child, indent + 1)
    end
end
