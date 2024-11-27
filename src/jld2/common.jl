struct CGNSJLD2IoHandler <: Bcube.AbstractIoHandler end

Bcube._filename_to_handler(::Val{:jld2}) = CGNSJLD2IoHandler()

"""
Wrapper for the JLD2.Group struct to keep track of the attributes associated to a Group
without needing the parent.
"""
struct Node{T}
    node::T
    name::String
    attrs::Dict{Symbol, String}
end

get_wrapped_node(n::Node) = n.node
get_attribute(n::Node, name::Symbol) = haskey(n.attrs, name) ? n.attrs[name] : nothing
get_cgns_type(n::Node) = get_attribute(n, :label)
get_data_type(n::Node) = get_attribute(n, :type)
get_cgns_name(n::Node) = get_attribute(n, :name)
get_name(n::Node) = n.name

function Node(root::JLD2.JLDFile)
    attrs = parse_attributes(root, "")
    return Node{typeof(root)}(root, "root", attrs)
end

function Node(parent, node_name::String)
    attrs = parse_attributes(parent, node_name)
    node = parent[node_name]
    return Node{typeof(node)}(parent[node_name], node_name, attrs)
end

Node(parent::Node, node_name::String) = Node(get_wrapped_node(parent), node_name)

function parse_attributes(parent, node_name)
    raw_attrs = JLD2.load_attributes(parent, node_name)
    attrs = Dict{Symbol, String}()
    attrs_keys = (:name, :label, :type)
    for (key, val) in raw_attrs
        if key in attrs_keys
            attrs[key] = first(split(val, "\0"))
        end
    end
    return attrs
end

"""
We could make this function type-stable by converting the "type" attribute to a type-parameter of `Node`
"""
function jld2_get_value(n::Node)
    data_type = get_data_type(n)
    data = get_wrapped_node(n)[" data"]
    if data_type == "C1"
        return String(UInt8.(data))
    elseif data_type in ("I4", "I8", "R4", "R8")
        return data
    else
        error("Datatype '$(data_type)' not handled")
    end
end

function jld2_get_child(parent; name = "", type = "")
    for child_name in keys(get_wrapped_node(parent))
        child = Node(parent, child_name)
        child_match(child, name, type) && (return child)
    end
end

function jld2_get_children(parent; name = "", type = "")
    filtered = filter(
        child_name -> child_match(Node(parent, child_name), name, type),
        keys(get_wrapped_node(parent)),
    )
    map(child_name -> Node(parent, child_name), filtered)
end
