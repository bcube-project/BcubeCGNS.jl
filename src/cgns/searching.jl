"""
`nmax`: max number of results
"""
function get_nodes_from_predicate(root, condition; depth = typemax(Int), nmax = typemax(Int), kwargs...)
    result = CGNSNode[]
    get_nodes_from_predicate!(result, root, condition, depth, nmax)
    return result
end

function get_nodes_from_predicate!(result, root, condition, depth, nmax)
    condition(root) && push!(result, root)
    (length(result) < nmax) || return
    if depth > 0
        foreach(get_children(root)) do child
            get_nodes_from_predicate!(result, child, condition, depth - 1, nmax)
            (length(result) < nmax) || return
        end
    end
end

get_node_from_predicate(root, condition; kwargs...) = first(get_nodes_from_predicate(root, condition; nmax = 1, kwargs...))

function get_nodes_from_label(root, label::Union{String, CGNSLabel}; kwargs...)
    L = label isa CGNSLabel ? label : STRING_TO_CGNS_LABEL[label]
    get_nodes_from_predicate(root, n -> n isa CGNSNode{L}; kwargs...)
end

function get_nodes_from_name(root, name::Union{String, Regex}; kwargs...)
    r = name isa Regex ? name : Regex(name)
    get_nodes_from_predicate(root, n -> occursin(r, get_name(n)); kwargs...)
end

function get_children_from_predicate(root, condition; kwargs...)
    get_nodes_from_predicate(root, condition; depth = 1, kwargs...)
end

function get_children_from_label(root, label; kwargs...)
    get_nodes_from_label(root, label; depth = 1, kwargs...)
end

function get_children_from_name(root, name; kwargs...)
    get_nodes_from_name(root, name; depth = 1, kwargs...)
end

get_child_from_predicate(root, condition; kwargs...) = first(get_children_from_predicate(root, condition; nmax = 1))
get_child_from_label(root, label; kwargs...) = first(get_children_from_label(root, label; nmax = 1))
get_child_from_name(root, name; kwargs...) = first(get_children_from_name(root, name; nmax = 1))