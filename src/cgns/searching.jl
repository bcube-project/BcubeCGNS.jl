"""
    get_nodes_from_predicate(
        root,
        condition;
        depth = typemax(Int),
        nmax = typemax(Int),
        kwargs...,
    )

Return nodes satisfying `condition`, including the root, in traversal order.

`depth` limits recursion and `nmax` limits the number of results.
"""
function get_nodes_from_predicate(
    root,
    condition;
    depth = typemax(Int),
    nmax = typemax(Int),
    kwargs...,
)
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

function get_nodes_from_label(root, label::Union{String, CGNSLabel}; kwargs...)
    L = label isa CGNSLabel ? label : STRING_TO_CGNS_LABEL[label]
    get_nodes_from_predicate(root, n -> n isa CGNSNode{L}; kwargs...)
end

function get_nodes_from_name(root, name::Union{String, Regex}; kwargs...)
    get_nodes_from_predicate(root, n -> matches_name(name, get_name(n)); kwargs...)
end

matches_name(name::AbstractString, candidate) = name == candidate
matches_name(name::Regex, candidate) = occursin(name, candidate)

function get_nodes_from_name_and_label(root, name, label; kwargs...)
    L = label isa CGNSLabel ? label : STRING_TO_CGNS_LABEL[label]
    get_nodes_from_predicate(
        root,
        n -> matches_name(name, get_name(n)) && (get_label(n) == L);
        kwargs...,
    )
end

function get_node_from_predicate(root, condition; kwargs...)
    first(get_nodes_from_predicate(root, condition; nmax = 1, kwargs...))
end

function get_node_from_label(root, label; kwargs...)
    first(get_nodes_from_label(root, label; nmax = 1, kwargs...))
end

function get_node_from_name_and_label(root, name, label; kwargs...)
    first(get_nodes_from_name_and_label(root, name, label; nmax = 1, kwargs...))
end

function get_children_from_predicate(root, condition; nmax = typemax(Int))
    result = CGNSNode[]
    for child in get_children(root)
        length(result) >= nmax && break
        condition(child) && push!(result, child)
    end
    return result
end

function get_children(root; name = nothing, label = nothing, kwargs...)
    (isnothing(name) && isnothing(label)) &&
        throw(ArgumentError("At least one of name and label is required"))
    if isnothing(name)
        return get_children_from_label(root, label; kwargs...)
    elseif isnothing(label)
        return get_children_from_name(root, name; kwargs...)
    else
        return get_children_from_name_and_label(root, name, label; kwargs...)
    end
end

function get_children_from_label(root, label; kwargs...)
    L = label isa CGNSLabel ? label : STRING_TO_CGNS_LABEL[label]
    get_children_from_predicate(root, n -> get_label(n) == L; kwargs...)
end

function get_children_from_name(root, name; kwargs...)
    get_children_from_predicate(root, n -> matches_name(name, get_name(n)); kwargs...)
end

function get_children_from_name_and_label(root, name, label; kwargs...)
    L = label isa CGNSLabel ? label : STRING_TO_CGNS_LABEL[label]
    get_children_from_predicate(
        root,
        n -> matches_name(name, get_name(n)) && get_label(n) == L;
        kwargs...,
    )
end

function get_child_from_predicate(root, condition; required = false, unique = false)
    nodes = get_children_from_predicate(root, condition; nmax = unique ? 2 : 1)
    required &&
        isempty(nodes) &&
        error("Missing child matching predicate in '$(get_name(root))'")
    unique &&
        length(nodes) > 1 &&
        error("Duplicate children matching predicate in '$(get_name(root))'")
    return first_or_nothing(nodes)
end

function get_child_from_label(root, label; kwargs...)
    get_child(root; label, kwargs...)
end

function get_child_from_name(root, name; kwargs...)
    get_child(root; name, kwargs...)
end

function get_child_from_name_and_label(root, name, label; kwargs...)
    get_child(root; name, label, kwargs...)
end

"""
    get_child(root; name = nothing, label = nothing, required = false, unique = false)

Return the first matching direct child, or `nothing` when absent.

Strings match exact names; Regex inputs request pattern matching. At least one
of `name` and `label` is required. Set `required` to reject missing matches and
`unique` to reject duplicates independently. The parent is never a candidate.
"""
function get_child(root; name = nothing, label = nothing, required = false, unique = false)
    nodes = get_children(root; name, label, nmax = unique ? 2 : 1)
    if (required && isempty(nodes)) || (unique && length(nodes) > 1)
        context = "child (name=$(repr(name)), label=$label) in '$(get_name(root))'"
        reason = isempty(nodes) ? "Missing" : "Duplicate"
        error("$reason $context")
    end
    return first_or_nothing(nodes)
end

function has_child(root; name = nothing, label = nothing, kwargs...)
    !isnothing(get_child(root; name, label, kwargs...))
end

function has_child(root, condition; kwargs...)
    !isnothing(get_child_from_predicate(root, condition; kwargs...))
end

first_or_nothing(x::Nothing) = x
first_or_nothing(x::AbstractVector) = length(x) > 0 ? first(x) : nothing
