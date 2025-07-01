struct CGNSHDF5IoHandler <: Bcube.AbstractIoHandler end

Bcube._filename_to_handler(::Union{Val{:cgns}, Val{:hdf5}, Val{:hdf}}) = CGNSHDF5IoHandler()

function get_child(parent; name = "", type = "")
    child_name = findfirst(child -> child_match(child, name, type), parent)
    return isnothing(child_name) ? nothing : parent[child_name]
end

function get_children(parent; name = "", type = "")
    child_names = findall(child -> child_match(child, name, type), parent)
    return map(child_name -> parent[child_name], child_names)
end

get_name(obj) = String(last(split(HDF5.name(obj), '/')))

get_data_type(obj) = attributes(obj)["type"][]

get_cgns_type(obj) = haskey(attributes(obj), "label") ? attributes(obj)["label"][] : nothing

function get_value(obj)
    data_type = get_data_type(obj)
    data = read(obj[" data"])
    if data_type == "C1"
        return String(UInt8.(data))
    elseif data_type in ("I4", "I8", "R4", "R8")
        return data
    else
        error("Datatype '$(data_type)' not handled")
    end
end

function get_cgns_base(obj)
    cgnsBases = get_children(obj; type = "CGNSBase_t")
    if length(cgnsBases) == 0
        error("Could not find any CGNSBase_t node in the file")
    elseif length(cgnsBases) > 1
        error("The file contains several CGNSBase_t nodes, only one base is supported")
    end
    return first(cgnsBases)
end

has_child(parent; kwargs...) = !isnothing(get_child(parent; kwargs...))