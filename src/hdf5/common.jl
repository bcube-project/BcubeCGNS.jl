struct CGNSHDF5IoHandler <: Bcube.AbstractIoHandler end

Bcube._filename_to_handler(::Union{Val{:cgns}, Val{:hdf5}, Val{:hdf}}) = CGNSHDF5IoHandler()

function get_cgns_base(obj)
    cgnsBases = get_children_from_label(obj, CGNS.CGNSBase_t)
    if length(cgnsBases) == 0
        error("Could not find any CGNSBase_t node in the file")
    elseif length(cgnsBases) > 1
        error("The file contains several CGNSBase_t nodes, only one base is supported")
    end
    return first(cgnsBases)
end
