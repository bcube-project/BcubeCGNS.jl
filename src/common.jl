const BCUBE_ENTITY_TO_CGNS_ENTITY = Dict(
    Bcube.Node_t => 2,
    Bcube.Bar2_t => 3,
    Bcube.Bar3_t => 4,
    Bcube.Tri3_t => 5,
    Bcube.Tri6_t => 6,
    Bcube.Tri9_t => 25,
    Bcube.Tri10_t => 26,
    Bcube.Quad4_t => 7,
    Bcube.Quad8_t => 8,
    Bcube.Quad9_t => 9,
    Bcube.Quad16_t => 27,
    Bcube.Tetra4_t => 10,
    Bcube.Tetra10_t => 11,
    Bcube.Pyra5_t => 12,
    Bcube.Penta6_t => 14,
    Bcube.Hexa8_t => 17,
    Bcube.Hexa27_t => 27,
)
const CGNS_ENTITY_TO_BCUBE_ENTITY = Dict(v => k for (k, v) in BCUBE_ENTITY_TO_CGNS_ENTITY)

function bcube_entity_to_cgns_entity(::T) where {T <: Bcube.AbstractEntityType}
    BCUBE_ENTITY_TO_CGNS_ENTITY[T]
end

cgns_entity_to_bcube_entity(code::Integer) = CGNS_ENTITY_TO_BCUBE_ENTITY[code]()

function child_match(child, name, type)
    if get_name(child) == name
        if length(name) > 0 && length(type) > 0
            (get_cgns_type(child) == type) && (return true)
        elseif length(name) > 0
            return true
        end
    end

    if get_cgns_type(child) == type
        if length(name) > 0 && length(type) > 0
            (get_name(child) == name) && (return true)
        elseif length(type) > 0
            return true
        end
    end

    return false
end