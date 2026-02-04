abstract type AbstractCGNSType end

abstract type AbstractCGNSTree <: AbstractCGNSType end

"""
`bases` is a Tuple-or-Array of `AbstractCGNSBase`
"""
struct CGNSTree <: AbstractCGNSTree
    bases::B
end

""" https://cgns.org/standard/SIDS/hierarchy.html#cgns-entry-level-structure-definition-cgnsbase-t """
abstract type AbstractCGNSBase <: AbstractCGNSType end

"""
`zones` is a Tuple-or-Array of `AbstractZone`
"""
struct CGNSBase{Z} <: AbstractCGNSBase
    zones::Z
end

""" https://cgns.org/standard/SIDS/hierarchy.html#zone-structure-definition-zone-t """
abstract type AbstractZone <: AbstractCGNSType end

"""
`gridCoordinates` is a Tuple of Array: (X, Y[, Z])
`elements` is a Tuple-or-Array of `AbstractElements`
"""
struct UnstructuredZone{GD, E, Z} <: AbstractZone
    gridCoordinates::GD
    elements::E
    zoneBC::ZoneBC
end

""" https://cgns.org/standard/SIDS/grid.html#elements-structure-definition-elements-t """
abstract type AbstractElements <: AbstractCGNSType end

"""
`E` is the Bcube element type (`Bar2_t`, `Tri3_t` etc)
`connectivity` is a flatten vector of the element-to-nodes connectivity
`range` represents the standard ElementRange node.
"""
struct CanonicalElements{E, C, R} <: AbstractElements
    connectivity::C
    range::R
end

""" https://cgns.org/standard/SIDS/boundary.html#zonal-boundary-condition-structure-definition-zonebc-t """
abstract type AbstractZoneBC <: AbstractCGNSType end

struct ZoneBC <: AbstractZoneBC end
