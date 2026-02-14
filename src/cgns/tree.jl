""" Top level abstract type for all abstract types defined here """
abstract type AbstractCGNSType end

abstract type AbstractCGNSNode <: AbstractCGNSType end

"""
    CGNSNode{L, V, C} <: AbstractCGNSNode

`L` is an CGNSLabel
`V` is usually a `Vector`
`C` is a vector of `CGNSNode` or `Nothing` if there is no child (I have not decided yet)
"""
struct CGNSNode{L, V, C} <: AbstractCGNSNode
    name::String
    value::V
    children::C
end

get_name(node::CGNSNode) = node.name
get_label(::CGNSNode{L}) where {L} = L
get_label_as_string(::CGNSNode{L}) where {L} = string(nameof(L))
get_children(node::CGNSNode) = node.children
get_value(node::CGNSNode) = node.value

function CGNSNode(name::String, value, children, label::String)
    _label = STRING_TO_CGNS_LABEL[label]
    return CGNSNode{_label, typeof(value), typeof(children)}(name, value, children)
end

# Note : I'm not sure the @enum is the good solution. Maybe an `AbstractCGNSLabel` and several `struct`
# would be better?
@enum CGNSLabel begin
    AdditionalExponents_t
    AdditionalFamilyName_t
    AdditionalUnits_t
    ArbitraryGridMotion_t
    Area_t
    AverageInterface_t
    Axisymmetry_t
    BaseIterativeData_t
    BC_t
    BCData_t
    BCDataSet_t
    BCProperty_t
    CGNSBase_t
    CGNSLibraryVersion_t
    CGNSTree_t
    ChemicalKineticsModel_t
    ConvergenceHistory_t
    DataArray_t
    DataConversion_t
    DataClass_t
    DimensionalExponents_t
    DimensionalUnits_t
    DiscreteData_t
    Descriptor_t
    Elements_t
    EMConductivityModel_t
    EMElectricFieldModel_t
    EMMagneticFieldModel_t
    Family_t
    FamilyBC_t
    FamilyBCDataSet_t
    FamilyName_t
    FlowEquationSet_t
    FlowSolution_t
    GasModel_t
    GeometryEntity_t
    GeometryFile_t
    GeometryFormat_t
    GeometryReference_t
    GoverningEquations_t
    Gravity_t
    GridConnectivity_t
    GridConnectivity1to1_t
    GridConnectivityProperty_t
    GridConnectivityType_t
    GridCoordinates_t
    GridLocation_t
    IndexArray_t
    IndexRange_t
    IntegralData_t
    ParticleBreakupModel_t
    ParticleCollisionModel_t
    ParticleCoordinates_t
    ParticleEquationSet_t
    ParticleForceModel_t
    ParticleGoverningEquations_t
    ParticlePhaseChangeModel_t
    ParticleSolution_t
    ParticleWallInteractionModel_t
    ParticleZone_t
    Periodic_t
    PointList_t
    Ordinal_t
    OversetHoles_t
    ReferenceState_t
    RigidGridMotion_t
    RotatingCoordinates_t
    Rind_t
    SimulationType_t
    ThermalConductivityModel_t
    ThermalRelaxationModel_t
    TurbulenceClosure_t
    TurbulenceModel_t
    UserDefinedData_t
    ViscosityModel_t
    WallFunction_t
    Zone_t
    ZoneBC_t
    ZoneGridConnectivity_t
    ZoneIterativeData_t
    ZoneSubRegion_t
    ZoneType_t
    Empty
end

const STRING_TO_CGNS_LABEL = Dict(
    "AdditionalExponents_t" => AdditionalExponents_t,
    "AdditionalFamilyName_t" => AdditionalFamilyName_t,
    "AdditionalUnits_t" => AdditionalUnits_t,
    "ArbitraryGridMotion_t" => ArbitraryGridMotion_t,
    "Area_t" => Area_t,
    "AverageInterface_t" => AverageInterface_t,
    "Axisymmetry_t" => Axisymmetry_t,
    "BaseIterativeData_t" => BaseIterativeData_t,
    "BC_t" => BC_t,
    "BCData_t" => BCData_t,
    "BCDataSet_t" => BCDataSet_t,
    "BCProperty_t" => BCProperty_t,
    "CGNSBase_t" => CGNSBase_t,
    "CGNSLibraryVersion_t" => CGNSLibraryVersion_t,
    "CGNSTree_t" => CGNSTree_t,
    "ChemicalKineticsModel_t" => ChemicalKineticsModel_t,
    "ConvergenceHistory_t" => ConvergenceHistory_t,
    "DataArray_t" => DataArray_t,
    "DataConversion_t" => DataConversion_t,
    "DataClass_t" => DataClass_t,
    "DimensionalExponents_t" => DimensionalExponents_t,
    "DimensionalUnits_t" => DimensionalUnits_t,
    "DiscreteData_t" => DiscreteData_t,
    "Descriptor_t" => Descriptor_t,
    "Elements_t" => Elements_t,
    "EMConductivityModel_t" => EMConductivityModel_t,
    "EMElectricFieldModel_t" => EMElectricFieldModel_t,
    "EMMagneticFieldModel_t" => EMMagneticFieldModel_t,
    "Family_t" => Family_t,
    "FamilyBC_t" => FamilyBC_t,
    "FamilyBCDataSet_t" => FamilyBCDataSet_t,
    "FamilyName_t" => FamilyName_t,
    "FlowEquationSet_t" => FlowEquationSet_t,
    "FlowSolution_t" => FlowSolution_t,
    "GasModel_t" => GasModel_t,
    "GeometryEntity_t" => GeometryEntity_t,
    "GeometryFile_t" => GeometryFile_t,
    "GeometryFormat_t" => GeometryFormat_t,
    "GeometryReference_t" => GeometryReference_t,
    "GoverningEquations_t" => GoverningEquations_t,
    "Gravity_t" => Gravity_t,
    "GridConnectivity_t" => GridConnectivity_t,
    "GridConnectivity1to1_t" => GridConnectivity1to1_t,
    "GridConnectivityProperty_t" => GridConnectivityProperty_t,
    "GridConnectivityType_t" => GridConnectivityType_t,
    "GridCoordinates_t" => GridCoordinates_t,
    "GridLocation_t" => GridLocation_t,
    "IndexArray_t" => IndexArray_t,
    "IndexRange_t" => IndexRange_t,
    "IntegralData_t" => IntegralData_t,
    "ParticleBreakupModel_t" => ParticleBreakupModel_t,
    "ParticleCollisionModel_t" => ParticleCollisionModel_t,
    "ParticleCoordinates_t" => ParticleCoordinates_t,
    "ParticleEquationSet_t" => ParticleEquationSet_t,
    "ParticleForceModel_t" => ParticleForceModel_t,
    "ParticleGoverningEquations_t" => ParticleGoverningEquations_t,
    "ParticlePhaseChangeModel_t" => ParticlePhaseChangeModel_t,
    "ParticleSolution_t" => ParticleSolution_t,
    "ParticleWallInteractionModel_t" => ParticleWallInteractionModel_t,
    "ParticleZone_t" => ParticleZone_t,
    "Periodic_t" => Periodic_t,
    "PointList_t" => PointList_t,
    "Ordinal_t" => Ordinal_t,
    "OversetHoles_t" => OversetHoles_t,
    "ReferenceState_t" => ReferenceState_t,
    "RigidGridMotion_t" => RigidGridMotion_t,
    "RotatingCoordinates_t" => RotatingCoordinates_t,
    "Rind_t" => Rind_t,
    "SimulationType_t" => SimulationType_t,
    "ThermalConductivityModel_t" => ThermalConductivityModel_t,
    "ThermalRelaxationModel_t" => ThermalRelaxationModel_t,
    "TurbulenceClosure_t" => TurbulenceClosure_t,
    "TurbulenceModel_t" => TurbulenceModel_t,
    "UserDefinedData_t" => UserDefinedData_t,
    "ViscosityModel_t" => ViscosityModel_t,
    "WallFunction_t" => WallFunction_t,
    "Zone_t" => Zone_t,
    "ZoneBC_t" => ZoneBC_t,
    "ZoneGridConnectivity_t" => ZoneGridConnectivity_t,
    "ZoneIterativeData_t" => ZoneIterativeData_t,
    "ZoneSubRegion_t" => ZoneSubRegion_t,
    "ZoneType_t" => ZoneType_t,
)