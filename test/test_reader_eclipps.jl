# Serialize literal trees directly, independent of the package mesh writer.
function write_literal_eclipps(filepath, tree)
    BcubeCGNS.h5open(filepath, "w") do file
        # Encode CGNS/HDF5 attributes explicitly; this is not a mesh-writer
        # roundtrip oracle. The literal connectivity above is the reference.
        attrs=BcubeCGNS.attributes(file)
        attrs["name"]="HDF5 MotherNode"
        attrs["label"]="Root Node of HDF5 File"
        function store(parent, node)
            group=BcubeCGNS.create_group(parent, reader_name(node))
            attrs=BcubeCGNS.attributes(group)
            attrs["name"]=reader_name(node)
            attrs["label"]=ReaderCGNS.enum_to_string(ReaderCGNS.get_label(node))
            value=reader_value(node)
            if !isnothing(value)
                attrs["type"]=if value isa String
                    "C1"
                elseif eltype(value) <: Integer
                    "I8"
                else
                    "R8"
                end
                group[" data"]=value isa String ? collect(codeunits(value)) : value
            end
            for child in reader_children(node)
                store(group, child)
            end
        end
        for base in reader_children(tree)
            store(file, base)
        end
    end
end

# Literal two-quad geometry, independent of both the package writer and Maia.
# Widths 1 and 2, height 2; surface IDs 3 and 4 are bottom edges oriented +x.
# The surface list reverses those IDs. Vertex order [3,2,1,6] retains an unused
# node, so connectivity, dataset ordering and compaction can be checked separately.
function eclipps_fixture(;
    structured = false,
    location = "Vertex",
    selection = nothing,
    datasets = nothing,
    surface = nothing,
)
    n = reader_node
    array(name, v) = n(name, v, [], "DataArray_t")
    if isnothing(selection)
        selection = if structured
            [3 2 1 3; 1 1 1 2]
        else
            reshape(location == "Vertex" ? [3, 2, 1, 6] : [3, 4], 1, :)
        end
    end
    if isnothing(surface)
        surface = structured ? [2, 1, 3, 1, 1, 1, 2, 1] : [4, 3]
    end
    if isnothing(datasets)
        vals = location == "Vertex" ? [30.0, 20.0, 10.0, 60.0] : [300.0, 400.0]
        datasets = [
            n(
                "IcingWallData:Gas",
                "BCWall",
                [
                    n("DirichletData", nothing, [array("Temperature", vals)], "BCData_t"),
                    n("NeumannData", nothing, [array("SkinFrictionX", [-2.0])], "BCData_t"),
                ],
                "BCDataSet_t",
            ),
        ]
    end
    bc = n(
        "Wall",
        "BCWall",
        [
            n("GridLocation", location, [], "GridLocation_t"),
            n("PointList", selection, [], "IndexArray_t"),
            n(
                "SurfacicElementList",
                nothing,
                [array("List", surface)],
                "UserDefinedData_t",
            ),
            datasets...,
        ],
        "BC_t",
    )
    grid = n(
        "GridCoordinates",
        nothing,
        [
            array("CoordinateX", [0.0, 1.0, 3.0, 0.0, 1.0, 3.0]),
            array("CoordinateY", [0.0, 0.0, 0.0, 2.0, 2.0, 2.0]),
            array("CoordinateZ", zeros(6)),
        ],
        "GridCoordinates_t",
    )
    sections = if structured
        []
    else
        [
            n(
                "Volume",
                Int32[7, 0],
                [
                    n("ElementRange", [1, 2], [], "IndexRange_t"),
                    array("ElementConnectivity", [1, 2, 5, 4, 2, 3, 6, 5]),
                ],
                "Elements_t",
            ),
            n(
                "Bottom",
                Int32[3, 0],
                [
                    n("ElementRange", [3, 4], [], "IndexRange_t"),
                    array("ElementConnectivity", [1, 2, 2, 3]),
                ],
                "Elements_t",
            ),
        ]
    end
    # A distractor BC precedes Wall; selection is by node name, never child index.
    other = n(
        "Other",
        "BCWall",
        [
            n(
                "PointList",
                structured ? reshape([1, 2], 2, 1) : reshape([4, 5], 1, :),
                [],
                "IndexArray_t",
            ),
        ],
        "BC_t",
    )
    zone = n(
        "Zone",
        structured ? [3 2 0; 2 1 0] : reshape([6, 2, 0], 1, :),
        [
            n("ZoneType", structured ? "Structured" : "Unstructured", [], "ZoneType_t"),
            grid,
            sections...,
            n("ZoneBC", nothing, [other, bc], "ZoneBC_t"),
        ],
        "Zone_t",
    )
    # Only one named reference quantity is supplied: this reader must not demand
    # the complete atmospheric/droplet input set listed for solvers in Table 1.
    base = n(
        "Base",
        [2, 3],
        [
            n(
                "ReferenceState",
                nothing,
                [array("VelocityMagnitudeReference", [42.0])],
                "ReferenceState_t",
            ),
            zone,
        ],
        "CGNSBase_t",
    )
    return n("CGNSTree", nothing, [base], "CGNSTree_t")
end

function eclipps_edit(tree, name, value; label = nothing)
    edit_reader_tree(
        tree,
        n -> if reader_name(n) == name
            reader_node(
                name,
                value,
                reader_children(n),
                if isnothing(label)
                    ReaderCGNS.enum_to_string(ReaderCGNS.get_label(n))
                else
                    label
                end,
            )
        else
            n
        end,
    )
end

@testset "Eclipps ordered extraction" begin
    for structured in (false, true)
        tree = eclipps_fixture(; structured)
        result = BcubeCGNS.extract_eclipps_tree(tree, ["Wall"])
        mesh = result.mesh
        # The returned node order is the BC order, including unused vertex 6.
        # Reversed surface ID order yields the long edge first, both oriented +x.
        xy = [collect(Bcube.get_coords(n)) for n in Bcube.get_nodes(mesh)]
        @test xy == [[3.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [3.0, 2.0, 0.0]]
        conn = Bcube.connectivities_indices(mesh, :c2n)
        @test [collect(conn[i]) for i in 1:2] == [[2, 1], [3, 2]]
        @test all(t -> t isa Bcube.Bar2_t, Bcube.cells(mesh))
        @test [xy[conn[i][2]][1]-xy[conn[i][1]][1] for i in 1:2] == [2.0, 1.0]
        @test [(xy[conn[i][1]]+xy[conn[i][2]])/2 for i in 1:2] == [[2.0, 0.0, 0.0], [0.5, 0.0, 0.0]]
        @test collect(Bcube.get_absolute_node_indices(mesh)) == 1:4
        @test collect(Bcube.get_absolute_cell_indices(mesh)) == 1:2
        gas = result.data["IcingWallData:Gas"]
        @test Bcube.get_values(gas["Temperature"]) == [30.0, 20.0, 10.0, 60.0]
        @test Bcube.get_location(gas["Temperature"]) isa Bcube.PointData
        @test Bcube.get_values(gas["SkinFrictionX"]) == fill(-2.0, 4)
        @test result.refState["VelocityMagnitudeReference"] == 42.0

        # Change to descending PointRange: all three bottom vertices, no unused
        # node. Values now have exactly three entries, in descending i order.
        ranged = edit_reader_tree(
            tree,
            n -> begin
                if reader_name(n) == "PointList" && size(reader_value(n), 2) == 4
                    v = structured ? [3 1; 1 1] : reshape([3, 1], 1, :)
                    return reader_node("PointRange", v, [], "IndexRange_t")
                elseif reader_name(n) == "Temperature"
                    return reader_node("Temperature", [30.0, 20.0, 10.0], [], "DataArray_t")
                end
                n
            end,
        )
        r = BcubeCGNS.extract_eclipps_tree(ranged, "Wall")
        @test Bcube.nnodes(r.mesh) == 3
        @test Bcube.get_values(r.data["IcingWallData:Gas"]["Temperature"]) ==
              [30.0, 20.0, 10.0]
    end

    # Face selection [3,4] differs from surface order [4,3]. Temperature must
    # follow surface cells, while referenced vertices compact in ascending order.
    tree = eclipps_fixture(; location = "FaceCenter")
    r = BcubeCGNS.extract_eclipps_tree(tree, "Wall")
    @test Bcube.nnodes(r.mesh) == 3
    @test [collect(Bcube.connectivities_indices(r.mesh, :c2n)[i]) for i in 1:2] == [[2, 3], [1, 2]]
    @test Bcube.get_values(r.data["IcingWallData:Gas"]["Temperature"]) == [400.0, 300.0]
    @test Bcube.get_location(r.data["IcingWallData:Gas"]["Temperature"]) isa Bcube.CellData
    # Downstream Bcube/VTK methods require typed mesh-data dictionaries;
    # Dict{String,Any} contents alone cannot satisfy that dispatch contract.
    @test valtype(typeof(r.data["IcingWallData:Gas"])) <: Bcube.MeshData
    @test valtype(typeof(r.data)) <: AbstractDict

    # A full explicit dataset permutation overrides inherited point ordering.
    # A second FaceCenter dataset on the Vertex BC uses its own element IDs.
    n=reader_node
    nodal=n(
        "Nodes",
        "BCWall",
        [
            n("PointList", reshape([1, 6, 2, 3], 1, :), [], "IndexArray_t"),
            n(
                "DirichletData",
                nothing,
                [n("Pressure", [10.0, 60.0, 20.0, 30.0], [], "DataArray_t")],
                "BCData_t",
            ),
        ],
        "BCDataSet_t",
    )
    face=n(
        "Faces",
        "BCWall",
        [
            n("GridLocation", "FaceCenter", [], "GridLocation_t"),
            n("PointRange", [3, 4], [], "IndexRange_t"),
            n(
                "NeumannData",
                nothing,
                [n("MassFlux", [3.0, 4.0], [], "DataArray_t")],
                "BCData_t",
            ),
        ],
        "BCDataSet_t",
    )
    r=BcubeCGNS.extract_eclipps_tree(eclipps_fixture(; datasets = [nodal, face]), "Wall")
    @test Bcube.get_values(r.data["Nodes"]["Pressure"]) == [30.0, 20.0, 10.0, 60.0]
    @test Bcube.get_values(r.data["Faces"]["MassFlux"]) == [4.0, 3.0]
    @test Bcube.get_location(r.data["Faces"]["MassFlux"]) isa Bcube.CellData

    # Three vertex values and three surface cells make length inference
    # ambiguous. Replace the volume with the triangle (1,2,4) and use its
    # complete CCW boundary; an explicit FaceCenter dataset stays cell-associated.
    equal = eclipps_fixture(;
        location = "FaceCenter",
        selection = reshape([3, 4, 5], 1, :),
        surface = [5, 4, 3],
        datasets = [
            n(
                "Equal",
                "BCWall",
                [
                    n(
                        "DirichletData",
                        nothing,
                        [n("Pressure", [3.0, 4.0, 5.0], [], "DataArray_t")],
                        "BCData_t",
                    ),
                ],
                "BCDataSet_t",
            ),
        ],
    )
    equal=edit_reader_tree(
        equal,
        n -> begin
            if reader_name(n)=="Bottom"
                return reader_node(
                    "Bottom",
                    Int32[3, 0],
                    [
                        reader_node("ElementRange", [3, 5], [], "IndexRange_t"),
                        reader_node(
                            "ElementConnectivity",
                            [1, 2, 2, 4, 4, 1],
                            [],
                            "DataArray_t",
                        ),
                    ],
                    "Elements_t",
                )
            elseif reader_name(n)=="Volume"
                return reader_node(
                    "Volume",
                    Int32[5, 0],
                    [
                        reader_node("ElementRange", [1, 1], [], "IndexRange_t"),
                        reader_node("ElementConnectivity", [1, 2, 4], [], "DataArray_t"),
                    ],
                    "Elements_t",
                )
            elseif reader_name(n)=="Zone"
                return reader_node(
                    "Zone",
                    reshape([6, 1, 0], 1, :),
                    reader_children(n),
                    "Zone_t",
                )
            end
            n
        end,
    )
    r=BcubeCGNS.extract_eclipps_tree(equal, "Wall")
    @test Bcube.nnodes(r.mesh)==ncells(r.mesh)==3
    @test Bcube.get_location(r.data["Equal"]["Pressure"]) isa Bcube.CellData
    @test Bcube.get_values(r.data["Equal"]["Pressure"]) == [5.0, 4.0, 3.0]
end

@testset "Eclipps diagnostics" begin
    tree=eclipps_fixture()
    # Selection cardinality and unknown names get stable API errors.
    @test_throws ArgumentError BcubeCGNS.extract_eclipps_tree(tree, ["Wall", "Other"])
    @test_throws ArgumentError BcubeCGNS.extract_eclipps_tree(tree, "Missing")
    # Each mutation isolates an invalid list, missing vertex, unsupported
    # location or wrong local array length. No solver fields are required.
    for bad in (
        eclipps_edit(tree, "List", [4, 99]), # missing surface element
        eclipps_edit(tree, "List", [4, 4]), # duplicate element
        eclipps_edit(tree, "List", [1]), # volume ID is not a surface ID
        eclipps_edit(tree, "PointList", reshape([3, 2, 6], 1, :)), # omits vertex 1
        eclipps_edit(tree, "GridLocation", "CellCenter"),
        eclipps_edit(tree, "Temperature", [1.0, 2.0]),
        eclipps_edit(tree, "PointList", reshape([3, 2, 1, 0], 1, :)), # invalid zero
        eclipps_edit(eclipps_fixture(; structured = true), "List", [1, 1, 2]), # incomplete corner tuples
        eclipps_edit(
            eclipps_fixture(; structured = true),
            "List",
            [2, 1, 4, 1, 1, 1, 2, 1],
        ), # i=4 out of bounds
    )
        @test_throws ErrorException BcubeCGNS.extract_eclipps_tree(bad, "Wall")
    end
    # One missing equation group and empty datasets are valid; duplicate names
    # across groups cannot be represented in the public flat variable dictionary.
    no_neumann=edit_reader_tree(
        tree,
        n -> if reader_name(n)=="IcingWallData:Gas"
            reader_node(
                reader_name(n),
                reader_value(n),
                filter(c -> reader_name(c)!="NeumannData", reader_children(n)),
                "BCDataSet_t",
            )
        else
            n
        end,
    )
    @test Set(
        keys(BcubeCGNS.extract_eclipps_tree(no_neumann, "Wall").data["IcingWallData:Gas"]),
    ) == Set(["Temperature"])
    # Neumann-only and no datasets are also useful partial exports; neither
    # implies a missing required solver-input error in a general reader.
    no_dirichlet=edit_reader_tree(
        tree,
        n -> if reader_name(n)=="IcingWallData:Gas"
            reader_node(
                reader_name(n),
                reader_value(n),
                filter(c -> reader_name(c)!="DirichletData", reader_children(n)),
                "BCDataSet_t",
            )
        else
            n
        end,
    )
    @test Set(
        keys(
            BcubeCGNS.extract_eclipps_tree(no_dirichlet, "Wall").data["IcingWallData:Gas"],
        ),
    ) == Set(["SkinFrictionX"])
    @test isempty(
        BcubeCGNS.extract_eclipps_tree(eclipps_fixture(; datasets = []), "Wall").data,
    )
    duplicate=edit_reader_tree(
        tree,
        n -> if reader_name(n)=="SkinFrictionX"
            reader_node("Temperature", [-2.0], [], "DataArray_t")
        else
            n
        end,
    )
    @test_throws ErrorException BcubeCGNS.extract_eclipps_tree(duplicate, "Wall")
    # A dataset cannot change location without saying which faces its values
    # refer to, nor provide a partial or foreign selection for a full mesh array.
    for children in (
        [reader_node("GridLocation", "FaceCenter", [], "GridLocation_t")],
        [reader_node("GridLocation", "EdgeCenter", [], "GridLocation_t")],
        [reader_node("PointList", reshape([1, 2], 1, :), [], "IndexArray_t")],
    )
        dataset=reader_node("Invalid", "BCWall", children, "BCDataSet_t")
        @test_throws ErrorException BcubeCGNS.extract_eclipps_tree(
            eclipps_fixture(; datasets = [dataset]),
            "Wall",
        )
    end
    # Embedded unsupported MIXED records still abort even when the selected
    # surface lies entirely in a different section and never references it.
    mixed=reader_node(
        "UnusedMixed",
        Int32[20, 0],
        [
            reader_node("ElementRange", [8, 8], [], "IndexRange_t"),
            reader_node("ElementConnectivity", [23, 1, 2, 3], [], "DataArray_t"),
        ],
        "Elements_t",
    )
    bad=edit_reader_tree(
        tree,
        n -> if reader_name(n)=="Zone"
            reader_node("Zone", reader_value(n), [reader_children(n)..., mixed], "Zone_t")
        else
            n
        end,
    )
    @test_throws ErrorException BcubeCGNS.extract_eclipps_tree(bad, "Wall")
end

@testset "Eclipps structured faces and public file API" begin
    n=reader_node
    # Face indices are the minimum corner tuples: (1,1), (2,1). The
    # surface list reverses that order, so local values must be reversed too.
    tree=eclipps_fixture(;
        structured = true,
        location = "FaceCenter",
        selection = [1 2; 1 1],
    )
    r=BcubeCGNS.extract_eclipps_tree(tree, "Wall")
    @test Bcube.get_values(r.data["IcingWallData:Gas"]["Temperature"]) == [400.0, 300.0]
    @test Bcube.get_location(r.data["IcingWallData:Gas"]["Temperature"]) isa Bcube.CellData
    @test [collect(Bcube.get_coords(x)) for x in Bcube.get_nodes(r.mesh)] == [[1.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    # A dataset can independently reorder those same face tuples; association
    # remains cell-based even on a Vertex BC with an extra unused node.
    dataset=n(
        "Faces",
        "BCWall",
        [
            n("GridLocation", "FaceCenter", [], "GridLocation_t"),
            n("PointList", [1 2; 1 1], [], "IndexArray_t"),
            n(
                "DirichletData",
                nothing,
                [n("Pressure", [10.0, 20.0], [], "DataArray_t")],
                "BCData_t",
            ),
        ],
        "BCDataSet_t",
    )
    r=BcubeCGNS.extract_eclipps_tree(
        eclipps_fixture(; structured = true, datasets = [dataset]),
        "Wall",
    )
    @test Bcube.get_values(r.data["Faces"]["Pressure"]) == [20.0, 10.0]
    # Generic HDF5 node serialization of the literal tree is independent of the
    # package's mesh writer. This exercises file opening and singleton selectors.
    mktempdir() do directory
        filepath=joinpath(directory, "literal-eclipps.cgns")
        write_literal_eclipps(filepath, tree)
        r=BcubeCGNS.extract_surf_from_eclipps(filepath, "Wall")
        @test ncells(r.mesh)==2
        @test Bcube.get_values(r.data["IcingWallData:Gas"]["Temperature"]) == [400.0, 300.0]
        @test ncells(BcubeCGNS.extract_surf_from_eclipps(filepath, ["Wall"]).mesh)==2
        @test_throws ArgumentError BcubeCGNS.extract_surf_from_eclipps(
            filepath,
            "Wall";
            timestep = 1,
        )
    end
end

@testset "Eclipps 3D corner tuples and orientation" begin
    n=reader_node
    for structured in (true, false)
        # One 1 x 2 x 3 hexahedron. The z=0 quad has the deliberately stored
        # +z orientation (not the volume's outward normal). Extraction preserves
        # that loop rather than silently flipping a solver-provided surface.
        grid=n(
            "GridCoordinates",
            nothing,
            [
                n(
                    "CoordinateX",
                    [0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0],
                    [],
                    "DataArray_t",
                ),
                n(
                    "CoordinateY",
                    [0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 2.0, 2.0],
                    [],
                    "DataArray_t",
                ),
                n(
                    "CoordinateZ",
                    [0.0, 0.0, 0.0, 0.0, 3.0, 3.0, 3.0, 3.0],
                    [],
                    "DataArray_t",
                ),
            ],
            "GridCoordinates_t",
        )
        points=structured ? [1 2 2 1; 2 2 1 1; 1 1 1 1] : reshape([3, 4, 2, 1], 1, :)
        corners=structured ? [1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1] : [2]
        bc=n(
            "Wall",
            "BCWall",
            [
                n("PointList", points, [], "IndexArray_t"),
                n(
                    "SurfacicElementList",
                    nothing,
                    [n("List", corners, [], "DataArray_t")],
                    "UserDefinedData_t",
                ),
                n(
                    "Gas",
                    "BCWall",
                    [
                        n(
                            "DirichletData",
                            nothing,
                            [n("Pressure", [30.0, 40.0, 20.0, 10.0], [], "DataArray_t")],
                            "BCData_t",
                        ),
                    ],
                    "BCDataSet_t",
                ),
            ],
            "BC_t",
        )
        elements=if structured
            []
        else
            [
                n(
                    "Volume",
                    [17, 0],
                    [
                        n("ElementRange", [1, 1], [], "IndexRange_t"),
                        n(
                            "ElementConnectivity",
                            [1, 2, 4, 3, 5, 6, 8, 7],
                            [],
                            "DataArray_t",
                        ),
                    ],
                    "Elements_t",
                ),
                n(
                    "Surface",
                    [7, 0],
                    [
                        n("ElementRange", [2, 2], [], "IndexRange_t"),
                        n("ElementConnectivity", [1, 2, 4, 3], [], "DataArray_t"),
                    ],
                    "Elements_t",
                ),
            ]
        end
        zone=n(
            "Zone",
            structured ? [2 1 0; 2 1 0; 2 1 0] : reshape([8, 1, 0], 1, :),
            [
                n("ZoneType", structured ? "Structured" : "Unstructured", [], "ZoneType_t"),
                grid,
                elements...,
                n("ZoneBC", nothing, [bc], "ZoneBC_t"),
            ],
            "Zone_t",
        )
        tree=n("CGNSTree", nothing, [n("Base", [3, 3], [zone], "CGNSBase_t")], "CGNSTree_t")
        r=BcubeCGNS.extract_eclipps_tree(tree, "Wall")
        # Omitted BC location defaults to Vertex. The returned quad follows
        # the original corner cycle through the reversed local node numbering.
        @test only(Bcube.cells(r.mesh)) isa Bcube.Quad4_t
        conn=collect(Bcube.connectivities_indices(r.mesh, :c2n)[1])
        @test conn == [4, 3, 2, 1]
        xy=[collect(Bcube.get_coords(x)) for x in Bcube.get_nodes(r.mesh)]
        @test xy == [[0.0, 2.0, 0.0], [1.0, 2.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
        p=xy[conn]
        @test sum(p[j][1]*p[mod1(j+1, 4)][2]-p[mod1(j+1, 4)][1]*p[j][2] for j in 1:4)/2 ==
              2.0
        @test sum(p)/4 == [0.5, 1.0, 0.0]
        @test Bcube.get_values(r.data["Gas"]["Pressure"]) == [30.0, 40.0, 20.0, 10.0]
        @test isnothing(r.refState)
    end
end

@testset "Eclipps dataset remapping and stored grids" begin
    n=reader_node
    # Vertex values on a FaceCenter BC need an explicit vertex selection;
    # its order is [3,1,2], whereas extracted nodes are sorted [1,2,3].
    dataset=n(
        "Nodes",
        "BCWall",
        [
            n("GridLocation", "Vertex", [], "GridLocation_t"),
            n("PointList", reshape([3, 1, 2], 1, :), [], "IndexArray_t"),
            n(
                "DirichletData",
                nothing,
                [n("Temperature", [30.0, 10.0, 20.0], [], "DataArray_t")],
                "BCData_t",
            ),
        ],
        "BCDataSet_t",
    )
    r=BcubeCGNS.extract_eclipps_tree(
        eclipps_fixture(; location = "FaceCenter", datasets = [dataset]),
        "Wall",
    )
    @test Bcube.get_values(r.data["Nodes"]["Temperature"]) == [10.0, 20.0, 30.0]
    @test Bcube.get_location(r.data["Nodes"]["Temperature"]) isa Bcube.PointData
    # Structured explicit vertex tuples exercise the same permutation without
    # treating the three-dimensional embedding as a third index direction.
    dataset=n(
        "Nodes",
        "BCWall",
        [
            n("PointList", [1 3 2 3; 1 2 1 1], [], "IndexArray_t"),
            n(
                "DirichletData",
                nothing,
                [n("Temperature", [10.0, 60.0, 20.0, 30.0], [], "DataArray_t")],
                "BCData_t",
            ),
        ],
        "BCDataSet_t",
    )
    r=BcubeCGNS.extract_eclipps_tree(
        eclipps_fixture(; structured = true, datasets = [dataset]),
        "Wall",
    )
    @test Bcube.get_values(r.data["Nodes"]["Temperature"]) == [30.0, 20.0, 10.0, 60.0]
    # Repeated Vertex entries are an observed legacy encoding. Keep all data
    # slots and attach connectivity to the final occurrence, as the old API did.
    tree=eclipps_fixture(; selection = reshape([3, 2, 1, 1], 1, :))
    r=BcubeCGNS.extract_eclipps_tree(tree, "Wall")
    @test Bcube.nnodes(r.mesh)==4
    @test collect(Bcube.connectivities_indices(r.mesh, :c2n)[2]) == [4, 2]
    @test Bcube.get_values(r.data["IcingWallData:Gas"]["Temperature"]) ==
          [30.0, 20.0, 10.0, 60.0]

    # Two stored grids have opposite child and pointer order. Time values are
    # deliberately decreasing: step one must use the translated grid, not the
    # first grid child or the earliest numeric time. BC datasets remain static.
    timed=edit_reader_tree(
        eclipps_fixture(),
        node -> begin
            if reader_name(node)=="Base"
                metadata=n(
                    "Steps",
                    [2],
                    [n("TimeValues", [9.0, 1.0], [], "DataArray_t")],
                    "BaseIterativeData_t",
                )
                return n(
                    "Base",
                    reader_value(node),
                    [reader_children(node)..., metadata],
                    "CGNSBase_t",
                )
            elseif reader_name(node)=="Zone"
                grid=only(
                    filter(x -> reader_name(x)=="GridCoordinates", reader_children(node)),
                )
                shifted=edit_reader_tree(
                    grid,
                    x -> if reader_name(x)=="CoordinateX"
                        n("CoordinateX", reader_value(x) .+ 100, [], "DataArray_t")
                    else
                        x
                    end,
                )
                shifted=n(
                    "Shifted",
                    nothing,
                    reader_children(shifted),
                    "GridCoordinates_t",
                )
                metadata=n(
                    "Pointers",
                    nothing,
                    [
                        n(
                            "GridCoordinatesPointers",
                            ["Shifted", "GridCoordinates"],
                            [],
                            "DataArray_t",
                        ),
                    ],
                    "ZoneIterativeData_t",
                )
                return n(
                    "Zone",
                    reader_value(node),
                    [reader_children(node)..., shifted, metadata],
                    "Zone_t",
                )
            end
            node
        end,
    )
    @test_throws ArgumentError BcubeCGNS.extract_eclipps_tree(timed, "Wall")
    a=BcubeCGNS.extract_eclipps_tree(timed, "Wall"; timestep = 1)
    b=BcubeCGNS.extract_eclipps_tree(timed, "Wall"; timestep = 2)
    @test Bcube.get_coords(first(Bcube.get_nodes(a.mesh)))[1]==103.0
    @test Bcube.get_coords(first(Bcube.get_nodes(b.mesh)))[1]==3.0
    @test Bcube.get_values(a.data["IcingWallData:Gas"]["Temperature"]) ==
          Bcube.get_values(b.data["IcingWallData:Gas"]["Temperature"])
end

# Logical 3x3 (or 3x3x3) zones expose skipped-cell and interior selections.
# All vertices appear in the BC in i-fastest order so returned connectivity
# equals the hand-computed linear IDs, without a remapping oracle.
function structured_surface_tree(d, surface)
    dims = fill(3, d)
    points = collect(Iterators.product((1:3 for _ in 1:d)...)) |> vec
    xyz = ([0.0, 1.0, 3.0], [0.0, 2.0, 5.0], [0.0, 3.0, 7.0])
    edit_reader_tree(
        eclipps_fixture(; structured = true, datasets = []),
        n -> begin
            label = ReaderCGNS.get_label(n)
            if label == ReaderCGNS.CGNSBase_t
                return reader_node(reader_name(n), [d, 3], reader_children(n), "CGNSBase_t")
            elseif label == ReaderCGNS.Zone_t
                return reader_node(
                    reader_name(n),
                    hcat(dims, dims .- 1, zeros(Int, d)),
                    reader_children(n),
                    "Zone_t",
                )
            elseif reader_name(n) in ("CoordinateX", "CoordinateY", "CoordinateZ")
                k = findfirst(
                    ==(reader_name(n)),
                    ("CoordinateX", "CoordinateY", "CoordinateZ"),
                )
                values =
                    k <= d ? [xyz[k][p[k]] for p in points] : zeros(length(points))
                return reader_node(reader_name(n), values, [], "DataArray_t")
            elseif reader_name(n) == "PointList"
                return reader_node("PointList", hcat(collect.(points)...), [], "IndexArray_t")
            elseif reader_name(n) == "List"
                return reader_node("List", surface, [], "DataArray_t")
            end
            n
        end,
    )
end
@testset "Structured Eclipps elementary topology" begin
    # Invalid 2D: diagonal, skipped vertex, interior horizontal edge.
    # Invalid 3D: crossed cycle, width-two face, and interior z=2 face.
    for (d, corners) in (
        (2, [1, 1, 3, 2]),
        (2, [1, 1, 3, 1]),
        (2, [1, 2, 2, 2]),
        (3, [1, 1, 1, 2, 2, 1, 2, 1, 1, 1, 2, 1]),
        (3, [1, 1, 1, 3, 1, 1, 3, 2, 1, 1, 2, 1]),
        (3, [1, 1, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2]),
    )
        @test_throws ErrorException BcubeCGNS.extract_eclipps_tree(
            structured_surface_tree(d, corners),
            "Wall",
        )
    end
    # Both orientations are valid; neither corner order nor signed orientation
    # may be normalized by topology validation. Coordinates are embedded in 3D.
    for (d, corners, expected, measure) in (
        (2, [1, 1, 2, 1], [1, 2], 1.0),
        (2, [2, 1, 1, 1], [2, 1], -1.0),
        (3, [1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1], [1, 2, 5, 4], 2.0),
        (3, [1, 1, 1, 1, 2, 1, 2, 2, 1, 2, 1, 1], [1, 4, 5, 2], -2.0),
    )
        mesh =
            BcubeCGNS.extract_eclipps_tree(structured_surface_tree(d, corners), "Wall").mesh
        conn = collect(Bcube.connectivities_indices(mesh, :c2n)[1])
        @test conn == expected
        p = [collect(Bcube.get_coords(Bcube.get_nodes(mesh)[i])) for i in conn]
        signed = if d == 2
            p[2][1] - p[1][1]
        else
            sum(p[j][1]*p[mod1(j+1, 4)][2] - p[mod1(j+1, 4)][1]*p[j][2] for j in 1:4)/2
        end
        @test signed == measure
    end
    # Exercise the malformed diagonal through the public file API as well.
    mktempdir() do dir
        path = joinpath(dir, "diagonal.cgns")
        write_literal_eclipps(path, structured_surface_tree(2, [1, 1, 3, 2]))
        @test_throws ErrorException BcubeCGNS.extract_surf_from_eclipps(path, "Wall")
    end
end
