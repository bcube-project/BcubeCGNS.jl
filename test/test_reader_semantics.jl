# Independent Maia fixture: two CCW quads of widths 1 and 2, height 2.
# Node and cell IDs deliberately differ, exposing field association/order errors.
# Exact connectivity and signed shoelace areas detect permutations/reflections
# that counts, centers, or a writer roundtrip alone would miss.
@testset "Reader semantics: two quads" begin
    result = read_file(joinpath(assets_path, "two-quads.cgns"); varnames = "*")
    mesh = result.mesh
    # File node order runs along X first, then the upper row; unequal widths
    # make a coordinate permutation observable independently of cell counts.
    xy = [collect(Bcube.get_coords(n)) for n in Bcube.get_nodes(mesh)]
    @test xy == [[0.0, 0.0], [1.0, 0.0], [3.0, 0.0], [0.0, 2.0], [1.0, 2.0], [3.0, 2.0]]
    @test Bcube.nnodes(mesh) == 6
    @test ncells(mesh) == 2
    @test all(t -> t isa Bcube.Quad4_t, Bcube.cells(mesh))
    # CCW corner cycles join the middle edge 2--5; reflection changes signed
    # areas even though the arithmetic centers would remain unchanged.
    c2n = Bcube.connectivities_indices(mesh, :c2n)
    @test [collect(c2n[i]) for i in 1:2] == [[1, 2, 5, 4], [2, 3, 6, 5]]
    centers = [sum(xy[j] for j in c2n[i]) / 4 for i in 1:2]
    @test centers == [[0.5, 1.0], [2.0, 1.0]]
    areas = map(1:2) do i
        p = xy[c2n[i]]
        sum(p[j][1] * p[mod1(j+1, 4)][2] - p[mod1(j+1, 4)][1] * p[j][2] for j in 1:4) / 2
    end
    @test areas == [2.0, 4.0]
    # Without Maia numbering metadata, original vertex and element IDs survive.
    @test collect(Bcube.get_absolute_node_indices(mesh)) == collect(1:6)
    @test collect(Bcube.get_absolute_cell_indices(mesh)) == [1, 2]
    # Left selects vertices 1 and 4, so exactly that exterior edge is tagged.
    f2n = Bcube.connectivities_indices(mesh, :f2n)
    @test [sort(collect(f2n[f])) for f in Bcube.boundary_faces(mesh, "Left")] == [[1, 4]]
    # Markers 1:6 belong to vertices, while 10/20 follow the two cells; wrappers
    # and values together expose accidental association by array length.
    nodal = result.data["Nodes"]["node_id"]
    cell = result.data["Cells"]["cell_id"]
    @test Bcube.get_location(nodal) isa Bcube.PointData
    @test Bcube.get_location(cell) isa Bcube.CellData
    @test Bcube.get_values(nodal) == collect(1.0:6.0)
    @test Bcube.get_values(cell) == [10.0, 20.0]
    # Filtering retains the empty Nodes solution; the default skips all data.
    selected = read_file(joinpath(assets_path, "two-quads.cgns"); varnames = ["cell_id"])
    @test isempty(selected.data["Nodes"])
    @test Set(keys(selected.data["Cells"])) == Set(["cell_id"])
    @test isnothing(read_file(joinpath(assets_path, "two-quads.cgns")).data)
end

# Reuse the independent file's tree for malformed/default encodings. Rebuilding
# nodes keeps each case isolated; the real-file public API checks remain above.
const ReaderCGNS = BcubeCGNS.CGNS
function reader_node(name, value, children, label)
    ReaderCGNS.CGNSNode(name, value, ReaderCGNS.CGNSNode[children...], label)
end
function edit_reader_tree(node, edit)
    children = ReaderCGNS.CGNSNode[
        edit_reader_tree(c, edit) for c in ReaderCGNS.get_children(node)
    ]
    copy = reader_node(
        ReaderCGNS.get_name(node),
        deepcopy(ReaderCGNS.get_value(node)),
        children,
        ReaderCGNS.enum_to_string(ReaderCGNS.get_label(node)),
    )
    return edit(copy)
end
reader_fixture() = ReaderCGNS.read_cgns_file(joinpath(assets_path, "two-quads.cgns"))
reader_name(n) = ReaderCGNS.get_name(n)
reader_children(n) = ReaderCGNS.get_children(n)
reader_value(n) = ReaderCGNS.get_value(n)

@testset "Reader validation and CGNS defaults" begin
    # Coordinate names must be exact: a metadata array containing Y is not Y.
    tree = edit_reader_tree(
        reader_fixture(),
        n -> begin
            if reader_name(n) == "GridCoordinates"
                push!(reader_children(n), reader_node("MetadataY", [99.0], [], "DataArray_t"))
            end
            n
        end,
    )
    @test Bcube.nnodes(BcubeCGNS.read_tree(tree).mesh) == 6
    # Renaming Y removes a required axis; shortening it violates ZoneVertexCount.
    tree = edit_reader_tree(
        reader_fixture(),
        n -> begin
            if reader_name(n) == "CoordinateY"
                return reader_node("OtherY", reader_value(n), [], "DataArray_t")
            end
            n
        end,
    )
    @test_throws ErrorException BcubeCGNS.read_tree(tree)
    tree = edit_reader_tree(
        reader_fixture(),
        n -> begin
            if reader_name(n) == "CoordinateY"
                reader_node("CoordinateY", [0.0], [], "DataArray_t")
            else
                n
            end
        end,
    )
    @test_throws ErrorException BcubeCGNS.read_tree(tree)

    # Omitted GridLocation means Vertex, even if a wrong-length array happens
    # to match the cell count. Inferring CellCenter hides malformed CGNS data.
    for sol in ("Nodes", "Cells")
        tree = edit_reader_tree(
            reader_fixture(),
            n -> begin
                reader_name(n) == sol && filter!(
                    c -> reader_name(c) != "GridLocation",
                    reader_children(n),
                )
                n
            end,
        )
        if sol == "Nodes"
            result = BcubeCGNS.read_tree(tree; varnames = "*")
            @test Bcube.get_location(result.data["Nodes"]["node_id"]) isa Bcube.PointData
        else
            @test_throws ErrorException BcubeCGNS.read_tree(tree; varnames = "*")
        end
    end
    tree = edit_reader_tree(
        reader_fixture(),
        n -> begin
            if (reader_name(n) == "GridLocation" && reader_value(n) == "CellCenter")
                reader_node("GridLocation", "FaceCenter", [], "GridLocation_t")
            else
                n
            end
        end,
    )
    @test_throws ErrorException BcubeCGNS.read_tree(tree; varnames = "*")

    # Malformed fixed-arity connectivity must be rejected before Mesh building.
    # First variant truncates two QUAD_4 records; second references node 7
    # outside the fixture's six vertices despite having the correct arity.
    for conn in ([1, 2, 5], [1, 2, 5, 4, 2, 3, 7, 5])
        tree = edit_reader_tree(
            reader_fixture(),
            n -> begin
                if reader_name(n) == "ElementConnectivity"
                    reader_node(reader_name(n), conn, [], "DataArray_t")
                else
                    n
                end
            end,
        )
        @test_throws ErrorException BcubeCGNS.read_tree(tree)
    end

    # Present sections in descending CGNS ID order; cell fields remain stored
    # in CGNS cell order. Sorting geometry and IDs together preserves alignment.
    tree = edit_reader_tree(
        reader_fixture(),
        n -> begin
            if reader_name(n) == "TwoQuads"
                filter!(c -> reader_name(c) != "Quads", reader_children(n))
                for (name, id, conn) in
                    (("Second", 2, [2, 3, 6, 5]), ("First", 1, [1, 2, 5, 4]))
                    push!(
                        reader_children(n),
                        reader_node(
                            name,
                            [7, 0],
                            [
                                reader_node("ElementRange", [id, id], [], "IndexRange_t"),
                                reader_node("ElementConnectivity", conn, [], "DataArray_t"),
                            ],
                            "Elements_t",
                        ),
                    )
                end
            end
            n
        end,
    )
    mesh = BcubeCGNS.read_tree(tree; varnames = "*").mesh
    @test collect(Bcube.get_absolute_cell_indices(mesh)) == [1, 2]
    c2n = Bcube.connectivities_indices(mesh, :c2n)
    @test [collect(c2n[i]) for i in 1:2] == [[1, 2, 5, 4], [2, 3, 6, 5]]
end

# Reuse the hand-checked two-quad geometry. Boundary IDs 3:5 are deliberately
# split into three descending-order sections. Selecting bottom-left, top-left
# and right edges includes all six vertices: node union would wrongly add the
# left and bottom-right/top-right edges. No package writer is used as an oracle.
function reader_bc_tree(location, ids; range = false, gap = false)
    edit_reader_tree(
        reader_fixture(),
        n -> begin
            if reader_name(n) == "TwoQuads"
                filter!(
                    c -> ReaderCGNS.get_label(c) != ReaderCGNS.ZoneBC_t,
                    reader_children(n),
                )
                for (id, conn) in ((5, [3, 6]), (4, [4, 5]), (3, [1, 2]))
                    gap && id == 4 && continue
                    push!(
                        reader_children(n),
                        reader_node(
                            "Edge$id",
                            [3, 0],
                            [
                                reader_node("ElementRange", [id, id], [], "IndexRange_t"),
                                reader_node("ElementConnectivity", conn, [], "DataArray_t"),
                            ],
                            "Elements_t",
                        ),
                    )
                end
                children = ReaderCGNS.CGNSNode[reader_node(
                    range ? "PointRange" : "PointList",
                    ids,
                    [],
                    range ? "IndexRange_t" : "IndexArray_t",
                ),]
                isnothing(location) || push!(
                    children,
                    reader_node("GridLocation", location, [], "GridLocation_t"),
                )
                push!(
                    reader_children(n),
                    reader_node(
                        "ZoneBC",
                        nothing,
                        [reader_node("Selected", "BCWall", children, "BC_t")],
                        "ZoneBC_t",
                    ),
                )
            end
            n
        end,
    )
end
function reader_selected_faces(tree)
    mesh = BcubeCGNS.read_tree(tree).mesh
    f2n = Bcube.connectivities_indices(mesh, :f2n)
    sort([sort(collect(f2n[f])) for f in Bcube.boundary_faces(mesh, "Selected")])
end
@testset "Unstructured BC selections" begin
    # IDs 3/4/5 are bottom-left/top-left/right; child and selection order
    # must not affect membership. Ranges traverse all three separate sections.
    expected = [[1, 2], [3, 6], [4, 5]]
    for location in ("FaceCenter", "CellCenter", "EdgeCenter")
        @test reader_selected_faces(reader_bc_tree(location, [5, 3, 4])) == expected
        @test reader_selected_faces(reader_bc_tree(location, [3, 5]; range = true)) ==
              expected
        @test reader_selected_faces(reader_bc_tree(location, [5, 3]; range = true)) ==
              expected
    end
    @test reader_selected_faces(reader_bc_tree("FaceCenter", [5, 3])) == [[1, 2], [3, 6]]
    # Vertex selections use node sets: reversed left/bottom pairs select one
    # edge each, with the same result when GridLocation is omitted.
    for location in (nothing, "Vertex")
        @test reader_selected_faces(reader_bc_tree(location, [4, 1])) == [[1, 4]]
        @test reader_selected_faces(reader_bc_tree(location, [2, 1]; range = true)) ==
              [[1, 2]]
    end
    # Gaps, out-of-range IDs, mixed dimensions and volume elements at a face
    # location must fail clearly, before they can silently erase/add boundaries.
    for (location, ids, range, gap) in (
        ("FaceCenter", [3, 5], true, true),
        ("FaceCenter", [9], false, false),
        ("Vertex", [0, 1], false, false),
        ("Vertex", [7], false, false),
        ("FaceCenter", [1, 3], false, false), # quad plus edge: mixed dimensions
        ("FaceCenter", [1], false, false), # volume-only selection must not disappear
    )
        @test_throws ErrorException BcubeCGNS.read_tree(
            reader_bc_tree(location, ids; range, gap),
        )
    end
    # BC resolution must never reorder the underlying PointList (dataset order).
    tree = reader_bc_tree("FaceCenter", [5, 3])
    before = deepcopy(tree)
    BcubeCGNS.read_tree(tree)
    lists(t) = [
        reader_value(c) for z in ReaderCGNS.get_children(BcubeCGNS.get_cgns_base(t)) for
        b in reader_children(z) if reader_name(b) == "ZoneBC" for
        bc in reader_children(b) for
        c in reader_children(bc) if reader_name(c) == "PointList"
    ]
    @test lists(tree) == lists(before) == [[5, 3]]
end

@testset "Explicit BC mesh associations" begin
    tree = reader_bc_tree("FaceCenter", [5, 3])
    mesh = BcubeCGNS.read_tree(tree).mesh
    # Numeric face tags must identify exactly the same faces as the named BC.
    f2bc = Bcube.connectivities_indices(mesh, :f2bc)
    @test [f for f in 1:length(f2bc) if !isempty(f2bc[f])] == Bcube.boundary_faces(mesh, "Selected")
    # A section may describe an interior edge, but it cannot be an exterior BC.
    tree = edit_reader_tree(
        tree,
        n -> begin
            if reader_name(n) == "Edge5"
                return reader_node(
                    "Edge5",
                    [3, 0],
                    [
                        reader_node("ElementRange", [5, 5], [], "IndexRange_t"),
                        reader_node("ElementConnectivity", [2, 5], [], "DataArray_t"),
                    ],
                    "Elements_t",
                )
            end
            n
        end,
    )
    @test_throws ErrorException BcubeCGNS.read_tree(tree)
end

# Several BCs expose the distinction between CGNS child order and Bcube's
# dictionary iteration order. Both named face lists and numeric tags must agree.
@testset "Multiple BC tag ordering" begin
    tree = edit_reader_tree(
        reader_bc_tree("FaceCenter", [5, 3]),
        n -> begin
            if reader_name(n) == "ZoneBC"
                push!(
                    reader_children(n),
                    reader_node(
                        "Top",
                        "BCWall",
                        [
                            reader_node("GridLocation", "FaceCenter", [], "GridLocation_t"),
                            reader_node("PointList", [4], [], "IndexArray_t"),
                        ],
                        "BC_t",
                    ),
                )
            end
            n
        end,
    )
    mesh = BcubeCGNS.read_tree(tree).mesh
    f2n = Bcube.connectivities_indices(mesh, :f2n)
    f2bc = Bcube.connectivities_indices(mesh, :f2bc)
    expected = Dict(:Selected => [[1, 2], [3, 6]], :Top => [[4, 5]])
    @test all(enumerate(keys(Bcube.boundary_faces(mesh)))) do (tag, name)
        faces = Bcube.boundary_faces(mesh, name)
        sort([sort(collect(f2n[f])) for f in faces]) == expected[name] && all(f -> collect(f2bc[f]) == [tag], faces)
    end
end

@testset "Dimensional field association" begin
    # Two volume quads (1:2) and two boundary bars (3:4) deliberately have
    # equal counts. Volume markers must never become edge data by coincidence.
    tree = reader_bc_tree("FaceCenter", [3, 4]; gap = false)
    tree = edit_reader_tree(
        tree,
        n -> begin
            filter!(c -> reader_name(c) != "Edge5", reader_children(n))
            n
        end,
    )
    @test_throws ErrorException BcubeCGNS.read_tree(
        tree;
        topodim = 1,
        varnames = ["cell_id"],
    )
    @test ncells(BcubeCGNS.read_tree(tree; topodim = 1).mesh) == 2
    result = BcubeCGNS.read_tree(tree; topodim = 1, varnames = ["node_id"])
    @test Bcube.get_values(result.data["Nodes"]["node_id"]) == collect(1.0:6.0)
    @test isempty(result.data["Cells"])
    @test Bcube.get_values(
        BcubeCGNS.read_tree(tree; varnames = ["cell_id"]).data["Cells"]["cell_id"],
    ) == [10.0, 20.0]
end

@testset "Exact required children and empty boundaries" begin
    # Legal section names may contain reserved child names. A similarly named
    # direct child must not stand in for the required exact-name child either.
    @testset "Section $token" for token in (
        "ElementRange",
        "ElementConnectivity",
        "ElementStartOffset",
    )
        tree = edit_reader_tree(
            reader_fixture(),
            n -> begin
                if ReaderCGNS.get_label(n) == ReaderCGNS.Elements_t
                    pushfirst!(
                        reader_children(n),
                        reader_node("OtherElementRange", [99, 99], [], "IndexRange_t"),
                    )
                    pushfirst!(
                        reader_children(n),
                        reader_node("OtherElementConnectivity", [99], [], "DataArray_t"),
                    )
                    return reader_node(
                        token * "Quads",
                        reader_value(n),
                        reader_children(n),
                        "Elements_t",
                    )
                end
                n
            end,
        )
        result = BcubeCGNS.read_tree(tree; varnames = "*")
        @test [collect(c) for c in Bcube.connectivities_indices(result.mesh, :c2n)] == [[1, 2, 5, 4], [2, 3, 6, 5]]
        @test collect(Bcube.get_absolute_cell_indices(result.mesh)) == [1, 2]
        @test Bcube.get_values(result.data["Cells"]["cell_id"]) == [10.0, 20.0]
    end
    # Empty and absent ZoneBC carry identical meaning, including dictionary
    # types required by Bcube's mesh constructor.
    @testset "Empty ZoneBC $empty" for empty in (false, true)
        tree = edit_reader_tree(
            reader_fixture(),
            n -> begin
                filter!(
                    c -> ReaderCGNS.get_label(c) != ReaderCGNS.ZoneBC_t,
                    reader_children(n),
                )
                if empty && ReaderCGNS.get_label(n) == ReaderCGNS.Zone_t
                    push!(reader_children(n), reader_node("ZoneBC", nothing, [], "ZoneBC_t"))
                end
                n
            end,
        )
        mesh = BcubeCGNS.read_tree(tree).mesh
        @test ncells(mesh) == 2
        @test isempty(Bcube.boundary_names(mesh))
    end
end

@testset "EdgeCenter exact tags" begin
    # Standard 2D edge selections use original element IDs, just like legacy
    # aliases; no union of vertices may invent the unselected bottom-right edge.
    mesh = BcubeCGNS.read_tree(reader_bc_tree("EdgeCenter", [5, 3])).mesh
    f2n = Bcube.connectivities_indices(mesh, :f2n)
    @test sort([sort(collect(f2n[f])) for f in Bcube.boundary_faces(mesh, "Selected")]) == [[1, 2], [3, 6]]
    f2bc = Bcube.connectivities_indices(mesh, :f2bc)
    @test [f for f in 1:length(f2bc) if !isempty(f2bc[f])] == sort(collect(Bcube.boundary_faces(mesh, "Selected")))
end

@testset "Boundary node reporting order" begin
    # The left edge is selected in reverse file order. Membership is unchanged,
    # but reporting must preserve first occurrence instead of sorting node IDs.
    mesh = BcubeCGNS.read_tree(reader_bc_tree("Vertex", [4, 1, 4])).mesh
    @test Bcube.boundary_nodes(mesh, :Selected) == [4, 1]
    @test reader_selected_faces(reader_bc_tree("Vertex", [4, 1, 4])) == [[1, 4]]
end
