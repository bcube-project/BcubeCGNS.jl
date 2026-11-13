# Independent Maia trees: x=(0,1,3), y=(0,2,5), z=(0,4), truncated
# to each dimension. CGNS uses i-fastest vertex/cell order. Hand-written
# connectivity and signed measures detect axis swaps and reflected cells.
function structured_fixture(d)
    ReaderCGNS.read_cgns_file(joinpath(assets_path, "structured-$(d)d.cgns"))
end
function structured_bc_tree(d, location, indices; range = false)
    edit_reader_tree(
        structured_fixture(d),
        n -> begin
            if reader_name(n) == "Left"
                children = ReaderCGNS.CGNSNode[reader_node(
                    range ? "PointRange" : "PointList",
                    indices,
                    [],
                    range ? "IndexRange_t" : "IndexArray_t",
                ),]
                isnothing(location) || push!(
                    children,
                    reader_node("GridLocation", location, [], "GridLocation_t"),
                )
                return reader_node("Selected", "BCWall", children, "BC_t")
            end
            n
        end,
    )
end
function structured_faces(tree, name = "Selected")
    mesh = BcubeCGNS.read_tree(tree).mesh
    f2n = Bcube.connectivities_indices(mesh, :f2n)
    sort([sort(collect(f2n[f])) for f in Bcube.boundary_faces(mesh, name)])
end
@testset "Structured topology and fields" begin
    # Node ID = i + 3(j-1) + 9(k-1); each layer has nine vertices.
    # Literal CCW lower corners and corresponding upper corners define positive cells.
    expected = (
        [[1, 2], [2, 3]],
        [[1, 2, 5, 4], [2, 3, 6, 5], [4, 5, 8, 7], [5, 6, 9, 8]],
        [
            [1, 2, 5, 4, 10, 11, 14, 13],
            [2, 3, 6, 5, 11, 12, 15, 14],
            [4, 5, 8, 7, 13, 14, 17, 16],
            [5, 6, 9, 8, 14, 15, 18, 17],
        ],
    )
    types = (Bcube.Bar2_t, Bcube.Quad4_t, Bcube.Hexa8_t)
    for d in 1:3
        result = read_file(joinpath(assets_path, "structured-$(d)d.cgns"); varnames = "*")
        mesh = result.mesh
        xy = [collect(Bcube.get_coords(n)) for n in Bcube.get_nodes(mesh)]
        axes = ([0.0, 1.0, 3.0], [0.0, 2.0, 5.0], [0.0, 4.0])[1:d]
        @test xy == vec([collect(p) for p in Iterators.product(axes...)])
        c2n = Bcube.connectivities_indices(mesh, :c2n)
        @test [collect(c2n[i]) for i in 1:ncells(mesh)] == expected[d]
        @test all(t -> t isa types[d], Bcube.cells(mesh))
        # Axis midpoints are (0.5,2), (1,3.5), (2); products follow cell i-fastest
        # order, independently of the reader's structured connectivity generation.
        centers = [sum(xy[j] for j in conn)/length(conn) for conn in expected[d]]
        mid = ([0.5, 2.0], [1.0, 3.5], [2.0])[1:d]
        @test centers == vec([collect(p) for p in Iterators.product(mid...)])
        # Orthogonal edge determinant is signed, and uses the tested Bcube
        # corner ordering (x edge 1->2, y edge 1->4, z edge 1->5).
        measures = map(expected[d]) do c
            u = xy[c[2]] - xy[c[1]]
            d == 1 && return u[1]
            v = xy[c[4]] - xy[c[1]]
            d == 2 && return u[1]*v[2] - u[2]*v[1]
            w = xy[c[5]] - xy[c[1]]
            u[1]*(v[2]*w[3]-v[3]*w[2]) - u[2]*(v[1]*w[3]-v[3]*w[1]) +
            u[3]*(v[1]*w[2]-v[2]*w[1])
        end
        @test measures == (
            d == 1 ? [1.0, 2.0] : d == 2 ? [2.0, 4.0, 3.0, 6.0] : [8.0, 16.0, 12.0, 24.0]
        )
        # Default absolute IDs and field markers follow the same linear order;
        # Cells/id is ten times the ID, distinct from nodal markers.
        @test collect(Bcube.get_absolute_node_indices(mesh)) == collect(1:length(xy))
        @test collect(Bcube.get_absolute_cell_indices(mesh)) ==
              collect(1:length(expected[d]))
        @test Bcube.get_values(result.data["Nodes"]["id"]) == collect(1.0:length(xy))
        @test Bcube.get_values(result.data["Cells"]["id"]) ==
              10 .* collect(1.0:length(expected[d]))
        # Left fixes i=1: endpoint in 1D, two j-edges in 2D, two j/k faces in 3D.
        @test structured_faces(structured_fixture(d), "Left") == (
            d == 1 ? [[1]] : d == 2 ? [[1, 4], [4, 7]] : [[1, 4, 10, 13], [4, 7, 13, 16]]
        )
    end
end

@testset "Structured BC indexing" begin
    # Vertex lists are columns of (i,j,k), not flattened node IDs. Default
    # Vertex, descending ranges and noncontiguous lists retain node semantics.
    @test structured_faces(structured_bc_tree(1, nothing, reshape([3], 1, 1))) == [[3]]
    @test structured_faces(structured_bc_tree(2, nothing, [1 1; 3 1]; range = true)) ==
          [[1, 4], [4, 7]]
    @test structured_faces(structured_bc_tree(3, "Vertex", [3 3 3 3; 1 2 1 2; 1 1 2 2])) ==
          [[3, 6, 12, 15]]
    # Directional locations index the minimum connecting vertex: normal index
    # reaches nv, tangential indices only nv-1. Disconnected patches must not
    # gain faces from the union of their vertices. Exercise both sides and all
    # axes, with ranges and explicitly reordered lists.
    cases = (
        ("IFaceCenter", [3 3; 2 1; 1 1], [[3, 6, 12, 15], [6, 9, 15, 18]]),
        ("JFaceCenter", [2 1; 3 3; 1 1], [[7, 8, 16, 17], [8, 9, 17, 18]]),
        ("KFaceCenter", [2 1; 2 1; 2 2], [[10, 11, 13, 14], [14, 15, 17, 18]]),
    )
    for (loc, ids, faces) in cases
        @test structured_faces(structured_bc_tree(3, loc, ids)) == faces
    end
    @test structured_faces(
        structured_bc_tree(3, "IFaceCenter", [1 1; 2 1; 1 1]; range = true),
    ) == [[1, 4, 10, 13], [4, 7, 13, 16]]
    @test structured_faces(
        structured_bc_tree(3, "JFaceCenter", [2 1; 1 1; 1 1]; range = true),
    ) == [[1, 2, 10, 11], [2, 3, 11, 12]]
    @test structured_faces(
        structured_bc_tree(3, "KFaceCenter", [1 2; 1 2; 1 1]; range = true),
    ) == [[1, 2, 4, 5], [2, 3, 5, 6], [4, 5, 7, 8], [5, 6, 8, 9]]
    # 2D directional aliases describe edges; standard EdgeCenter ranges infer
    # the unique constant exterior normal axis. Ambiguous lists are rejected.
    @test structured_faces(structured_bc_tree(2, "IFaceCenter", [3 3; 2 1])) ==
          [[3, 6], [6, 9]]
    @test structured_faces(
        structured_bc_tree(2, "JFaceCenter", [2 1; 3 3]; range = true),
    ) == [[7, 8], [8, 9]]
    @test structured_faces(structured_bc_tree(2, "EdgeCenter", [3 3; 1 2]; range = true)) ==
          [[3, 6], [6, 9]]
    # Malformed variants isolate bounds, tuple shape, exterior membership and
    # location rules; rejection must happen before constructing a misleading BC.
    for (d, loc, ids, range) in (
        (2, "Vertex", [0 1; 1 1], false), # indices are one-based
        (2, "Vertex", reshape([1, 2], 1, 2), false), # missing j row
        (3, "IFaceCenter", [1; 3; 1;;], false), # j exceeds tangential cell bound 2
        (3, "IFaceCenter", [2; 1; 1;;], false), # interior face
        (2, "KFaceCenter", [1; 1;;], false), # no K axis in 2D
        (3, "CellCenter", [1; 1; 1;;], false), # volume location is not a face BC
        (3, "FaceCenter", [1; 1; 1;;], false), # no direction
        (2, "EdgeCenter", [1; 1;;], false), # ambiguous corner
        (2, "Vertex", zeros(Int, 2, 0), false), # empty selection
        (2, "Vertex", ones(Int, 2, 3), true), # range needs two endpoints
    )
        @test_throws ErrorException BcubeCGNS.read_tree(
            structured_bc_tree(d, loc, ids; range),
        )
    end
    # Selecting both complete K planes includes every vertex, but must not
    # accidentally tag the eight side faces as a node-union implementation does.
    @test structured_faces(
        structured_bc_tree(3, "KFaceCenter", [1 2; 1 2; 1 2]; range = true),
    ) == [
        [1, 2, 4, 5],
        [2, 3, 5, 6],
        [4, 5, 7, 8],
        [5, 6, 8, 9],
        [10, 11, 13, 14],
        [11, 12, 14, 15],
        [13, 14, 16, 17],
        [14, 15, 17, 18],
    ]
    # Reordered right-side tuples must remain unchanged for later dataset mapping.
    ids = [3 3; 2 1; 1 1]
    tree = structured_bc_tree(3, "IFaceCenter", ids)
    before = copy(ids)
    structured_faces(tree)
    @test ids == before
    # Descending i runs fastest: (2,2), (1,2), then the j=1 row.
    # Inspect the resolver directly because sorted face sets cannot detect order.
    selection = reader_node("PointRange", [2 1; 2 1; 1 1], [], "IndexRange_t")
    @test BcubeCGNS.structured_bc_indices(selection, [3, 3, 2], "Selected") ==
          [(2, 2, 1), (1, 2, 1), (2, 1, 1), (1, 1, 1)]
    # Inconsistent dimensions must fail before coordinate or Mesh construction.
    tree = edit_reader_tree(
        structured_fixture(2),
        n -> if reader_name(n) == "Grid"
            reader_node("Grid", [3 1 0; 3 2 0], reader_children(n), "Zone_t")
        else
            n
        end,
    )
    @test_throws ErrorException BcubeCGNS.read_tree(tree)
end

# Family-based public names are preserved. These disjoint patches select
# bottom-left and lower-right exterior edges, not the edge connecting them.
# Resolving the union of their vertices first would invent bottom-right [2,3].
function shared_family_tree(structured, explicit; overlap = false)
    tree = structured ? structured_fixture(2) : reader_bc_tree("FaceCenter", [3, 5])
    edit_reader_tree(
        tree,
        n -> begin
            if ReaderCGNS.get_label(n) == ReaderCGNS.CGNSBase_t
                push!(reader_children(n), reader_node("Wall", nothing, [], "Family_t"))
            elseif ReaderCGNS.get_label(n) == ReaderCGNS.ZoneBC_t
                patches = ReaderCGNS.CGNSNode[]
                for (name, ids, loc) in (
                    (
                        "BottomPatch",
                        if structured
                            (explicit ? reshape([1, 1], 2, 1) : [1 2; 1 1])
                        else
                            (explicit ? [3] : [1, 2])
                        end,
                        if explicit
                            (structured ? "JFaceCenter" : "EdgeCenter")
                        else
                            "Vertex"
                        end,
                    ),
                    (
                        "RightPatch",
                        if structured
                            (explicit ? reshape([3, 1], 2, 1) : [3 3; 1 2])
                        else
                            (explicit ? [5] : [3, 6])
                        end,
                        if explicit
                            (structured ? "IFaceCenter" : "EdgeCenter")
                        else
                            "Vertex"
                        end,
                    ),
                )
                    push!(
                        patches,
                        reader_node(
                            name,
                            "BCWall",
                            [
                                reader_node("PointList", ids, [], "IndexArray_t"),
                                reader_node("GridLocation", loc, [], "GridLocation_t"),
                                reader_node("FamilyName", "Wall", [], "FamilyName_t"),
                            ],
                            "BC_t",
                        ),
                    )
                end
                if overlap
                    ids = structured ? [1 2; 1 1] : [1, 2]
                    push!(
                        patches,
                        reader_node(
                            "Other",
                            "BCWall",
                            [reader_node("PointList", ids, [], "IndexArray_t")],
                            "BC_t",
                        ),
                    )
                end
                return reader_node("ZoneBC", nothing, patches, "ZoneBC_t")
            end
            n
        end,
    )
end
@testset "Shared BC families" begin
    @testset "structured=$structured explicit=$explicit" for structured in (false, true),
        explicit in (false, true)

        mesh = BcubeCGNS.read_tree(shared_family_tree(structured, explicit)).mesh
        f2n = Bcube.connectivities_indices(mesh, :f2n)
        @test Set(Bcube.boundary_names(mesh)) == Set([:Wall])
        @test sort([sort(collect(f2n[f])) for f in Bcube.boundary_faces(mesh, "Wall")]) == [[1, 2], [3, 6]]
        @test sort(collect(mesh.bc_nodes.Wall)) == [1, 2, 3, 6]
        f2bc = Bcube.connectivities_indices(mesh, :f2bc)
        @test [f for f in 1:length(f2bc) if !isempty(f2bc[f])] == Bcube.boundary_faces(mesh, "Wall")
    end
    # Across different output names, preserve precedence in Bcube's named-tuple
    # order: the first output boundary owns a shared face, never both names.
    mesh = BcubeCGNS.read_tree(shared_family_tree(false, false; overlap = true)).mesh
    f2n = Bcube.connectivities_indices(mesh, :f2n)
    names = collect(keys(mesh.bc_nodes))
    owner = first(names)
    @test all(
        name ->
            ([1, 2] in [sort(collect(f2n[f])) for f in Bcube.boundary_faces(mesh, name)]) == (name == owner),
        names,
    )
end
