# Maia fixture: the two-quad rectangle, with the right quad split along 2->6.
# Section child order is Later (IDs 3:5), First (IDs 1:2). First contains
# QUAD_4/TRI_3; Later contains TRI_3/BAR_2/BAR_2. Cell fields follow IDs 1:3,
# not child order, and the boundary edges must not become cells. Expectations
# are hand-computed; neither the package writer nor Maia is needed at test time.
mixed_fixture() = ReaderCGNS.read_cgns_file(joinpath(assets_path, "mixed-quads-tris.cgns"))

@testset "Exact required MIXED offsets" begin
    # A legal section name containing the offset identifier must not shadow
    # its child; misleading array names must not supply the offsets either.
    tree = edit_reader_tree(
        mixed_fixture(),
        n -> begin
            if reader_name(n) == "Later"
                pushfirst!(
                    reader_children(n),
                    reader_node("OtherElementStartOffset", [99], [], "DataArray_t"),
                )
                return reader_node(
                    "ElementStartOffsetLater",
                    reader_value(n),
                    reader_children(n),
                    "Elements_t",
                )
            end
            n
        end,
    )
    result = BcubeCGNS.read_tree(tree; varnames = "*")
    @test collect(Bcube.get_absolute_cell_indices(result.mesh)) == [1, 2, 3]
    @test Bcube.get_values(result.data["Cells"]["cell_id"]) == [10.0, 20.0, 30.0]
end
function mixed_edit_section(edit; section = "Later")
    edit_reader_tree(mixed_fixture(), n -> reader_name(n) == section ? edit(n) : n)
end
# Replace just the encoded records; keep geometry, fields and BCs reusable.
function mixed_records(conn; offsets = nothing)
    mixed_edit_section(
        n -> begin
            filter!(
                c -> reader_name(c) ∉ ("ElementConnectivity", "ElementStartOffset"),
                reader_children(n),
            )
            push!(
                reader_children(n),
                reader_node("ElementConnectivity", conn, [], "DataArray_t"),
            )
            isnothing(offsets) || push!(
                reader_children(n),
                reader_node("ElementStartOffset", offsets, [], "DataArray_t"),
            )
            n
        end,
    )
end
@testset "MIXED topology, fields and boundaries" begin
    result = read_file(joinpath(assets_path, "mixed-quads-tris.cgns"); varnames = "*")
    mesh = result.mesh
    xy = [collect(Bcube.get_coords(n)) for n in Bcube.get_nodes(mesh)]
    @test xy == [[0.0, 0.0], [1.0, 0.0], [3.0, 0.0], [0.0, 2.0], [1.0, 2.0], [3.0, 2.0]]
    @test typeof.(Bcube.cells(mesh)) == [Bcube.Quad4_t, Bcube.Tri3_t, Bcube.Tri3_t]
    c2n = Bcube.connectivities_indices(mesh, :c2n)
    # The left quad has area 2; splitting the width-2 right quad gives two
    # CCW triangles of area 2 each. Their vertex means give the listed centers.
    @test [collect(c2n[i]) for i in 1:ncells(mesh)] == [[1, 2, 5, 4], [2, 3, 6], [2, 6, 5]]
    @test [sum(xy[j] for j in c2n[i])/length(c2n[i]) for i in 1:3] ≈ [[0.5, 1.0], [7/3, 2/3], [5/3, 4/3]]
    @test [
        sum(
            xy[c[j]][1]*xy[c[mod1(j+1, length(c))]][2]-xy[c[mod1(j+1, length(c))]][1]*xy[c[j]][2]
            for j in eachindex(c)
        )/2 for c in c2n
    ] == [2.0, 2.0, 2.0]
    # Compact cell markers follow original IDs 1:3, skipping boundary IDs 4:5;
    # reversed section order must not permute geometry or field association.
    @test collect(Bcube.get_absolute_cell_indices(mesh)) == [1, 2, 3]
    @test collect(Bcube.get_absolute_node_indices(mesh)) == collect(1:6)
    @test Bcube.get_values(result.data["Cells"]["cell_id"]) == [10.0, 20.0, 30.0]
    @test Bcube.get_location(result.data["Cells"]["cell_id"]) isa Bcube.CellData
    @test Bcube.get_values(result.data["Nodes"]["node_id"]) == collect(1.0:6.0)
    @test structured_faces(mixed_fixture(), "Ends") == [[1, 4], [3, 6]]
    # Legacy MIXED encodings omit offsets and must decode the same record stream.
    legacy = edit_reader_tree(
        mixed_fixture(),
        n -> begin
            filter!(c -> reader_name(c) != "ElementStartOffset", reader_children(n))
            n
        end,
    )
    @test structured_faces(legacy, "Ends") == [[1, 4], [3, 6]]
    @test Bcube.get_values(
        BcubeCGNS.read_tree(legacy; varnames = "*").data["Cells"]["cell_id"],
    ) == [10.0, 20.0, 30.0]
    # A range resolves both boundary records after the triangle of different arity.
    tree = edit_reader_tree(
        mixed_fixture(),
        n -> if reader_name(n) == "PointList"
            reader_node("PointRange", [5, 4], [], "IndexRange_t")
        else
            n
        end,
    )
    @test structured_faces(tree, "Ends") == [[1, 4], [3, 6]]
end
@testset "MIXED complete-section validation" begin
    # Replace the final edge by unsupported embedded types: PYRA_14, BAR_4,
    # nested MIXED, NGON, or an unknown code. The fixed-arity cases have complete
    # records outside the requested cell dimension. Remove BCs as well: every
    # record must be checked even when no cell or boundary selects it.
    for (code, arity, modern) in (
        (13, 14, false),
        (24, 4, false),
        (20, 0, false),
        (22, 3, false),
        (999, 0, false),
        (24, 4, true),
    )
        tree = mixed_records(
            vcat([5, 2, 6, 5, 3, 1, 4], code, ones(Int, arity));
            offsets = modern ? [0, 4, 7, 8+arity] : nothing,
        )
        tree = edit_reader_tree(
            tree,
            n -> begin
                filter!(
                    c -> ReaderCGNS.get_label(c) != ReaderCGNS.ZoneBC_t,
                    reader_children(n),
                )
                n
            end,
        )
        err = try
            BcubeCGNS.read_tree(tree)
            nothing
        catch e
            e
        end
        @test err isa ErrorException &&
              occursin("Later", sprint(showerror, err)) &&
              occursin("element 5", sprint(showerror, err)) &&
              occursin("$code", sprint(showerror, err))
    end
    # Modern offsets must agree exactly with the type-prefixed records. Detect
    # nonzero origins, nonmonotone starts, wrong arity/count and trailing data.
    for offsets in ([1, 5, 8, 11], [0, 4, 3, 10], [0, 4, 7, 9], [0, 5, 7, 10], [0, 4, 7])
        tree = mixed_edit_section(
            n -> begin
                for (i, c) in enumerate(reader_children(n))
                    reader_name(c) == "ElementStartOffset" && (
                        reader_children(n)[i] = reader_node(
                            "ElementStartOffset",
                            offsets,
                            [],
                            "DataArray_t",
                        )
                    )
                end
                n
            end,
        )
        @test_throws ErrorException BcubeCGNS.read_tree(tree)
    end
    # Distinct stream failures: truncated triangle, missing final edge,
    # trailing integer, and a complete edge referencing nonexistent vertex 7.
    for conn in (
        [5, 2, 6],
        [5, 2, 6, 5, 3, 1, 4],
        [5, 2, 6, 5, 3, 1, 4, 3, 3, 6, 1],
        [5, 2, 6, 5, 3, 1, 4, 3, 3, 7],
    )
        tree = mixed_records(conn)
        @test_throws ErrorException BcubeCGNS.read_tree(tree)
    end
    # Swapping IDs 3 (triangle) and 4 (edge) leaves compact cell fields in
    # sorted cell order, with original IDs 1,2,4. Boundary-first/interspersed
    # numbering is valid reader compatibility behavior (also used by NACA).
    tree = mixed_records([3, 1, 4, 5, 2, 6, 5, 3, 3, 6]; offsets = [0, 3, 7, 10])
    tree = edit_reader_tree(
        tree,
        n -> if reader_name(n) == "PointList"
            reader_node("PointList", [5, 3], [], "IndexArray_t")
        else
            n
        end,
    )
    result = BcubeCGNS.read_tree(tree; varnames = "*")
    @test collect(Bcube.get_absolute_cell_indices(result.mesh)) == [1, 2, 4]
    @test Bcube.get_values(result.data["Cells"]["cell_id"]) == [10.0, 20.0, 30.0]
    @test structured_faces(tree, "Ends") == [[1, 4], [3, 6]]
    # A BC selecting a triangle and an edge from one section is not a face patch.
    tree = edit_reader_tree(
        mixed_fixture(),
        n -> if reader_name(n) == "PointList"
            reader_node("PointList", [3, 4], [], "IndexArray_t")
        else
            n
        end,
    )
    @test_throws ErrorException BcubeCGNS.read_tree(tree)
end
