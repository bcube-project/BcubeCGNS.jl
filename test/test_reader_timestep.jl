# Independent Maia fixture reuses the two-quad topology. Stored step 1 points
# to Zulu (x+10, cell_id+100), step 2 to Alpha (original geometry/data).
# Children are in the opposite order and TimeValues=[9,2]: neither sorting
# times/names nor taking the first grid/merging solutions is correct.
timestep_fixture() = ReaderCGNS.read_cgns_file(joinpath(assets_path, "timesteps.cgns"))
function timestep_edit(name, edit)
    edit_reader_tree(timestep_fixture(), n -> reader_name(n) == name ? edit(n) : n)
end
function timestep_remove(name)
    edit_reader_tree(
        timestep_fixture(),
        n -> begin
            filter!(c -> reader_name(c) != name, reader_children(n))
            n
        end,
    )
end
function timestep_value(name, value)
    timestep_edit(
        name,
        n -> reader_node(
            name,
            value,
            reader_children(n),
            ReaderCGNS.enum_to_string(ReaderCGNS.get_label(n)),
        ),
    )
end
@testset "Stored timestep selection" begin
    file = joinpath(assets_path, "timesteps.cgns")
    # Explicit selection is required for both mesh-only and field reads.
    @test_throws ArgumentError read_file(file)
    @test_throws ArgumentError read_file(file; varnames = "*")
    for k in 1:2
        # Stored pointer order pairs x+10 with marker+100, then original data;
        # topology and element IDs remain fixed across these stored grids.
        result = read_file(file; timestep = k, varnames = "*")
        name, offset = k == 1 ? ("Zulu", 10.0) : ("Alpha", 0.0)
        xy = [collect(Bcube.get_coords(n)) for n in Bcube.get_nodes(result.mesh)]
        @test xy == [
            [x+offset, y] for (x, y) in
            ((0.0, 0.0), (1.0, 0.0), (3.0, 0.0), (0.0, 2.0), (1.0, 2.0), (3.0, 2.0))
        ]
        @test Set(keys(result.data)) == Set(["Solution" * name])
        @test Bcube.get_values(result.data["Solution" * name]["cell_id"]) ==
              [10.0, 20.0] .+ 10offset
        @test [collect(c) for c in Bcube.connectivities_indices(result.mesh, :c2n)] == [[1, 2, 5, 4], [2, 3, 6, 5]]
        @test collect(Bcube.get_absolute_cell_indices(result.mesh)) == [1, 2]
        @test isnothing(read_file(file; timestep = k).data)
    end
    # Bounds and integer type are strict; Bool must not act as integer step 1.
    for index in (0, -1, 3, 1.0, true, "1")
        @test_throws ArgumentError read_file(file; timestep = index)
    end
    @test_throws ArgumentError BcubeCGNS.read_tree(reader_fixture(); timestep = 1)
    # Metadata is found by label, not by the conventional node name.
    for name in
        ("StoredSteps", "StepMap", "GridCoordinatesPointers", "FlowSolutionPointers")
        # select_timestep itself is intentionally unchanged by this review.
        expected_error = name == "StoredSteps" ? ArgumentError : ErrorException
        @test_throws expected_error BcubeCGNS.read_tree(timestep_remove(name); timestep = 1)
    end
    # A TimeAccurate marker cannot substitute for missing iterative metadata.
    marker_only = edit_reader_tree(
        timestep_fixture(),
        n -> begin
            filter!(
                c -> !(
                    ReaderCGNS.get_label(c) in (
                        ReaderCGNS.BaseIterativeData_t,
                        ReaderCGNS.ZoneIterativeData_t,
                    )
                ),
                reader_children(n),
            )
            n
        end,
    )
    @test_throws ArgumentError BcubeCGNS.read_tree(marker_only)
    # Corrupt step count, per-step array length, exact pointer target or empty
    # name independently; none may silently select the first available child.
    for (name, value) in (
        ("StoredSteps", [0]),
        ("StoredSteps", [2.0]),
        ("TimeValues", [9.0]),
        ("IterationValues", [90]),
        ("FlowSolutionPointers", ["SolutionZulu"]),
        ("GridCoordinatesPointers", ["Missing", "GridAlpha"]),
        ("FlowSolutionPointers", ["SolutionZuluSuffix", "SolutionAlpha"]),
        ("FlowSolutionPointers", ["", "SolutionAlpha"]),
    )
        expected_error = if name in ("StoredSteps", "TimeValues", "IterationValues")
            ArgumentError
        else
            ErrorException
        end
        @test_throws expected_error BcubeCGNS.read_tree(
            timestep_value(name, value);
            timestep = 1,
        )
    end
    # Null means no stored solution; never fall back to an unrelated solution.
    tree = timestep_value("FlowSolutionPointers", ["Null", "SolutionAlpha"])
    @test isempty(BcubeCGNS.read_tree(tree; timestep = 1, varnames = "*").data)
    @test_throws ErrorException BcubeCGNS.read_tree(
        timestep_value("GridCoordinatesPointers", ["Null", "GridAlpha"]);
        timestep = 1,
    )
    # Original cell IDs need not be 1:ZoneCellCount, including at a timestep.
    shifted = timestep_edit(
        "Quads",
        n -> begin
            for c in reader_children(n)
                reader_name(c) == "ElementRange" && (reader_value(c) .+= 10)
            end
            n
        end,
    )
    shifted_result = BcubeCGNS.read_tree(shifted; timestep = 1, varnames = "*")
    @test collect(Bcube.get_absolute_cell_indices(shifted_result.mesh)) == [11, 12]
    @test Bcube.get_values(shifted_result.data["SolutionZulu"]["cell_id"]) == [110.0, 120.0]
    # An inactive zone at the selected step cannot yield this zone's mesh.
    absent_zone = timestep_edit(
        "StoredSteps",
        n -> begin
            push!(
                reader_children(n),
                reader_node("ZonePointers", ["Null" "Missing"], [], "DataArray_t"),
            )
            n
        end,
    )
    @test_throws ErrorException BcubeCGNS.read_tree(absent_zone; timestep = 1)
    # Validate only the selected full-zone solution, including its array length.
    tree = timestep_edit(
        "SolutionZulu",
        n -> begin
            filter!(c -> reader_name(c) != "cell_id", reader_children(n))
            push!(reader_children(n), reader_node("cell_id", [1.0], [], "DataArray_t"))
            n
        end,
    )
    @test_throws ErrorException BcubeCGNS.read_tree(tree; timestep = 1, varnames = "*")
    @test Bcube.get_values(
        BcubeCGNS.read_tree(tree; timestep = 2, varnames = "*").data["SolutionAlpha"]["cell_id"],
    ) == [10.0, 20.0]
    # A pointer-selected Vertex solution keeps node order and location semantics.
    vertex = timestep_value("cell_id", collect(1.0:6.0))
    vertex = edit_reader_tree(
        vertex,
        n -> if reader_name(n) == "GridLocation"
            reader_node("GridLocation", "Vertex", [], "GridLocation_t")
        else
            n
        end,
    )
    nodal =
        BcubeCGNS.read_tree(vertex; timestep = 1, varnames = "*").data["SolutionZulu"]["cell_id"]
    @test Bcube.get_values(nodal) == collect(1.0:6.0)
    @test Bcube.get_location(nodal) isa Bcube.PointData
    # Selection does not relax the exclusion of general FaceCenter solutions.
    face = timestep_edit(
        "SolutionZulu",
        n -> begin
            filter!(
                c -> ReaderCGNS.get_label(c) != ReaderCGNS.GridLocation_t,
                reader_children(n),
            )
            push!(
                reader_children(n),
                reader_node("GridLocation", "FaceCenter", [], "GridLocation_t"),
            )
            n
        end,
    )
    face_error = try
        BcubeCGNS.read_tree(face; timestep = 1, varnames = "*")
        nothing
    catch e
        e
    end
    @test face_error isa ErrorException && occursin(
        "FlowSolution 'SolutionZulu': unsupported GridLocation",
        sprint(showerror, face_error),
    )
    # Active motion cannot be applied; Null at step 2 requires no transform.
    motion = timestep_edit(
        "StepMap",
        n -> begin
            push!(
                reader_children(n),
                reader_node(
                    "RigidGridMotionPointers",
                    ["Motion", "Null"],
                    [],
                    "DataArray_t",
                ),
            )
            n
        end,
    )
    @test_throws ArgumentError BcubeCGNS.read_tree(motion; timestep = 1)
    @test isnothing(BcubeCGNS.read_tree(motion; timestep = 2).data)
    # Selecting a step is non-mutating and can be repeated in any order.
    tree = timestep_fixture()
    before = deepcopy(
        reader_value(first(ReaderCGNS.get_nodes_from_name(tree, "FlowSolutionPointers"))),
    )
    BcubeCGNS.read_tree(tree; timestep = 2, varnames = "*")
    @test Bcube.get_values(
        BcubeCGNS.read_tree(tree; timestep = 1, varnames = "*").data["SolutionZulu"]["cell_id"],
    ) == [110.0, 120.0]
    @test reader_value(
        first(ReaderCGNS.get_nodes_from_name(tree, "FlowSolutionPointers")),
    ) == before

    # Constant grids/solutions do not require pointers, even with iterative
    # metadata. A single stored step still requires explicit selection.
    constant = edit_reader_tree(
        reader_fixture(),
        n -> begin
            if ReaderCGNS.get_label(n) == ReaderCGNS.CGNSBase_t
                push!(reader_children(n), reader_node("Steps", [1], [], "BaseIterativeData_t"))
            elseif ReaderCGNS.get_label(n) == ReaderCGNS.Zone_t
                filter!(c -> reader_name(c) != "Nodes", reader_children(n))
            end
            n
        end,
    )
    @test_throws ArgumentError BcubeCGNS.read_tree(constant)
    @test Bcube.get_values(
        BcubeCGNS.read_tree(constant; timestep = 1, varnames = "*").data["Cells"]["cell_id"],
    ) == [10.0, 20.0]
    # Time selection must retain every element section for whole-MIXED checks,
    # including unsupported lower-dimensional/unreferenced records.
    tree = edit_reader_tree(
        timestep_fixture(),
        n -> begin
            if ReaderCGNS.get_label(n) == ReaderCGNS.Zone_t
                push!(
                    reader_children(n),
                    reader_node(
                        "Unsupported",
                        [20, 0],
                        [
                            reader_node("ElementRange", [99, 99], [], "IndexRange_t"),
                            reader_node(
                                "ElementConnectivity",
                                [22, 1, 2, 3],
                                [],
                                "DataArray_t",
                            ),
                        ],
                        "Elements_t",
                    ),
                )
            end
            n
        end,
    )
    err = try
        BcubeCGNS.read_tree(tree; timestep = 1)
        nothing
    catch e
        e
    end
    @test err isa ErrorException && occursin("unsupported embedded", sprint(showerror, err))
end

# Reuse opposite child/pointer order and marker/grid pairing from the independent
# two-step fixture. Activation metadata is added literally, not by the writer.
function timestep_activation(counts, pointers)
    timestep_edit(
        "StoredSteps",
        n -> begin
            isnothing(counts) || push!(
                reader_children(n),
                reader_node("NumberOfZones", counts, [], "DataArray_t"),
            )
            isnothing(pointers) || push!(
                reader_children(n),
                reader_node("ZonePointers", pointers, [], "DataArray_t"),
            )
            n
        end,
    )
end
@testset "Stored-step zone activation" begin
    @testset "Zone target $target" for target in ("TwoQuads", "Base/TwoQuads")
        tree = timestep_activation([1, 1], [target target])
        for step in 1:2
            result = BcubeCGNS.read_tree(tree; timestep = step, varnames = "*")
            name = step == 1 ? "SolutionZulu" : "SolutionAlpha"
            @test Bcube.get_values(result.data[name]["cell_id"]) ==
                  (step == 1 ? [110.0, 120.0] : [10.0, 20.0])
            @test collect(Bcube.get_coords(first(Bcube.get_nodes(result.mesh)))) ==
                  (step == 1 ? [10.0, 0.0] : [0.0, 0.0])
        end
    end
    # Qualified names contain exactly one slash and two legal CGNS names;
    # pointers into another base/zone cannot silently return the available mesh.
    for target in (
        "Other/TwoQuads",
        "Base/Other",
        "/Base/TwoQuads",
        "Base//TwoQuads",
        "Base/TwoQuads/More",
        "Base/",
        "/TwoQuads",
        "",
    )
        @test_throws ErrorException BcubeCGNS.read_tree(
            timestep_activation([1, 1], [target "TwoQuads"]);
            timestep = 1,
        )
    end
    # Counts are a one-dimensional integer array, one per stored step. Only
    # zero/one active zone is representable; Bool and floating counts are invalid.
    for counts in ([1.0, 1.0], [true, true], [1], [-1, 1], [2, 1], [1 1])
        @test_throws ErrorException BcubeCGNS.read_tree(
            timestep_activation(counts, ["TwoQuads" "TwoQuads"]);
            timestep = 1,
        )
    end
    # Contradictions are errors even when the selected step itself looks valid.
    for (counts, pointers) in (
        ([0, 1], ["TwoQuads" "TwoQuads"]),
        ([1, 1], ["Null" "TwoQuads"]),
        ([1, 0], ["TwoQuads" "TwoQuads"]),
    )
        @test_throws ErrorException BcubeCGNS.read_tree(
            timestep_activation(counts, pointers);
            timestep = 1,
        )
    end
    # A genuinely inactive step cannot produce a mesh; the next active step
    # remains readable. With no pointers, a count of one uses the single zone.
    tree = timestep_activation([0, 1], ["Null" "Base/TwoQuads"])
    @test_throws ErrorException BcubeCGNS.read_tree(tree; timestep = 1)
    @test ncells(BcubeCGNS.read_tree(tree; timestep = 2).mesh) == 2
    @test_throws ErrorException BcubeCGNS.read_tree(
        timestep_activation([0, 1], nothing);
        timestep = 1,
    )
    @test ncells(
        BcubeCGNS.read_tree(timestep_activation([1, 1], nothing); timestep = 1).mesh,
    ) == 2
    @test ncells(
        BcubeCGNS.read_tree(
            timestep_activation(nothing, ["Base/TwoQuads" "TwoQuads"]);
            timestep = 1,
        ).mesh,
    ) == 2
end

@testset "Stored-step qualified-name bounds" begin
    # Base-level pointers allow 65 bytes: two maximum-length CGNS names and '/'.
    # Zonal grid/solution pointers still obey their separate 32-byte contract.
    basename, zonename = repeat("B", 32), repeat("Z", 32)
    tree = timestep_activation([1, 1], fill(basename * "/" * zonename, 1, 2))
    tree = edit_reader_tree(
        tree,
        n -> begin
            label = ReaderCGNS.get_label(n)
            if label in (ReaderCGNS.CGNSBase_t, ReaderCGNS.Zone_t)
                return reader_node(
                    label == ReaderCGNS.CGNSBase_t ? basename : zonename,
                    reader_value(n),
                    reader_children(n),
                    ReaderCGNS.enum_to_string(label),
                )
            end
            n
        end,
    )
    @test Bcube.get_values(
        BcubeCGNS.read_tree(tree; timestep = 1, varnames = "*").data["SolutionZulu"]["cell_id"],
    ) == [110.0, 120.0]
    # Wrong slot/step dimensions and overlong path tokens are malformed even
    # when the remaining metadata names the supported zone.
    for pointers in (
        reshape(["TwoQuads"], 1, 1),
        fill("TwoQuads", 2, 2),
        [repeat("B", 33) * "/TwoQuads" "TwoQuads"],
    )
        @test_throws ErrorException BcubeCGNS.read_tree(
            timestep_activation([1, 1], pointers);
            timestep = 1,
        )
    end
end

@testset "Latest stored step and metadata-free solutions" begin
    # TimeValues=[9,2]: latest means the last stored pointer, not maximum time.
    latest = read_file(
        joinpath(assets_path, "timesteps.cgns");
        timestep = :latest,
        varnames = "*",
    )
    @test Set(keys(latest.data)) == Set(["SolutionAlpha"])
    @test Bcube.get_values(latest.data["SolutionAlpha"]["cell_id"]) == [10.0, 20.0]
    @test collect(Bcube.get_coords(first(Bcube.get_nodes(latest.mesh)))) == [0.0, 0.0]
    # The independent two-quad fixture has both Vertex and CellCenter solutions
    # and no iterative metadata. Neither solution may disappear in the fallback.
    allsolutions = BcubeCGNS.read_tree(reader_fixture(); timestep = nothing, varnames = "*")
    @test Set(keys(allsolutions.data)) == Set(["Nodes", "Cells"])
    fallback = @test_logs (:warn, r"No iterative metadata") BcubeCGNS.read_tree(
        reader_fixture();
        timestep = :latest,
        varnames = "*",
    )
    @test Bcube.get_values(fallback.data["Nodes"]["node_id"]) == collect(1.0:6.0)
    @test Bcube.get_values(fallback.data["Cells"]["cell_id"]) == [10.0, 20.0]
    @test_logs (:warn, r"No iterative metadata") BcubeCGNS.read_tree(
        reader_fixture();
        timestep = :latest,
    )
    @test_throws ArgumentError BcubeCGNS.read_tree(
        timestep_value("StoredSteps", [0]);
        timestep = :latest,
    )
    @test_throws ErrorException BcubeCGNS.read_tree(
        timestep_activation([1, 0], ["TwoQuads" "Null"]);
        timestep = :latest,
    )
    @test_throws ArgumentError read_file(
        joinpath(assets_path, "timesteps.cgns");
        timestep = :unknown,
    )
end
