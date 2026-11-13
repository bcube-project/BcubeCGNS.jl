# Literal children distinguish exact names, regex patterns, the parent and a
# grandchild. Duplicate names exercise existence and uniqueness independently.
@testset "CGNS child lookup" begin
    C = BcubeCGNS.CGNS
    child(name, children = C.CGNSNode[]) =
        C.CGNSNode(name, nothing, children, "DataArray_t")
    exact = child("ElementRange")
    suffix = child("ElementRangeExtra")
    root = child("ElementRange", [exact, suffix, child("Branch", [child("Nested")])])
    @test C.get_children_from_name(root, "ElementRange") == [exact]
    @test C.get_children_from_name(root, r"^ElementRange") == [exact, suffix]
    @test_throws ArgumentError C.get_child(root)
    @test C.get_child(root; name = "Nested") === nothing
    @test C.get_child(root; name = "ElementRange", required = true, unique = true) === exact
    @test C.get_child_from_name_and_label(
        root,
        "ElementRange",
        "DataArray_t";
        required = true,
        unique = true,
    ) === exact
    @test_throws ErrorException C.get_child(root; name = "Missing", required = true)
    @test C.get_child(root; name = "Missing", unique = true) === nothing
    push!(C.get_children(root), child("ElementRange"))
    @test C.get_child(root; name = "ElementRange") === exact
    @test_throws ErrorException C.get_child(root; name = "ElementRange", unique = true)
    @test_throws ErrorException C.get_child_from_predicate(
        root,
        n -> C.get_name(n) == "ElementRange";
        unique = true,
    )
    # Recursive node search still includes the root; string names remain exact.
    @test length(C.get_nodes_from_name(root, "ElementRange")) == 3
    @test length(C.get_nodes_from_name(root, r"^ElementRange")) == 4
end

@testset "Zone dimension inference" begin
    C = BcubeCGNS.CGNS
    # Both layouts use 32-bit metadata; the fixed-size result uses machine Int
    # independently of zone kind and avoids a runtime-length Tuple conversion.
    for (kind, sizes) in (
        ("Unstructured", reshape(Int32[6, 2, 0], 1, 3)),
        ("Structured", Int32[3 2 0; 2 1 0]),
    )
        zone = C.CGNSNode(
            "Zone",
            sizes,
            [C.CGNSNode("ZoneType", kind, C.CGNSNode[], "ZoneType_t")],
            "Zone_t",
        )
        @test (@inferred BcubeCGNS.get_zone_dims(zone)) === (6, 2, 0)
    end
end
