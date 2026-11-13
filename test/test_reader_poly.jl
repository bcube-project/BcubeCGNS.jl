# Independent serial Maia fixture: two unit cubes along X, eleven NGON faces
# (including shared face 3), NFACE cells 12:13. Negative reference -3 reverses
# that shared loop for the right cell. End BCs are faces 5 and 9. All expected
# coordinates, signed unit volumes, face sets and markers are hand checked.
poly_fixture() = ReaderCGNS.read_cgns_file(joinpath(assets_path, "poly-cubes.cgns"))
function poly_edit(edit; name = "Cells")
    edit_reader_tree(poly_fixture(), n -> reader_name(n)==name ? edit(n) : n)
end
@testset "NGON/NFACE cubes" begin
    r = read_file(joinpath(assets_path, "poly-cubes.cgns"); varnames = "*")
    m = r.mesh
    @test typeof.(Bcube.cells(m)) == [Bcube.Hexa8_t, Bcube.Hexa8_t]
    c = Bcube.connectivities_indices(m, :c2n)
    @test [collect(v) for v in c] == [[1, 2, 5, 4, 7, 8, 11, 10], [2, 3, 6, 5, 8, 9, 12, 11]]
    p = [collect(Bcube.get_coords(n)) for n in Bcube.get_nodes(m)]
    @test p == [[x, y, z] for z in (0.0, 1.0) for y in (0.0, 1.0) for x in (0.0, 1.0, 2.0)]
    @test [sum(p[v] for v in c[i])/8 for i in 1:2] == [[0.5, 0.5, 0.5], [1.5, 0.5, 0.5]]
    # The oriented corner triple equals the signed volume for these affine cubes.
    @test [
        sum(
            (p[c[i][2]]-p[c[i][1]]) .*
            BcubeCGNS.cross(p[c[i][4]]-p[c[i][1]], p[c[i][5]]-p[c[i][1]]),
        ) for i in 1:2
    ] == [1.0, 1.0]
    # Face IDs 1:11 precede cells 12:13. Compact fields still have two values;
    # Ends selects only the x=0/x=2 faces, not the shared face or side faces.
    @test collect(Bcube.get_absolute_cell_indices(m)) == [12, 13]
    @test Bcube.get_values(r.data["CellValues"]["marker"]) == [12.0, 13.0]
    @test structured_faces(poly_fixture(), "Ends") == [[1, 4, 7, 10], [3, 6, 9, 12]]
    # Arbitrary ordering of face references must not affect the deterministic Bcube assignment or cell field order.
    t = poly_edit(
        n -> begin
            a=reader_value(ReaderCGNS.get_child(n; name = "ElementConnectivity"))
            a[1:6] = reverse(a[1:6])
            a[7:12]=reverse(a[7:12])
            n
        end,
    )
    @test [
        collect(v) for v in Bcube.connectivities_indices(BcubeCGNS.read_tree(t).mesh, :c2n)
    ] == [collect(v) for v in c]
end

# Build tiny literal CGNS trees, independently of the package writer. Face loops
# below are outward on a reference tetrahedron, square pyramid and right prism.
# A tetra's signed determinant is +1 (volume 1/6), pyramid volume 1/3 and prism
# volume 1/2. Exact connectivity catches reflection even when centers agree.
function poly_literal(points, loops; volume = true)
    data(name, v, label = "DataArray_t") = reader_node(name, v, [], label)
    section(name, code, ids, records) = reader_node(
        name,
        [code, 0],
        [
            data("ElementRange", ids, "IndexRange_t"),
            data("ElementConnectivity", vcat(records...)),
            data("ElementStartOffset", vcat(0, cumsum(length.(records)))),
        ],
        "Elements_t",
    )
    children = [
        data("ZoneType", "Unstructured", "ZoneType_t"),
        reader_node(
            "GridCoordinates",
            nothing,
            [
                data("Coordinate"*string(a), [p[i] for p in points]) for
                (i, a) in enumerate("XYZ"[1:length(first(points))])
            ],
            "GridCoordinates_t",
        ),
        section("Faces", 22, [1, length(loops)], loops),
    ]
    volume && push!(
        children,
        section(
            "Cells",
            23,
            [length(loops)+1, length(loops)+1],
            [collect(1:length(loops))],
        ),
    )
    zone=reader_node(
        "Zone",
        reshape([length(points), volume ? 1 : length(loops), 0], 1, 3),
        children,
        "Zone_t",
    )
    reader_node(
        "CGNSTree",
        nothing,
        [
            reader_node(
                "Base",
                [volume ? 3 : 2, length(first(points))],
                [zone],
                "CGNSBase_t",
            ),
        ],
        "CGNSTree_t",
    )
end
@testset "Representable polyhedral shapes" begin
    for (points, loops, T, expected_volume) in (
        (
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [[1, 3, 2], [1, 2, 4], [2, 3, 4], [3, 1, 4]],
            Bcube.Tetra4_t,
            1/6,
        ),
        (
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.5, 0.5, 1.0],
            ],
            [[1, 4, 3, 2], [1, 2, 5], [2, 3, 5], [3, 4, 5], [4, 1, 5]],
            Bcube.Pyra5_t,
            1/3,
        ),
        (
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 1.0],
                [0.0, 1.0, 1.0],
            ],
            [[1, 3, 2], [4, 5, 6], [1, 2, 5, 4], [2, 3, 6, 5], [3, 1, 4, 6]],
            Bcube.Penta6_t,
            1/2,
        ),
    )
        m=BcubeCGNS.read_tree(poly_literal(points, loops)).mesh
        @test only(Bcube.cells(m)) isa T
        @test collect(Bcube.connectivities_indices(m, :c2n)[1]) == collect(1:length(points))
        # Divergence theorem on hand-written face loops, triangulated from the
        # first corner; use the returned cell ordering, not source node order.
        q = points[collect(Bcube.connectivities_indices(m, :c2n)[1])]
        determinant(a, b, c) =
            a[1]*(b[2]*c[3]-b[3]*c[2]) - a[2]*(b[1]*c[3]-b[3]*c[1]) +
            a[3]*(b[1]*c[2]-b[2]*c[1])
        volume = sum(
            determinant(q[f[1]], q[f[i]], q[f[i + 1]])/6 for f in loops for
            i in 2:(length(f) - 1)
        )
        @test volume ≈ expected_volume
    end
    # Octahedron: every face is a valid triangle, but the shell is not any
    # existing Bcube linear entity. Counts alone must never trigger conversion.
    octa = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ]
    shell = [
        [1, 2, 5],
        [2, 3, 5],
        [3, 4, 5],
        [4, 1, 5],
        [2, 1, 6],
        [3, 2, 6],
        [4, 3, 6],
        [1, 4, 6],
    ]
    @test_throws ErrorException BcubeCGNS.read_tree(poly_literal(octa, shell))
    @test_throws ErrorException BcubeCGNS.read_tree(poly_literal(octa, shell[1:7]))
    # Without NFACE, NGON records are planar cells: a CCW square and triangle.
    # Subsequent variants reverse orientation, add a fifth corner, or cross edges.
    points=[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [2.0, 0.0]]
    m=BcubeCGNS.read_tree(poly_literal(points, [[1, 2, 3, 4], [2, 5, 3]]; volume = false)).mesh
    @test typeof.(Bcube.cells(m)) == [Bcube.Quad4_t, Bcube.Tri3_t]
    @test [collect(v) for v in Bcube.connectivities_indices(m, :c2n)] == [[1, 2, 3, 4], [2, 5, 3]]
    @test_throws ErrorException BcubeCGNS.read_tree(
        poly_literal(points, [[1, 4, 3, 2]]; volume = false),
    )
    @test_throws ErrorException BcubeCGNS.read_tree(
        poly_literal(points, [[1, 2, 5, 3, 4]]; volume = false),
    )
    @test_throws ErrorException BcubeCGNS.read_tree(
        poly_literal(points, [[1, 3, 2, 4]]; volume = false),
    )
end
@testset "Malformed NGON/NFACE rejection" begin
    # Each edit starts from the same valid cubes and isolates reference,
    # orientation, distinct-vertex or modern-offset validation before conversion.
    for (section, child, edit) in (
        ("Cells", "ElementConnectivity", a -> (a[1]=99)), # missing face
        ("Cells", "ElementConnectivity", a -> (a[1]=0)), # zero is not a face ID
        ("Cells", "ElementConnectivity", a -> (a[1]=2)), # duplicate face
        ("Cells", "ElementConnectivity", a -> (a[1]=-1)), # one reversed face
        ("Cells", "ElementConnectivity", a -> (a[11]=3)), # same-sign shared face
        ("Cells", "ElementConnectivity", a -> (a[1:6] .*= -1)), # inward shell
        ("Faces", "ElementConnectivity", a -> (a[1]=a[2])), # repeated corner
        ("Faces", "ElementConnectivity", a -> (a[1]=99)), # vertex out of bounds
        ("Faces", "ElementStartOffset", a -> (a[2]=3)), # splits quad into invalid loops
        ("Cells", "ElementStartOffset", a -> (a[1]=1)), # offsets must start at zero
        ("Faces", "ElementStartOffset", a -> pop!(a)), # missing terminal offset
    )
        t=poly_edit(
            n -> begin
                edit(reader_value(ReaderCGNS.get_child(n; name = child)))
                n
            end;
            name = section,
        )
        @test_throws ErrorException BcubeCGNS.read_tree(t)
    end
    # Removing offsets is not a supported legacy polygon encoding.
    t=poly_edit(
        n -> begin
            filter!(c -> reader_name(c)!="ElementStartOffset", reader_children(n))
            n
        end;
        name = "Faces",
    )
    @test_throws ErrorException BcubeCGNS.read_tree(t)
    # Raise one upper corner: adjacent quads become warped, outside planar support.
    t=edit_reader_tree(
        poly_fixture(),
        n -> begin
            reader_name(n)=="CoordinateZ" && (reader_value(n)[11]=1.5)
            n
        end,
    )
    @test_throws ErrorException BcubeCGNS.read_tree(t)
end
