@testset "CGNSMesh" begin
    function common_check(filepath, n_subdomains)
        bcube_mesh = read_mesh(filepath)
        dΩ_bcube = Measure(CellDomain(bcube_mesh), 1)
        vol_bcube = sum(Bcube.compute(∫(PhysicalFunction(x -> 1.0))dΩ_bcube))

        tree = BcubeCGNS.CGNS.read_cgns_file(filepath)
        cgns_mesh = BcubeCGNS.CGNSMesh(tree)
        Ω_cgns = CellDomain(cgns_mesh)
        @test length(Ω_cgns.subdomains) == n_subdomains
        dΩ_cgns = Measure(Ω_cgns, 1)
        vol_cgns = sum(Bcube.compute(∫(PhysicalFunction(x -> 1.0))dΩ_cgns))

        @test ncells(bcube_mesh) == ncells(cgns_mesh)
        @test vol_bcube ≈ vol_cgns
    end

    @testset "Rectangle mesh (3×4)" begin
        common_check(joinpath(assets_path, "rectangle_mesh.cgns"), 1)
    end

    @testset "NACA12 Multielts" begin
        common_check(
            joinpath(assets_path, "naca12-eclipps-unstructured-multielts.cgns"),
            28,
        )
    end
end