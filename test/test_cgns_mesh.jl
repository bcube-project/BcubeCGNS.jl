@testset "CGNSMesh" begin
    @testset "Rectangle mesh (3×4)" begin
        bcube_mesh = read_mesh(filepath)
        dΩ_bcube = Measure(CellDomain(bcube_mesh), 1)
        vol_bcube = Bcube.compute(∫(PhysicalFunction(x -> 1.0))dΩ_bcube)

        tree = BcubeCGNS.CGNS.read_cgns_file(filepath)
        cgns_mesh = BcubeCGNS.CGNSMesh(tree)
        Ω_cgns = CellDomain(cgns_mesh)
        @test length(Ω_cgns) == 1
        dΩ_cgns = Measure(Ω_cgns, 1)
        vol_cgns = Bcube.compute(∫(PhysicalFunction(x -> 1.0))dΩ_cgns)

        @test ncells(bcube_mesh) == ncells(cgns_mesh)
        @test vol_bcube ≈ vol_cgns
    end
end