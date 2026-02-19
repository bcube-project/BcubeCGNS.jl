@testset "Write" begin
    filepath = joinpath(tempdir, "output.cgns")

    mesh = rectangle_mesh(3, 2)
    U = TrialFESpace(FunctionSpace(:Lagrange, 1), mesh)
    u = FEFunction(U, 1.0 * collect(1:get_ndofs(U)))
    # write_file(filepath, mesh, Dict("Temperature" => u))
    write_file(filepath, mesh, Dict("Temperature" => u), 1)

    result = read_file(filepath; varnames = "*")
    @test get_values(result.data["FlowSolutionVertex#1"]["Temperature"]) ==
          [1.0, 2.0, 5.0, 3.0, 4.0, 6.0]

    set_dof_values!(u, 2.0 * collect(1:get_ndofs(U)))
    write_file(filepath, mesh, Dict("Temperature" => u), 2; collection_append = true)

    result = read_file(filepath; varnames = "*")
    @test get_values(result.data["FlowSolutionVertex#2"]["Temperature"]) ==
          2 .* [1.0, 2.0, 5.0, 3.0, 4.0, 6.0]
end