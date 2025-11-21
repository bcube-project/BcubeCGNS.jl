@testset "Read" begin
    @testset "Unstructured" begin
        @testset "CGNS base = 1 2" begin
            basename = "naca12-surfBase-unstructured"

            # read
            filepath = joinpath(@__DIR__, "assets", "$(basename).cgns")
            result = read_file(filepath; varnames = "*")

            # write
            filepath = joinpath(tempdir, "$(basename).vtu")
            write_file(filepath, result.mesh, result.data["FlowSolution"])

            # test
            @test fname2sum[basename] == bytes2hex(open(sha1, filepath))
        end

        @testset "Eclipps" begin
            @testset "CGNS base = 3 3" begin
                basename = "naca12-eclipps-unstructured"

                # read
                filepath = joinpath(@__DIR__, "assets", "$(basename).cgns")
                result = BcubeCGNS.extract_surf_from_eclipps(filepath, "WALL")
                mesh = result.mesh
                data = result.data["IcingWallData:Gas"]

                # interpolate node-centered data on Lagrange P1
                fs = FunctionSpace(:Lagrange, 1)
                U = TrialFESpace(fs, mesh)
                I = invperm(vec(build_node_to_idof(mesh, U)))
                data = Dict(
                    key => FEFunction(U, get_values(value)[I]) for (key, value) in data
                )

                # write
                filepath = joinpath(tempdir, "$(basename).vtu")
                write_file(filepath, mesh, data)

                # test
                @test fname2sum[basename] == bytes2hex(open(sha1, filepath))
            end
        end
    end
end