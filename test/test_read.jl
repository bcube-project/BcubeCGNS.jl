@testset "Read" begin
    @testset "Unstructured" begin
        @testset "CGNS base = 1 2" begin
            # Lineic CGNS file with cell-centered data
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
                @testset "NACA0012 ECLIPPS unstructured node-centered" begin
                    # ECLIPPS unstructured file of a 3D mesh with a NACA0012
                    # BC data are node-centered
                    basename = "naca12-eclipps-unstructured"

                    # read
                    filepath = joinpath(@__DIR__, "assets", "$(basename).cgns")
                    result = BcubeCGNS.extract_surf_from_eclipps(filepath, "WALL")
                    mesh = result.mesh
                    data = result.data["IcingWallData:Gas"]

                    # interpolate node-centered data on Lagrange P1
                    data = Dict(
                        key => Bcube.convert_to_lagrange_P1(mesh, value) for
                        (key, value) in data
                    )

                    # write
                    filepath = joinpath(tempdir, "$(basename).vtu")
                    write_file(filepath, mesh, data)

                    # test
                    @test fname2sum[basename] == bytes2hex(open(sha1, filepath))
                end

                @testset "NACA0012 ECLIPPS unstructured node-centered multi elts" begin
                    # ECLIPPS unstructured file of a 3D mesh with a NACA0012
                    # BC data are node-centered
                    basename = "naca12-eclipps-unstructured-multielts"

                    # read
                    filepath = joinpath(@__DIR__, "assets", "$(basename).cgns")
                    result = BcubeCGNS.extract_surf_from_eclipps(filepath, "WALL")
                    mesh = result.mesh
                    data = result.data["IcingWallData:Gas"]

                    # interpolate node-centered data on Lagrange P1
                    data = Dict(
                        key => Bcube.convert_to_lagrange_P1(mesh, value) for
                        (key, value) in data
                    )

                    # write
                    filepath = joinpath(tempdir, "$(basename).vtu")
                    write_file(filepath, mesh, data)

                    # test
                    @test fname2sum[basename] == bytes2hex(open(sha1, filepath))
                end

                @testset "RG-15 ECLIPPS unstructured face-centered" begin
                    # ECLIPPS unstructured file of a 3D mesh with a RG-15
                    # BC data are cell-centered
                    basename = "RG-15-eclipps-unstructured"

                    # read
                    filepath = joinpath(@__DIR__, "assets", "$(basename).cgns")
                    result = BcubeCGNS.extract_surf_from_eclipps(filepath, "AEROFOIL")
                    mesh = result.mesh
                    data = result.data["IcingWallData:Gas"]

                    # write
                    filepath = joinpath(tempdir, "$(basename).vtu")
                    write_file(filepath, mesh, data)

                    # test
                    @test fname2sum[basename] == bytes2hex(open(sha1, filepath))
                end
            end
        end
    end
end
