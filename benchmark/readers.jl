# Optional warmed-up reader benchmarks; no dependency beyond the package's own
# environment. Run this same script in separate, sequential Julia processes
# with --project pointing to each checkout. Inputs always come from assets_dir.
using Bcube, BcubeCGNS
using HDF5, JLD2
using Statistics: median
using Printf: @printf

function benchmark_readers(assets_dir, samples)
    samples >= 3 || error("Use at least three samples")
    # All six fixtures present on main before reader coverage development.
    fixtures = [
        "rectangle_mesh",
        "naca12-surfBase-unstructured",
        "naca12-eclipps-unstructured",
        "naca12-eclipps-unstructured-multielts",
        "RG-15-eclipps-unstructured",
        "naca12-eclipps-structured",
    ]
    println(
        "# Julia=",
        VERSION,
        " Bcube=",
        pkgversion(Bcube),
        " threads=",
        Threads.nthreads(),
    )
    println("# package=", pathof(BcubeCGNS), " samples=", samples)
    checkout = dirname(dirname(pathof(BcubeCGNS)))
    println("# commit=", readchomp(`git -C $checkout rev-parse HEAD`))
    println("# HDF5=", pkgversion(HDF5), " JLD2=", pkgversion(JLD2), " CPU=", Sys.CPU_NAME)
    println(
        "fixture,operation,status,median_ms,min_ms,median_bytes,median_allocations,cells,nodes",
    )
    for fixture in fixtures
        file = joinpath(assets_dir, fixture * ".cgns")
        for mode in ("mesh", "all_fields")
            measure_reader(
                fixture,
                mode,
                () -> read_file(file; varnames = mode == "mesh" ? nothing : "*"),
                samples,
            )
        end
        # These are the original Eclipps reader workloads. Its API always reads
        # BC datasets, so there is no equivalent mesh-only extraction option.
        if occursin("eclipps", fixture)
            bc = if startswith(fixture, "RG-15")
                "AEROFOIL"
            elseif fixture == "naca12-eclipps-structured"
                "WALL1"
            else
                "WALL"
            end
            measure_reader(
                fixture,
                "surface_fields",
                () -> BcubeCGNS.extract_surf_from_eclipps(file, bc),
                samples,
            )
        end
    end
end

function measure_reader(fixture, operation, read, samples)
    try
        result = read() # Compile and warm HDF5/OS caches outside measurements.
        timings, bytes, allocations = Float64[], Int[], Int[]
        for _ in 1:samples
            GC.gc()
            sample = @timed read()
            push!(timings, sample.time * 1000)
            push!(bytes, sample.bytes)
            stats = sample.gcstats
            push!(
                allocations,
                stats.malloc + stats.realloc + stats.poolalloc + stats.bigalloc,
            )
        end
        @printf(
            "%s,%s,ok,%.3f,%.3f,%.0f,%.0f,%d,%d\n",
            fixture,
            operation,
            median(timings),
            minimum(timings),
            median(bytes),
            median(allocations),
            ncells(result.mesh),
            nnodes(result.mesh)
        )
    catch e
        # Failures remain visible and comparable; never time a different fallback.
        println(fixture, ",", operation, ",unsupported,,,,,,")
        println(stderr, fixture, " / ", operation, ": ", sprint(showerror, e))
    end
    flush(stdout)
end

assets_dir =
    isempty(ARGS) ? normpath(joinpath(@__DIR__, "..", "test", "assets")) : abspath(ARGS[1])
samples = length(ARGS) < 2 ? 10 : parse(Int, ARGS[2])
benchmark_readers(assets_dir, samples)
