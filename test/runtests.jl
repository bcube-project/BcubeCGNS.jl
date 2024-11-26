using Test
using BcubeCGNS
using Bcube

"""
Custom way to "include" a file to print infos.
"""
function custom_include(path)
    filename = split(path, "/")[end]
    print("Running test file " * filename * "...")
    include(path)
    println("done.")
end

# This dir will be removed at the end of the tests
tempdir = mktempdir()

# Fake test
@test true