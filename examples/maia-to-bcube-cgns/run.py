from mpi4py import MPI
import maia
import maia.pytree as PT
from juliacall import Main as jl
import os
print(os.environ["PYTHON_JULIACALL_PROJECT"])

dist_tree = maia.factory.generate_dist_block(3, 'Poly', MPI.COMM_WORLD)
maia.algo.dist.convert_ngon_to_elements(dist_tree, MPI.COMM_WORLD)

jl.include("run.jl")
jl.translate(dist_tree)
# jl.test(dist_tree)
