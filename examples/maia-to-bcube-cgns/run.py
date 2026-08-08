from mpi4py import MPI
import maia
from juliacall import Main as jl

# Build a distributed tree with Maia, and convert NGONS to elements
dist_tree = maia.factory.generate_dist_block(3, 'Poly', MPI.COMM_WORLD)
maia.algo.dist.convert_ngon_to_elements(dist_tree, MPI.COMM_WORLD)

# Build the corresponding Bcube mesh from this tree
jl.include("run.jl")
mesh = jl.build_bcube_mesh(dist_tree)

# Compute and display the gravity center
jl.compute_gravity_center(mesh)
# > x = [0.5, 0.5, 0.5]

# Translate, with maia, the tree
zone = maia.pytree.get_all_Zone_t(dist_tree)[0]
maia.algo.transform_affine(zone, translation=[3,2,0])

# Compute and display the gravity center, note that mesh is not "updated" before!
jl.compute_gravity_center(mesh)
# > x = [3.5, 2.5, 0.5]

