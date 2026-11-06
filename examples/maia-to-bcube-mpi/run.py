from mpi4py import MPI
from maia.factory import generate_dist_block
from maia.factory import partition_dist_tree
from maia.transfer import part_tree_to_dist_tree_all
from maia.algo import compute_elements_center
from maia.io import write_tree, dist_tree_to_file
import maia.pytree as PT
import maia.pytree.maia as MT
from maia_pdm_partition_extension import extend_partition_3d
from juliacall import Main as jl
import numpy as np

# MPI infos
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Build (or read) dist. mesh
eltype = 'HEXA_8'
Lx = 1.
Ly = 1.
Lz = 1.
dist_tree = generate_dist_block((21,2,2), eltype, comm, origin=(-Lx/2, -Ly/2, -Lz/2))

# Initial solution
compute_elements_center(dist_tree, "CellCenter", comm)
for zone in PT.iter_all_Zone_t(dist_tree):
    geom = PT.get_child_from_predicate(zone, lambda n: PT.get_name(n).startswith("Geometry") and PT.get_label(n) == "DiscreteData_t")
    x, y, z = [PT.get_value(n) for n in PT.get_children_from_label(geom, "DataArray_t")]
    u = np.zeros_like(x)
    u[np.where(np.abs(x) < 0.3)[0]] = 1.
    PT.new_FlowSolution('FS', loc='CellCenter', fields={'u': u}, parent=zone)

# Partition the tree
part_tree = partition_dist_tree(dist_tree, comm, save_all_connectivities=True, dump_pdm_output=True, data_transfer=["FlowSolution_t"])

# Compute partition extension
extend_partition_3d(part_tree, comm)

# Build ghost tag => part
ghost_tag2part = {}
pred = lambda n: PT.get_name(n).startswith("GHOST_") and PT.get_label(n) == "Elements_t"
for elt_node in PT.get_nodes_from_predicate(part_tree, pred):
    cell_ext_ln_to_gn = PT.get_value(MT.get_GlobalNumbering(elt_node, 'Element'))
    cell_to_rank = PT.get_value(MT.get_GlobalNumbering(elt_node, 'Rank'))
    cell_to_part = cell_to_rank + 1
    ghost_tag2part.update(zip(cell_ext_ln_to_gn, cell_to_part))

# Transfer the mesh and initial field to Bcube, and setup simulation
jl.include("run.jl")
bcube_io = jl.tree_to_bcube(part_tree, ghost_tag2part, comm.py2f())
bcube_solver = jl.setup_simulation(bcube_io, np.array([1., 0., 0.]))

# Time loop!
u_part = PT.get_value(PT.get_node_from_name_and_label(part_tree, "u", "DataArray_t")) # Value of `u` in the CGNS part. tree
for i in range(10):
    jl.step_forward(bcube_solver)
    print(u_part) # CGNS node is modified in memory by Bcube!

# Write the final solution to CGNS
part_tree_to_dist_tree_all(dist_tree, part_tree, comm)
dist_tree_to_file(dist_tree, "tmp/final.cgns", comm)
