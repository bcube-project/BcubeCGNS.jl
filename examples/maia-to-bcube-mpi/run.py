from mpi4py import MPI
from maia.factory import generate_dist_block
from maia.factory import partition_dist_tree
import maia.pytree as PT
import maia.pytree.maia as MT
from maia_pdm_partition_extension import extend_partition_3d
from juliacall import Main as jl

# MPI infos
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
my_part = rank + 1

# Build (dist) mesh
eltype = 'HEXA_8'
dist_tree = generate_dist_block((5,2,2), eltype, comm, origin = (0,10,20))

# Partition the tree
part_tree = partition_dist_tree(dist_tree, comm, save_all_connectivities=True, dump_pdm_output=True)

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

# Transfer the mesh to Bcube
jl.include("run.jl")
mesh = jl.build_bcube_dmesh(part_tree, ghost_tag2part, comm)
