from mpi4py import MPI
comm = MPI.COMM_WORLD
from maia.io import write_tree
import maia.pytree as PT
import maia.pytree.maia as MT
from cmaia.algo import combine_to_hexa
import Pypdm.Pypdm as PDM
import numpy as np

def print_dbg(messages, comm = MPI.COMM_WORLD, only_root = False, flush=True, decorate_all=False, decorate_first=False):
    rank = comm.Get_rank()

    _messages = messages.copy()
    if decorate_first:
        _messages[0] = f"{[rank]} " + _messages[0]
    elif decorate_all:
        _messages = map(lambda msg: f"{[rank]} {msg}", messages)

    if only_root and rank == 0:
        for msg in _messages:
                print(msg, flush=flush)
    else:
        for r in range(comm.Get_size()):
            if r == rank:
                for msg in _messages:
                    print(msg, flush=flush)
            comm.barrier()

def concatenate_connectivities_3d(data, data_ext):
    # Unpack
    cell_ln_to_gn = data["cell_ln_to_gn"]
    face_ln_to_gn = data["face_ln_to_gn"]
    vtx_ln_to_gn = data["vtx_ln_to_gn"]

    cell_face_idx = data["cell_face_idx"]
    cell_face = data["cell_face"]
    face_vtx_idx = data["face_vtx_idx"]
    face_vtx = data["face_vtx"]

    vtx_coord = data["vtx_coord"]

    cell_ext_ln_to_gn = data_ext["cell_ln_to_gn"]
    face_ext_ln_to_gn = data_ext["face_ln_to_gn"]
    vtx_ext_ln_to_gn = data_ext["vtx_ln_to_gn"]

    cell_face_ext_idx = data_ext["cell_face_idx"]
    cell_face_ext = data_ext["cell_face"]
    face_vtx_ext_idx = data_ext["face_vtx_idx"]
    face_vtx_ext = data_ext["face_vtx"]

    vtx_coord_ext = data_ext["vtx_coord"]

    n_cell = len(cell_ln_to_gn)
    n_face = len(face_ln_to_gn)
    n_vtx = len(vtx_ln_to_gn)

    n_cell_ext = len(cell_ext_ln_to_gn)
    n_face_ext = len(face_ext_ln_to_gn)
    n_vtx_ext = len(vtx_ext_ln_to_gn)

    total_n_cell = n_cell + n_cell_ext
    total_n_face = n_face + n_face_ext
    total_n_vtx = n_vtx + n_vtx_ext

    # Cells
    total_cell_ln_to_gn = np.concatenate(
        (cell_ln_to_gn, cell_ext_ln_to_gn), axis=0, dtype=PDM.npy_pdm_gnum_dtype)

    total_cell_face_idx = np.zeros(total_n_cell + 1, dtype=np.intc)
    for i in range(n_cell + 1):
        total_cell_face_idx[i] = cell_face_idx[i]
    for i in range(n_cell_ext + 1):
        total_cell_face_idx[n_cell +
                            i] = cell_face_idx[n_cell] + cell_face_ext_idx[i]

    total_cell_face = np.concatenate(
        (cell_face, cell_face_ext), axis=0, dtype=np.intc)

    # Faces
    total_face_ln_to_gn = np.concatenate(
        (face_ln_to_gn, face_ext_ln_to_gn), axis=0, dtype=PDM.npy_pdm_gnum_dtype)

    total_face_vtx_idx = np.zeros(total_n_face + 1, dtype=np.intc)
    for i in range(n_face + 1):
        total_face_vtx_idx[i] = face_vtx_idx[i]
    for i in range(n_face_ext + 1):
        total_face_vtx_idx[n_face + i] = face_vtx_idx[n_face] + face_vtx_ext_idx[i]

    total_face_vtx = np.concatenate(
        (face_vtx, face_vtx_ext), axis=0, dtype=np.intc)

    # Vertices
    total_vtx_ln_to_gn = np.concatenate(
        (vtx_ln_to_gn, vtx_ext_ln_to_gn), axis=0, dtype=PDM.npy_pdm_gnum_dtype)

    total_vtx_coord = np.concatenate((vtx_coord, vtx_coord_ext), axis=0, dtype=np.double)


    # Pack
    data_total ={}
    data_total["cell_ln_to_gn"] = total_cell_ln_to_gn
    data_total["face_ln_to_gn"] = total_face_ln_to_gn
    data_total["vtx_ln_to_gn"] = total_vtx_ln_to_gn

    data_total["cell_face_idx"] = total_cell_face_idx
    data_total["cell_face"] = total_cell_face
    data_total["face_vtx_idx"] = total_face_vtx_idx
    data_total["face_vtx"] = total_face_vtx

    data_total["vtx_coord"] = total_vtx_coord

    print_dbg([
        f"len(face_vtx) = {len(face_vtx)}",
        f"len(face_vtx_ext) = {len(face_vtx_ext)}",
        f"len(total_face_vtx) = {len(total_face_vtx)}",
        ]
        , comm)

    return data_total

def concatenate_connectivities_2d(data, data_ext):
    # Unpack
    face_ln_to_gn = data["face_ln_to_gn"]
    edge_ln_to_gn = data["edge_ln_to_gn"]
    vtx_ln_to_gn = data["vtx_ln_to_gn"]

    face_edge_idx = data["face_edge_idx"]
    face_edge = data["face_edge"]
    edge_vtx_idx = data["edge_vtx_idx"]
    edge_vtx = data["edge_vtx"]

    vtx_coord = data["vtx_coord"]

    face_ext_ln_to_gn = data_ext["face_ln_to_gn"]
    edge_ext_ln_to_gn = data_ext["edge_ln_to_gn"]
    vtx_ext_ln_to_gn = data_ext["vtx_ln_to_gn"]

    face_edge_ext_idx = data_ext["face_edge_idx"]
    face_edge_ext = data_ext["face_edge"]
    edge_vtx_ext_idx = data_ext["edge_vtx_idx"]
    edge_vtx_ext = data_ext["edge_vtx"]

    vtx_coord_ext = data_ext["vtx_coord"]

    n_face = len(face_ln_to_gn)
    n_edge = len(edge_ln_to_gn)
    n_vtx = len(vtx_ln_to_gn)

    n_face_ext = len(face_ext_ln_to_gn)
    n_edge_ext = len(edge_ext_ln_to_gn)
    n_vtx_ext = len(vtx_ext_ln_to_gn)

    total_n_face = n_face + n_face_ext
    total_n_edge = n_edge + n_edge_ext
    total_n_vtx = n_vtx + n_vtx_ext

    # Cells
    total_face_ln_to_gn = np.concatenate(
        (face_ln_to_gn, face_ext_ln_to_gn), axis=0, dtype=PDM.npy_pdm_gnum_dtype)

    total_face_edge_idx = np.zeros(total_n_face + 1, dtype=np.intc)
    for i in range(n_face + 1):
        total_face_edge_idx[i] = face_edge_idx[i]
    for i in range(n_face_ext + 1):
        total_face_edge_idx[n_face +
                            i] = face_edge_idx[n_face] + face_edge_ext_idx[i]

    total_face_edge = np.concatenate(
        (face_edge, face_edge_ext), axis=0, dtype=np.intc)

    # Faces
    total_edge_ln_to_gn = np.concatenate(
        (edge_ln_to_gn, edge_ext_ln_to_gn), axis=0, dtype=PDM.npy_pdm_gnum_dtype)

    total_edge_vtx_idx = np.zeros(total_n_edge + 1, dtype=np.intc)
    for i in range(n_edge + 1):
        total_edge_vtx_idx[i] = edge_vtx_idx[i]
    for i in range(n_edge_ext + 1):
        total_edge_vtx_idx[n_edge + i] = edge_vtx_idx[n_edge] + edge_vtx_ext_idx[i]

    total_edge_vtx = np.concatenate(
        (edge_vtx, edge_vtx_ext), axis=0, dtype=np.intc)

    # Vertices
    total_vtx_ln_to_gn = np.concatenate(
        (vtx_ln_to_gn, vtx_ext_ln_to_gn), axis=0, dtype=PDM.npy_pdm_gnum_dtype)

    total_vtx_coord = np.concatenate((vtx_coord, vtx_coord_ext), axis=0, dtype=np.double)


    # Pack
    data_total ={}
    data_total["face_ln_to_gn"] = total_face_ln_to_gn
    data_total["edge_ln_to_gn"] = total_edge_ln_to_gn
    data_total["vtx_ln_to_gn"] = total_vtx_ln_to_gn

    data_total["face_edge_idx"] = total_face_edge_idx
    data_total["face_edge"] = total_face_edge
    data_total["edge_vtx_idx"] = total_edge_vtx_idx
    data_total["edge_vtx"] = total_edge_vtx

    data_total["vtx_coord"] = total_vtx_coord

    print_dbg([
        f"len(edge_vtx) = {len(edge_vtx)}",
        f"len(edge_vtx_ext) = {len(edge_vtx_ext)}",
        f"len(total_edge_vtx) = {len(total_edge_vtx)}",
        ]
        , comm)

    return data_total

def append_extension_to_tree_3d(zone, data_ext, data_total, cell_to_rank = None):
    # Coordinates
    x, y, z = data_total["vtx_coord"].reshape(-1, 3).T
    gc = PT.get_child_from_label(zone, "GridCoordinates_t")
    for (suffix, arr) in zip(("X", "Y", "Z"), (x, y, z)):
        elt_node = PT.get_child_from_name_and_label(gc, f"Coordinate{suffix}", "DataArray_t")
        PT.set_value(elt_node, arr)

    # Alias
    cell_ext_ln_to_gn = data_ext["cell_ln_to_gn"]
    total_cell_ln_to_gn = data_total["cell_ln_to_gn"]
    n_local_cells = len(total_cell_ln_to_gn)
    n_ghost_cells = len(cell_ext_ln_to_gn)
    n_owned_cells = n_local_cells - n_ghost_cells

    # Cells
    # BIG ASSUMPTION HERE: same element type for all ghosts. It's not impossible to deal with mixed elements, but would first require a sort.
    # However, note that most likely a part-to-part would be necessary to identify the ghost "element-type" because here we only have NGON
    # infos.
    cell_face_ext_idx = data_ext["cell_face_idx"]
    cell_face_ext = data_ext["cell_face"]
    total_face_vtx = data_total["face_vtx"]
    if n_ghost_cells == 0:
        return
    n_faces_by_cell = cell_face_ext_idx[1]
    if n_faces_by_cell == 6: # HEXA_8
        eltype = "HEXA_8"
        n_vtx_by_cell = 8
        n_vtx_by_face = 4
        # In order to call `combine_to_hexa`, we need to build
        # * np_face_vtx_n : number of vertex per face (here, 4)
        # * np_face_vtx : flattened connectivity face->vtx. It must describe
        #   exactly the faces that are "converted" to classical elements, it cannot
        #   be a "bigger" set of faces. In other words, element 1 is described by the 6 first
        #   faces of this array, element 2 by the 6 following etc
        # * np_cell_face : face orientation (so only the sign is important)
        econn = np.zeros(n_ghost_cells*n_vtx_by_cell, dtype = cell_face_ext.dtype)
        np_face_vtx_n = n_vtx_by_face*np.ones(n_ghost_cells*n_faces_by_cell, dtype = cell_face_ext_idx.dtype)
        total_face_vtx_reshaped = np.reshape(total_face_vtx, (-1, 4))
        np_face_vtx = total_face_vtx_reshaped[np.abs(cell_face_ext) - 1, :].flatten() - 1 # 0-based for python
        np_cell_face = cell_face_ext # Only ghost cells, (we just need the sign so no need for 0-based)

        combine_to_hexa(np_face_vtx_n, np_face_vtx, np_cell_face, econn)
        econn += 1

    elif n_faces_by_cell == 4: # TETRA_4
        raise ValueError(f"n_faces_by_cell = {n_faces_by_cell} not implemented yet")
    else:
        raise ValueError(f"n_faces_by_cell = {n_faces_by_cell} not implemented yet")

    # Search the range
    ranges = PT.Zone.get_elt_range_per_dim(zone)
    max_range = max(r[1] for r in ranges)
    erange = [max_range + 1, max_range + n_ghost_cells]

    # Create Element node
    elt_node = PT.new_Elements('GHOST_'+eltype, type=eltype, erange=erange, econn=econn, parent=zone)
    gnum = {'Element' : cell_ext_ln_to_gn}
    if cell_to_rank is not None:
        gnum['Rank'] = cell_to_rank[n_owned_cells::]
    MT.new_GlobalNumbering(gnum, parent = elt_node)

    # Fix Zone infos
    n_vtx_total = len(x)
    n_cell_total = len(total_cell_ln_to_gn)
    PT.set_value(zone, [[n_vtx_total, n_cell_total, 0]])
    PT.set_value(MT.get_GlobalNumbering(zone, 'Vertex'), data_total['vtx_ln_to_gn'])
    PT.set_value(MT.get_GlobalNumbering(zone, 'Cell'), total_cell_ln_to_gn)

    # Fix FlowSolution (increase size and place zeros)
    for fs in PT.get_children_from_label(zone, "FlowSolution_t"):
        gl = PT.Container.GridLocation(fs)
        new_length = n_vtx_total if gl == "Vertex" else n_cell_total
        for darray in PT.get_children_from_label(fs, "DataArray_t"):
            val = PT.get_value(darray)
            new_val = np.zeros(new_length, dtype = val.dtype)
            new_val[0:len(val)] = val
            PT.set_value(darray, new_val)


def append_extension_to_tree_2d(zone, data_ext, data_total):
    ## Coordinates
    x, y, z = data_total["vtx_coord"].reshape(-1, 3).T
    gc = PT.get_child_from_label(zone, "GridCoordinates_t")
    for (suffix, arr) in zip(("X", "Y", "Z"), (x, y, z)):
        node = PT.get_child_from_name_and_label(gc, f"Coordinate{suffix}", "DataArray_t")
        PT.set_value(node, arr)

    ## Cells
    ## BIG ASSUMPTION HERE: same element type for all ghosts. It's not impossible to deal with mixed elements, but would first require a sort.
    ## However, note that most likely a part-to-part would be necessary to identify the ghost "element-type" because here we only have NGON
    ## infos.
    n_ghost_faces = len(data_ext["face_ln_to_gn"])
    face_edge_ext_idx = data_ext["face_edge_idx"]
    face_edge_ext = data_ext["face_edge"]
    total_edge_vtx = data_total["edge_vtx"]
    if n_ghost_faces == 0:
        return
    n_edges_by_cell = face_edge_ext_idx[1]
    if n_edges_by_cell == 4: # QUAD_4
        eltype = "QUAD_4"
        n_vtx_by_face = 4
        n_vtx_by_edge = 2
        # In order to call `combine_to_hexa`, we need to build
        # * np_edge_vtx_n : number of vertex per face (here, 4)
        # * np_edge_vtx : flattened connectivity face->vtx. It must describe
        #   exactly the faces that are "converted" to classical elements, it cannot
        #   be a "bigger" set of faces. In other words, element 1 is described by the 6 first
        #   faces of this array, element 2 by the 6 following etc
        # * np_face_edge : face orientation (so only the sign is important)
        econn = np.zeros(n_ghost_faces*n_vtx_by_face, dtype = face_edge_ext.dtype)
        np_edge_vtx_n = n_vtx_by_edge*np.ones(n_ghost_faces*n_edges_by_cell, dtype = face_edge_ext_idx.dtype)
        total_edge_vtx_reshaped = np.reshape(total_edge_vtx, (-1, 4))
        np_edge_vtx = total_edge_vtx_reshaped[np.abs(face_edge_ext) - 1, :].flatten() - 1 # 0-based for python
        np_face_edge = face_edge_ext # Only ghost cells, (we just need the sign so no need for 0-based)

        raise ValueError("in maia in 2D, it seems that connectivity is built from the array `edge_vtx`")
        econn += 1
    else:
        raise ValueError(f"n_faces_by_cell = {n_edges_by_cell} not implemented yet")

    # Search the range
    ranges = PT.Zone.get_elt_range_per_dim(zone)
    max_range = max(r[1] for r in ranges)
    erange = [max_range + 1, max_range + n_ghost_faces]

    PT.new_Elements('GHOST_'+eltype, type=eltype, erange=erange, econn=econn, parent=zone)

    # Fix Zone infos
    PT.set_value(zone, [[len(x), len(data_total["face_ln_to_gn"]), 0]])

def compute_cell_to_rank(cell_ln_to_gn, cell_ext_ln_to_gn, comm):
    "Compute the cell->rank connectivity, i.e. identify the rank owning each cell."
    rank = comm.Get_rank()
    # ln_to_gn = np.concatenate((cell_ln_to_gn, cell_ext_ln_to_gn))
    # cell_to_rank = rank * np.ones(len(ln_to_gn))
    cl_to_cg = cell_ln_to_gn
    pl_to_cg = cell_ext_ln_to_gn
    pl_to_pg = cell_ext_ln_to_gn

    g_gnum1 = [pl_to_pg]
    g_gnum2 = [cl_to_cg]
    part1_to_part2 = [pl_to_cg]
    part1_to_part2_idx = [np.arange(len(g) + 1).astype(np.intc) for g in g_gnum1]
    ptp = PDM.PartToPart(comm, g_gnum1, g_gnum2, part1_to_part2_idx, part1_to_part2)

    part2_data = [rank * np.ones(len(cell_ln_to_gn), dtype=np.intc)]
    request = ptp.reverse_iexch(
        PDM._PDM_MPI_COMM_KIND_P2P,
        PDM._PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
        # PDM._PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM ,
        part2_data,
    )
    part1_stride, part1_data = ptp.reverse_wait(request)
    print_dbg(["", f"[{rank}] part1_data", part1_data], comm)

    cell_to_rank = np.concatenate((part2_data[0], part1_data[0]))
    return cell_to_rank


def extend_partition_3d(part_tree, comm):
    """
    Tree is modified in-place.
    """
    rank = comm.Get_rank()

    # Retrieve connectivities from tree
    ppart = PT.get_node_from_name_and_label(part_tree, ":CGNS#Ppart", "UserDefinedData_t")
    pdm_output = {PT.get_name(n): PT.get_value(n) for n in  PT.get_nodes_from_label(ppart, "DataArray_t")}

    # Create part extension object
    extend_type = PDM.PartExtension.VTX
    depth = 1
    n_domain = 1
    n_part = 1
    part_ext = PDM.PartExtension(n_domain,
                                np.array([n_part]).astype(np.intc),
                                extend_type,
                                depth,
                                comm)

    # Set connectivities and coords
    i_domain = 0
    i_part = 0

    cell_face_idx = pdm_output["np_cell_face_idx"]
    cell_face = pdm_output["np_cell_face"]
    part_ext.connectivity_set(i_domain, i_part, PDM._PDM_CONNECTIVITY_TYPE_CELL_FACE, cell_face_idx, cell_face)

    face_vtx_idx = pdm_output["np_face_vtx_idx"]
    face_vtx = pdm_output["np_face_vtx"]
    part_ext.connectivity_set(i_domain, i_part, PDM._PDM_CONNECTIVITY_TYPE_FACE_VTX, face_vtx_idx, face_vtx)

    vtx_coord = pdm_output["np_vtx_coord"]
    part_ext.vtx_coord_set(i_domain, i_part, vtx_coord)

    # Set numberings
    cell_ln_to_gn = pdm_output["np_cell_ln_to_gn"]
    part_ext.ln_to_gn_set(i_domain, i_part, PDM._PDM_MESH_ENTITY_CELL, cell_ln_to_gn)

    face_ln_to_gn = pdm_output["np_face_ln_to_gn"]
    part_ext.ln_to_gn_set(i_domain, i_part, PDM._PDM_MESH_ENTITY_FACE, face_ln_to_gn)

    vtx_ln_to_gn = pdm_output["np_vtx_ln_to_gn"]
    part_ext.ln_to_gn_set(i_domain, i_part, PDM._PDM_MESH_ENTITY_VTX,  vtx_ln_to_gn)

    # Bound graph
    vtx_part_bound_proc_idx = pdm_output["np_vtx_part_bound_proc_idx"]
    vtx_part_bound_part_idx = pdm_output["np_vtx_part_bound_part_idx"]
    vtx_part_bound = pdm_output["np_vtx_part_bound"]
    part_ext.part_bound_graph_set(i_domain,
                                i_part,
                                PDM._PDM_MESH_ENTITY_VTX,
                                vtx_part_bound_proc_idx,
                                vtx_part_bound_part_idx,
                                vtx_part_bound)

    # Compute
    part_ext.compute()

    # Post process result

    ## Cell

    ### `cell_ext_ln_to_gn` corresponds to ghost_cell_local_to_global (1-based)
    cell_ext_ln_to_gn = part_ext.ln_to_gn_get(i_domain, i_part, PDM._PDM_MESH_ENTITY_CELL)

    ### `cell_face_ext_idx` + `cell_face_ext` is a sparse representation of the connectivity ghost-cell-to-face
    ### `cell_face_ext` is a flattened list of faces (signed to indicate orientation). `cell_face_ext_idx` indicates
    ### for each ghost cell `icell`, the indices to pick in `cell_face_ext` : faces indices of ghost cell `icell` are
    ### `[cell_face_ext_idx[icell], cell_face_ext_idx[icell+1]]`.
    cell_face_ext_idx, cell_face_ext = part_ext.connectivity_get(i_domain, i_part, PDM._PDM_CONNECTIVITY_TYPE_CELL_FACE)

    ## Face

    ### `face_ext_ln_to_gn` corresponds to ghost_face_local_to_global (1-based)
    face_ext_ln_to_gn = part_ext.ln_to_gn_get(i_domain, i_part, PDM._PDM_MESH_ENTITY_FACE)

    ### `face_vtx_ext_idx` + `face_vtx_ext` is a sparse representation of the connectivity ghost-face-to-vtx, see previous comment.
    ### Warning: not all faces from the previous `cell_face_ext` connectivity are listed here, but only ghost faces.
    face_vtx_ext_idx, face_vtx_ext = part_ext.connectivity_get(i_domain, i_part, PDM._PDM_CONNECTIVITY_TYPE_FACE_VTX)

    ## Vertices

    ### `vtx_ext_ln_to_gn` corresponds to ghost_vtx_local_to_global (1-based)
    vtx_ext_ln_to_gn = part_ext.ln_to_gn_get(i_domain, i_part, PDM._PDM_MESH_ENTITY_VTX)

    ### `vtx_coord_ext` is the coordinates of "ghost" nodes, i.e nodes used in the previous connectivities
    ### Arranged as (x1,y1,z1, x2,y2,z2, x3,y3,z3, ...)
    vtx_coord_ext = part_ext.vtx_coord_get(i_domain, i_part)

    # Print
    print_dbg(
        [
            "", f"[{rank}] Infos",
            f"vtx_ln_to_gn",
            vtx_ln_to_gn,
            f"cell_face",
            cell_face,
            f"face_vtx l={len(face_vtx)}",
            face_vtx,
            f"cell_ext_ln_to_gn",
            cell_ext_ln_to_gn,
            f"cell_face_ext_idx",
            cell_face_ext_idx,
            f"cell_face_ext",
            cell_face_ext,
            f"face_vtx_ext_idx",
            face_vtx_ext_idx,
            f"face_vtx_ext",
            face_vtx_ext,
            f"vtx_ext_ln_to_gn",
            vtx_ext_ln_to_gn,
            f"vtx_coord_ext",
            vtx_coord_ext
        ],
    comm)

    # Append extension to CGNS
    zone = PT.get_node_from_label(part_tree, "Zone_t")
    data = {}
    data["cell_ln_to_gn"] = cell_ln_to_gn
    data["face_ln_to_gn"] = face_ln_to_gn
    data["vtx_ln_to_gn"] = vtx_ln_to_gn
    data["cell_face_idx"] = cell_face_idx
    data["cell_face"] = cell_face
    data["face_vtx_idx"] = face_vtx_idx
    data["face_vtx"] = face_vtx
    data["vtx_coord"] = vtx_coord
    data_ext = {}
    data_ext["cell_ln_to_gn"] = cell_ext_ln_to_gn
    data_ext["face_ln_to_gn"] = face_ext_ln_to_gn
    data_ext["vtx_ln_to_gn"] = vtx_ext_ln_to_gn
    data_ext["cell_face_idx"] = cell_face_ext_idx
    data_ext["cell_face"] = cell_face_ext
    data_ext["face_vtx_idx"] = face_vtx_ext_idx
    data_ext["face_vtx"] = face_vtx_ext
    data_ext["vtx_coord"] = vtx_coord_ext

    # Concatenate
    data_tot = concatenate_connectivities_3d(data, data_ext)

    # Compute cell_to_rank and append it to
    cell_to_rank = compute_cell_to_rank(cell_ln_to_gn, cell_ext_ln_to_gn, comm)

    # Append to cgns
    append_extension_to_tree_3d(zone, data_ext, data_tot, cell_to_rank)

    write_tree(part_tree, f"part{rank}.cgns")


def extend_partition_2d(part_tree, comm):
    """
    Tree is modified in-place. Doesn't work yet. Implementation is inspired from sonics:cgns_to_partition_extension.py
    """
    rank = comm.Get_rank()

    # Retrieve connectivities from tree
    ppart = PT.get_node_from_name_and_label(part_tree, ":CGNS#Ppart", "UserDefinedData_t")
    pdm_output = {PT.get_name(n): PT.get_value(n) for n in  PT.get_nodes_from_label(ppart, "DataArray_t")}

    # Create part extension object
    extend_type = PDM.PartExtension.VTX
    depth = 1
    n_domain = 1
    n_part = 1
    part_ext = PDM.PartExtension(n_domain,
                                np.array([n_part]).astype(np.intc),
                                extend_type,
                                depth,
                                comm)

    # Set connectivities and coords
    i_domain = 0
    i_part = 0

    face_edge_idx = pdm_output["np_face_edge_idx"]
    face_edge = pdm_output["np_face_edge"]
    part_ext.connectivity_set(i_domain, i_part, PDM._PDM_CONNECTIVITY_TYPE_FACE_EDGE, face_edge_idx, face_edge)

    edge_vtx = pdm_output["np_edge_vtx"]
    part_ext.connectivity_set(i_domain, i_part, PDM._PDM_CONNECTIVITY_TYPE_EDGE_VTX , None, edge_vtx)

    face_vtx_idx = pdm_output["np_face_vtx_idx"]
    face_vtx = pdm_output["np_face_vtx"]
    part_ext.connectivity_set(i_domain, i_part, PDM._PDM_CONNECTIVITY_TYPE_FACE_VTX, face_vtx_idx, face_vtx)

    vtx_coord = pdm_output["np_vtx_coord"]
    part_ext.vtx_coord_set(i_domain, i_part, vtx_coord)

    # Set numberings
    face_ln_to_gn = pdm_output["np_face_ln_to_gn"]
    part_ext.ln_to_gn_set(i_domain, i_part, PDM._PDM_MESH_ENTITY_FACE, face_ln_to_gn)

    edge_ln_to_gn = pdm_output["np_edge_ln_to_gn"]
    part_ext.ln_to_gn_set(i_domain, i_part, PDM._PDM_MESH_ENTITY_EDGE, edge_ln_to_gn)

    vtx_ln_to_gn = pdm_output["np_vtx_ln_to_gn"]
    part_ext.ln_to_gn_set(i_domain, i_part, PDM._PDM_MESH_ENTITY_VTX, vtx_ln_to_gn)

    # Bound graph
    vtx_part_bound_proc_idx = pdm_output["np_vtx_part_bound_proc_idx"]
    vtx_part_bound_part_idx = pdm_output["np_vtx_part_bound_part_idx"]
    vtx_part_bound = pdm_output["np_vtx_part_bound"]
    part_ext.part_bound_graph_set(i_domain,
                                i_part,
                                PDM._PDM_MESH_ENTITY_VTX,
                                vtx_part_bound_proc_idx,
                                vtx_part_bound_part_idx,
                                vtx_part_bound)

    print_dbg(["passage 1"])

    print_dbg(
        [
            "", f"[{rank}] Infos",
            f"vtx_ln_to_gn",
            vtx_ln_to_gn,
            f"face_ln_to_gn",
            face_ln_to_gn,
            f"face_vtx_idx",
            face_vtx_idx,
            f"face_vtx",
            face_vtx,
            f"vtx_ln_to_gn",
            vtx_ln_to_gn,
            f"vtx_coord",
            vtx_coord
        ],
    comm)

    # Compute
    part_ext.compute()

    print_dbg(["passage 2"])

    # Post process result, for comments see the 3d version

    ## Face

    face_ext_ln_to_gn = part_ext.ln_to_gn_get(i_domain,
                                            i_part,
                                            PDM._PDM_MESH_ENTITY_FACE)

    face_vtx_ext_idx, face_vtx_ext = part_ext.connectivity_get(i_domain,
                                                                i_part,
                                                                PDM._PDM_CONNECTIVITY_TYPE_FACE_VTX)

    ## Vertices

    ### `vtx_ext_ln_to_gn` corresponds to ghost_vtx_local_to_global (1-based)
    vtx_ext_ln_to_gn = part_ext.ln_to_gn_get(i_domain,
                                            i_part,
                                            PDM._PDM_MESH_ENTITY_VTX)

    ### `vtx_coord_ext` is the coordinates of "ghost" nodes, i.e nodes used in the previous connectivities
    ### Arranged as (x1,y1,z1, x2,y2,z2, x3,y3,z3, ...)
    vtx_coord_ext = part_ext.vtx_coord_get(i_domain,
                                        i_part)

    # Print
    print_dbg(
        [
            "", f"[{rank}] Infos",
            f"vtx_ln_to_gn",
            vtx_ln_to_gn,
            f"face_ext_ln_to_gn",
            face_ext_ln_to_gn,
            f"face_vtx_ext_idx",
            face_vtx_ext_idx,
            f"face_vtx_ext",
            face_vtx_ext,
            f"vtx_ext_ln_to_gn",
            vtx_ext_ln_to_gn,
            f"vtx_coord_ext",
            vtx_coord_ext
        ],
    comm)
    raise ValueError("dbg maia ")

    # Append extension to CGNS
    zone = PT.get_node_from_label(part_tree, "Zone_t")
    data = {}
    data["face_ln_to_gn"] = face_ln_to_gn
    data["edge_ln_to_gn"] = edge_ln_to_gn
    data["vtx_ln_to_gn"] = vtx_ln_to_gn
    data["face_edge_idx"] = face_edge_idx
    data["face_edge"] = face_edge
    data["edge_vtx_idx"] = edge_vtx_idx
    data["edge_vtx"] = edge_vtx
    data["vtx_coord"] = vtx_coord
    data_ext = {}
    data_ext["face_ln_to_gn"] = face_ext_ln_to_gn
    data_ext["edge_ln_to_gn"] = edge_ext_ln_to_gn
    data_ext["vtx_ln_to_gn"] = vtx_ext_ln_to_gn
    data_ext["face_edge_idx"] = face_edge_ext_idx
    data_ext["face_edge"] = face_edge_ext
    data_ext["edge_vtx_idx"] = edge_vtx_ext_idx
    data_ext["edge_vtx"] = edge_vtx_ext
    data_ext["vtx_coord"] = vtx_coord_ext

    data_tot = concatenate_connectivities_2d(data, data_ext)
    append_extension_to_tree_2d(zone, data_ext, data_tot)

    write_tree(part_tree, f"part{rank}.cgns")

    # Compute cell_to_rank
    cell_to_rank = compute_cell_to_rank(face_ln_to_gn, face_ext_ln_to_gn, comm)
    PT.new_DataArray(name="cell_to_rank", value=cell_to_rank, parent=ppart)

def extend_partition_2d_backup(part_tree, comm):
    """
    Tree is modified in-place.
    """
    rank = comm.Get_rank()

    # Retrieve connectivities from tree
    ppart = PT.get_node_from_name_and_label(part_tree, ":CGNS#Ppart", "UserDefinedData_t")
    pdm_output = {PT.get_name(n): PT.get_value(n) for n in  PT.get_nodes_from_label(ppart, "DataArray_t")}

    # Create part extension object
    extend_type = PDM.PartExtension.VTX
    depth = 1
    n_domain = 1
    n_part = 1
    part_ext = PDM.PartExtension(n_domain,
                                np.array([n_part]).astype(np.intc),
                                extend_type,
                                depth,
                                comm)

    # Set connectivities and coords
    i_domain = 0
    i_part = 0

    face_edge_idx = pdm_output["np_face_edge_idx"]
    face_edge = pdm_output["np_face_edge"]
    part_ext.connectivity_set(i_domain,
                            i_part,
                            PDM._PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                            face_edge_idx,
                            face_edge)

    edge_vtx_idx = pdm_output["np_edge_vtx_idx"]
    edge_vtx = pdm_output["np_edge_vtx"]
    part_ext.connectivity_set(i_domain,
                            i_part,
                            PDM._PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                            edge_vtx_idx,
                            edge_vtx)

    vtx_coord = pdm_output["np_vtx_coord"]
    part_ext.vtx_coord_set(i_domain,
                        i_part,
                        vtx_coord)

    # Set numberings
    face_ln_to_gn = pdm_output["np_face_ln_to_gn"]
    part_ext.ln_to_gn_set(i_domain,
                        i_part,
                        PDM._PDM_MESH_ENTITY_FACE,
                        face_ln_to_gn)

    edge_ln_to_gn = pdm_output["np_edge_ln_to_gn"]
    part_ext.ln_to_gn_set(i_domain,
                        i_part,
                        PDM._PDM_MESH_ENTITY_EDGE,
                        edge_ln_to_gn)

    vtx_ln_to_gn = pdm_output["np_vtx_ln_to_gn"]
    part_ext.ln_to_gn_set(i_domain,
                        i_part,
                        PDM._PDM_MESH_ENTITY_VTX,
                        vtx_ln_to_gn)

    # Bound graph
    vtx_part_bound_proc_idx = pdm_output["np_vtx_part_bound_proc_idx"]
    vtx_part_bound_part_idx = pdm_output["np_vtx_part_bound_part_idx"]
    vtx_part_bound = pdm_output["np_vtx_part_bound"]
    part_ext.part_bound_graph_set(i_domain,
                                i_part,
                                PDM._PDM_MESH_ENTITY_VTX,
                                vtx_part_bound_proc_idx,
                                vtx_part_bound_part_idx,
                                vtx_part_bound)

    # Compute
    part_ext.compute()

    # Post process result, for comments see the 3d version

    ## Face

    face_ext_ln_to_gn = part_ext.ln_to_gn_get(i_domain,
                                            i_part,
                                            PDM._PDM_MESH_ENTITY_FACE)

    face_edge_ext_idx, face_edge_ext = part_ext.connectivity_get(i_domain,
                                                                i_part,
                                                                PDM._PDM_CONNECTIVITY_TYPE_FACE_EDGE)
    face_vtx_ext_idx, face_vtx_ext = part_ext.connectivity_get(i_domain,
                                                                i_part,
                                                                PDM._PDM_CONNECTIVITY_TYPE_FACE_VTX)

    ## Edge

    edge_ext_ln_to_gn = part_ext.ln_to_gn_get(i_domain,
                                            i_part,
                                            PDM._PDM_MESH_ENTITY_EDGE)

    edge_vtx_ext_idx, edge_vtx_ext = part_ext.connectivity_get(i_domain,
                                                            i_part,
                                                            PDM._PDM_CONNECTIVITY_TYPE_EDGE_VTX)

    ## Vertices

    ### `vtx_ext_ln_to_gn` corresponds to ghost_vtx_local_to_global (1-based)
    vtx_ext_ln_to_gn = part_ext.ln_to_gn_get(i_domain,
                                            i_part,
                                            PDM._PDM_MESH_ENTITY_VTX)

    ### `vtx_coord_ext` is the coordinates of "ghost" nodes, i.e nodes used in the previous connectivities
    ### Arranged as (x1,y1,z1, x2,y2,z2, x3,y3,z3, ...)
    vtx_coord_ext = part_ext.vtx_coord_get(i_domain,
                                        i_part)

    # Print
    print_dbg(
        [
            "", f"[{rank}] Infos",
            f"vtx_ln_to_gn",
            vtx_ln_to_gn,
            f"face_edge",
            face_edge,
            f"edge_vtx l={len(edge_vtx)}",
            edge_vtx,
            f"face_ext_ln_to_gn",
            face_ext_ln_to_gn,
            f"face_edge_ext_idx",
            face_edge_ext_idx,
            f"face_edge_ext",
            face_edge_ext,
            f"face_vtx_ext_idx",
            face_vtx_ext_idx,
            f"face_edge_ext",
            face_edge_ext,
            f"edge_vtx_ext_idx",
            edge_vtx_ext_idx,
            f"edge_vtx_ext",
            edge_vtx_ext,
            f"vtx_ext_ln_to_gn",
            vtx_ext_ln_to_gn,
            f"vtx_coord_ext",
            vtx_coord_ext
        ],
    comm)
    raise ValueError("dbg maia ")

    # Append extension to CGNS
    zone = PT.get_node_from_label(part_tree, "Zone_t")
    data = {}
    data["face_ln_to_gn"] = face_ln_to_gn
    data["edge_ln_to_gn"] = edge_ln_to_gn
    data["vtx_ln_to_gn"] = vtx_ln_to_gn
    data["face_edge_idx"] = face_edge_idx
    data["face_edge"] = face_edge
    data["edge_vtx_idx"] = edge_vtx_idx
    data["edge_vtx"] = edge_vtx
    data["vtx_coord"] = vtx_coord
    data_ext = {}
    data_ext["face_ln_to_gn"] = face_ext_ln_to_gn
    data_ext["edge_ln_to_gn"] = edge_ext_ln_to_gn
    data_ext["vtx_ln_to_gn"] = vtx_ext_ln_to_gn
    data_ext["face_edge_idx"] = face_edge_ext_idx
    data_ext["face_edge"] = face_edge_ext
    data_ext["edge_vtx_idx"] = edge_vtx_ext_idx
    data_ext["edge_vtx"] = edge_vtx_ext
    data_ext["vtx_coord"] = vtx_coord_ext

    data_tot = concatenate_connectivities_2d(data, data_ext)
    append_extension_to_tree_2d(zone, data_ext, data_tot)

    write_tree(part_tree, f"part{rank}.cgns")

    # Compute cell_to_rank
    cell_to_rank = compute_cell_to_rank(face_ln_to_gn, face_ext_ln_to_gn, comm)
    PT.new_DataArray(name="cell_to_rank", value=cell_to_rank, parent=ppart)
