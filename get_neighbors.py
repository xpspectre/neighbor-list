import random
import numpy as np
from utils import gen_random_atoms, dump_atoms
from sklearn.metrics.pairwise import paired_distances
from scipy.spatial.distance import pdist, squareform, cdist
import itertools
import time

neighbor_mask = np.array(list(itertools.product([-1, 0, 1], [-1, 0, 1], [-1, 0, 1])))  # do this once for speed


def get_neighbors(cell, bounds):
    """Get list of neighbor cells (27 including cell).
    bounds is a 2-tuble of 3-tuples, 1st dim = (lo, hi), 2nd dim = (x,y,z)."""
    neighbors = cell + neighbor_mask
    keep = np.logical_and(np.greater_equal(neighbors, bounds[0]), np.less_equal(neighbors, bounds[1]))  # broadcasts so each of the 3 elements of bounds[i] goes down each col of neighbors
    keep = np.all(keep, axis=1)
    return neighbors[keep, :]


def simple_neighbors(pos, cutoff, max_neighbors=0):
    """O(N^2) get list of each atom's neighbors within cutoff using naive method. Mainly for comparison.
    N is the total number of atoms. Neighbor list returned as list of positions in the original pos array that are
    within cutoff.
    Use max_neighbors to specify the max_neighbors closest atoms. All are within cutoff; if there aren't enough, the
        list will be short. max_neighbors=0 indicates don't sort and take all neighbors within cutoff."""
    N = pos.shape[0]
    n = np.arange(N, dtype=np.int32)

    # Calculate distances of all atoms to all atoms (actually just the lower triangle) and then filter
    # O(N^2) in both time and space but simple
    dists = squareform(pdist(pos))

    neighbors = []
    for i in range(N):
        atom_i_dists = dists[i, :]
        within_cutoff = atom_i_dists < cutoff
        atom_i_neighbors = np.flatnonzero(within_cutoff)

        # Limit to max neighbors if specified; drop self-distance
        if max_neighbors > 0:
            atom_i_dists = atom_i_dists[within_cutoff]
            sort_ind = np.argsort(atom_i_dists)
            atom_i_neighbors = atom_i_neighbors[sort_ind][1:max_neighbors+1]
        else:
            atom_i_neighbors = atom_i_neighbors[atom_i_neighbors != i]

        neighbors.append(atom_i_neighbors)

    return neighbors


def cell_list_neighbors(pos, cutoff, max_neighbors=0):
    """O(N) get list of each atom's neighbors within cutoff using Yip and Elber (1989) cell list algorithm.
    N is the total number of atoms. Neighbor list returned as list of positions in the original pos array that are
    within cutoff.
    Use max_neighbors to specify the max_neighbors closest atoms. All are within cutoff; if there aren't enough, the
        list will be short. max_neighbors=0 indicates don't sort and take all neighbors within cutoff."""
    N = pos.shape[0]

    # Set cells to enclose all atoms
    lb = np.amin(pos, axis=0)
    ub = np.amax(pos, axis=0)
    size = ub - lb
    # Add in a little slack to the edges
    slack_mult = 1.01
    lb -= (slack_mult - 1.0) * size
    ub += (slack_mult - 1.0) * size

    # Set cell size slightly larger than cutoff
    #   This guarantees only neighboring cells need to be searched
    cutoff_mult = 1.01
    cell_size = cutoff * cutoff_mult

    # Divide space into cells, leaving extra space on the upper side if needed
    edges = []  # cell boundaries in each dim
    bounds = []  # lo and hi index of each dim's cells, inclusive
    for i in range(3):
        bounds_i = np.arange(lb[i], ub[i] + cell_size, cell_size)
        edges.append(bounds_i)
        bounds.append((0, bounds_i.size - 2))
    bounds = tuple(zip(*bounds))  # turn into form expected by get_neighbors

    # Distribute atoms to cells
    cell_pos = np.zeros((N,3), dtype=np.int64)
    for i in range(3):
        cell_pos[:, i] = np.digitize(pos[:, i], edges[i])  # in smallest bin is index 1
    cell_pos -= 1  # in smallest bin is index 0

    # Store as cell:atoms dict
    cell_atom_map = {}
    cell_pos_ = map(tuple, cell_pos)
    for i, cell in enumerate(cell_pos_):
        if cell in cell_atom_map:
            cell_atom_map[cell].append(i)
        else:
            cell_atom_map[cell] = [i]

    # Calculate neighbors cell-by-cell
    neighbors = {}
    for cell, atoms in cell_atom_map.items():
        neighbor_cells = list(map(tuple, get_neighbors(cell, bounds)))
        neighbor_atoms = set()
        for neighbor_cell in neighbor_cells:
            if neighbor_cell in cell_atom_map:
                neighbor_atoms.update(cell_atom_map[neighbor_cell])
        neighbor_atoms = np.sort(np.array(list(neighbor_atoms), dtype=np.int64))

        neighbor_pos = pos[neighbor_atoms, :]
        atom_pos = pos[atoms, :]
        dists = cdist(atom_pos, neighbor_pos)
        for i, atom in enumerate(atoms):
            atom_i_dists = dists[i, :]
            within_cutoff = atom_i_dists < cutoff
            atom_i_neighbors = neighbor_atoms[within_cutoff]

            # Limit to max neighbors if specified; drop self-distance
            if max_neighbors > 0:
                atom_i_dists = atom_i_dists[within_cutoff]
                sort_ind = np.argsort(atom_i_dists)
                atom_i_neighbors = atom_i_neighbors[sort_ind][1:max_neighbors + 1]
            else:
                atom_i_neighbors = atom_i_neighbors[atom_i_neighbors != atom]

            neighbors[atom] = atom_i_neighbors

    # Convert to list by original index
    neighbors_ = []
    for i in range(N):
        neighbors_.append(neighbors[i])

    return neighbors_


if __name__ == '__main__':
    random.seed(0)
    N = 1000  # cell list becomes faster than simple somewhere between 1000 and 5000
    names, pos = gen_random_atoms(N)

    dump_atoms(pos, 'test_atoms.txt')

    cutoff = 0.5

    # Simple neighbors calc
    start = time.time()
    neighbors = simple_neighbors(pos, cutoff, max_neighbors=12)
    end = time.time()
    print(str(end - start) + 's for simple neighbors')

    # Cell list neighbors calc
    start = time.time()
    neighbors_2 = cell_list_neighbors(pos, cutoff, max_neighbors=12)
    end = time.time()
    print(str(end - start) + 's for cell list neighbors')

    # Compare results
    assert len(neighbors) == len(neighbors_2)
    for i in range(N):
        assert np.all(neighbors[i] == neighbors_2[i])

    # x = np.array([1,2,3])
    # y = np.array([1,2,3])
    # p = cartesian_prod(x, y)

    # cell = [1,1,1]
    # cell = [0, 1, 1]
    # # bounds = ((0,3), (0,3), (0,3))
    # bounds = ((0,0,0), (3,3,1))
    # neighbors = get_neighbors(cell, bounds)

    1
