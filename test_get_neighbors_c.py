# Quick and dirty tests of Cython get neighbors

from get_neighbors_c import simple_neighbors_c, cell_list_neighbors_c
import random
import time
from utils import gen_random_atoms

if __name__ == '__main__':
    random.seed(0)
    N = 1000  # cell list becomes faster than simple somewhere between 1000 and 5000
    names, pos = gen_random_atoms(N)

    cutoff = 0.5

    # Simple neighbors calc
    start = time.time()
    neighbors = simple_neighbors_c(pos, cutoff)
    end = time.time()
    print(str(end - start) + 's for simple neighbors using C++ function')

    # Cell list neighbors calc
    start = time.time()
    neighbors = cell_list_neighbors_c(pos, cutoff)
    end = time.time()
    print(str(end - start) + 's for cell list neighbors using C++ function')

    # print(neighbors)

    1
