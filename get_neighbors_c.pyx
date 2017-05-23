# distutils: language = c++
# distutils: sources = neighbor/neighbor/neighbor.cpp

from libcpp.vector cimport vector
import numpy as np
cimport numpy as np

cdef extern from 'neighbor/neighbor/neighbor.h':
    vector[vector[int]] simple_neighbors(vector[vector[double]] pos, double cutoff)
    vector[vector[int]] cell_list_neighbors(vector[vector[double]] pos, double cutoff)

def simple_neighbors_c(np.ndarray pos, double cutoff):
    """C++ wrapper for simple neighbors algorithm"""
    # Convert 2-D numpy array pos into vector of vector of int pos_c
    N = pos.shape[0]
    cdef vector[vector[double]] pos_c
    for i in range(N):
        pos_c.push_back(pos[i, :])

    # Run C++ function
    cdef vector[vector[int]] neighbors_c = simple_neighbors(pos_c, cutoff)

    # Convert vector of vector of int neighbors_c into list of 1-D numpy arrays neighbors
    neighbors = []
    for i in range(N):
        neighbors.append(np.array(neighbors_c[i], dtype=np.int64))

    return neighbors


def cell_list_neighbors_c(np.ndarray pos, double cutoff):
    """C++ wrapper for cell list neighbors algorithm"""
    # Convert 2-D numpy array pos into vector of vector of int pos_c
    N = pos.shape[0]
    cdef vector[vector[double]] pos_c
    for i in range(N):
        pos_c.push_back(pos[i, :])

    # Run C++ function
    cdef vector[vector[int]] neighbors_c = cell_list_neighbors(pos_c, cutoff)

    # Convert vector of vector of int neighbors_c into list of 1-D numpy arrays neighbors
    neighbors = []
    for i in range(N):
        neighbors.append(np.array(neighbors_c[i], dtype=np.int64))

    return neighbors