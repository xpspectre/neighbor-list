import numpy as np
import random


def gen_random_atoms(N):
    """Generate test set of atoms. Output N x 1 list of strings of atom names, N x 3 numpy array of x,y,z coords.
    Draw from uniformly distributed random positions in 3-D space within [-1,1] in each dimension."""
    lb = -1
    ub = 1
    names = []
    pos = np.zeros((N,3), dtype=np.float64)
    for i in range(N):
        names.append('A{i}'.format(i=i))
        pos[i, 0] = random.uniform(lb, ub)
        pos[i, 1] = random.uniform(lb, ub)
        pos[i, 2] = random.uniform(lb, ub)

    return names, pos


def dump_atoms(pos, filename):
    """Dump atom coordinates to text file, each line as 3 coords, for testing C implementation"""
    with open(filename, 'w') as f:
        for atom in pos:
            f.write('{}\t{}\t{}\n'.format(atom[0], atom[1], atom[2]))

if __name__ == '__main__':
    N = 20
    names, pos = gen_random_atoms(N)
