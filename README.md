# neighbor-list

Algorithms for finding neighboring atoms within some cutoff. Useful enough that I put it in a separate repo.

Used Python 3. Requires numpy and scipy.

## C++ implementation

The `neighbor` subdir contains a MSVC project with fast, C++ implementations of the algorithms.

To compile using GCC, run `g++ -std=c++11 neighbor.cpp -o <whatever>`

To Cythonize, run `python setup.py build_ext -i`.

Note: It doesn't look like the C++ version is faster right now.