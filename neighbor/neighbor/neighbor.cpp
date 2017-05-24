#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <math.h>
#include <ctime>
#include "neighbor.h"

using namespace std;

vector<double> arange(double start, double stop, double step) {
	// Build range of values starting at stop, incrementing by step, and stopping after the last value is greater than or equal to stop (can overhang)
	vector<double> range;
	range.push_back(start);
	double val = start; // accumulates next val
	while (true) {
		val += step;
		range.push_back(val);
		if (val >= stop) {
			return range;
		}
	}
}

int assign_bin(double val, const vector<double>& bin_edges) {
	// Assign val to histogram bin, returning bin index. Leftmost bin is 0. Rightmost bin is n_bin_edges-2.
	// Assumes bin_edge are sorted. Assumes val falls within bin_edges.
	// Note: A binary search would be more appropriate if this were operating on larger bin_edges vectors, but linear search is OK for this.
	for (int i = 1; i < bin_edges.size(); ++i) {  // start at 2nd edge since assume val > leftmost bin_edge
		if (val < bin_edges[i]) {
			return i - 1;
		}
	}
	// Should not fall thru to here
	return -1;
}

vector<vector<int> > get_neighbor_cells(const vector<int>& cell) {
	// Get list of neighboring cells. Don't worry about out-of-bound cells - these will be ignored in the next step.
	vector<vector<int> > neighbors;
	int i, j, k;
	for (i = -1; i <= 1; ++i) {
		for (j = -1; j <= 1; ++j) {
			for (k = -1; k <= 1; ++k) {
				vector<int> neighbor = cell;
				neighbor[0] += i;
				neighbor[1] += j;
				neighbor[2] += k;
				neighbors.push_back(neighbor);
			}
		}
	}
	return neighbors;
}

vector<vector<double> > cdist(const vector<vector<double> >& x, const vector<vector<double> >& y) {
	// Calculate distance between points in x and points in y, outputting a matrix where x is 1st dim, y is 2nd dim
	// Assumes 3-D points
	int nx = x.size();
	int ny = y.size();
	int i, j;
	vector<vector<double> > dists(nx, vector<double>(ny));
	for (i = 0; i < nx; ++i) {
		for (j = 0; j < ny; ++j) {
			dists[i][j] = sqrt(pow(x[i][0] - y[j][0], 2) + pow(x[i][1] - y[j][1], 2) + pow(x[i][2] - y[j][2], 2));
		}
	}
	return dists;
}

vector<vector<int> > simple_neighbors(const vector<vector<double> >& pos, double cutoff) {
    // O(N^2) get list of each atom's neighbors within cutoff using naive method
    int N = pos.size();

    // Calculate distances of all atoms to all atoms
    int i, j;
    vector<vector<double> > dists(N, vector<double>(N));
    for (i = 0; i < N; ++i) {
        for (j = 0; j < i; ++j) {
            double dist = sqrt(pow(pos[i][0] - pos[j][0], 2) + pow(pos[i][1] - pos[j][1], 2) + pow(pos[i][2] - pos[j][2], 2));
            dists[i][j] = dist;
            dists[j][i] = dist;
        }
    }

    // Filter for neighboring atoms within cutoff
    vector<vector<int> > neighbors;
    for (i = 0; i < N; ++i) {
        vector<double> atom_i_dists = dists[i];
        vector<int> atom_i_neighbors;
        for (j = 0; j < N; ++j) {
            if (j != i && atom_i_dists[j] < cutoff) {
                atom_i_neighbors.push_back(j);
            }
        }
        neighbors.push_back(atom_i_neighbors);
    }

    return neighbors;
}

vector<vector<int> > cell_list_neighbors(const vector<vector<double> >& pos, double cutoff) {
	// O(n) get list of each atom's neighbors within cutoff using Yip and Elber (1989) cell list algorithm.
	int i, j;
	int N = pos.size();

	// Setup cells to enclose all atoms
	// Get lower and upper bounds on cells
	vector<double> lb = { INFINITY, INFINITY, INFINITY };
	vector<double> ub = { -INFINITY, -INFINITY, -INFINITY };
	double val;
	for (i = 0; i < N; ++i) {
		for (j = 0; j < 3; ++j) {
			val = pos[i][j];
			if (val < lb[j]) {
				lb[j] = val;
			} else if (val > ub[j]) {
				ub[j] = val;
			}
		}
	}
	vector<double> size(3);
	for (j = 0; j < 3; ++j) {
		size[j] = ub[j] - lb[j];
	}

	// Add in a little slack to the edges
	double slack_mult = 1.01;
	for (j = 0; j < 3; ++j) {
		lb[j] -= (slack_mult - 1.0) * size[j];
		ub[j] += (slack_mult - 1.0) * size[j];
	}

	// Set cell size slightly larger than cutoff
	//	This guarantees only neighboring cells need to be searched
	double cutoff_mult = 1.01;
	double cell_size = cutoff * cutoff_mult;

	// Divide space into cells, leaving extra space on the upper side if needed
	vector<vector<double> > edges(3); // cell boundaries in each dim
	for (j = 0; j < 3; ++j) {
		edges[j] = arange(lb[j], ub[j] + cell_size, cell_size);;
	}

	// Distribute atoms to cells
	vector<vector<int> > cell_pos(N, vector<int>(3));
	for (i = 0; i < N; ++i) {
		for (j = 0; j < 3; ++j) {
			cell_pos[i][j] = assign_bin(pos[i][j], edges[j]);
		}
	}
	
	// Store cell:atoms map
	map<vector<int>, vector<int> > cell_atom_map; // cell pos in x,y,z : atom indices
	vector<int> key;
	for (i = 0; i < N; ++i) {
		key = cell_pos[i];
		if (cell_atom_map.count(key) == 0) { // cell not present, make a new list of atoms
			cell_atom_map[key] = vector<int>{i};
		} else { // present, add to list of atoms already in cell
			cell_atom_map[key].push_back(i);
		}
	}

	// Calculate neighbors cell by cell
	map<int, vector<int> > neighbors_map;
	vector<int> cell;
	vector<int> atoms;
	vector<int> neighbor_cell;
	vector<int> neighbor_atom;
	for (auto const& x : cell_atom_map) {
		cell = x.first;
		atoms = x.second;

		// Get indices of neighbors
		vector<vector<int> > neighbor_cells = get_neighbor_cells(cell);
		set<int> neighbor_atoms_set; // all atoms in cell and neighboring cells
		for (i = 0; i < neighbor_cells.size(); ++i) {
			neighbor_cell = neighbor_cells[i];
			if (cell_atom_map.count(neighbor_cell) == 1) {
				neighbor_atom = cell_atom_map[neighbor_cell];
				for (j = 0; j < neighbor_atom.size(); ++j) {
					neighbor_atoms_set.insert(neighbor_atom[j]);
				}
			}
		}
		vector<int> neighbor_atoms(neighbor_atoms_set.begin(), neighbor_atoms_set.end()); // copy into a vector; sorted by default

		// Get positions of cell atoms and neighbor atoms
		vector<vector<double> > atom_pos;
		vector<vector<double> > neighbor_pos;
		for (i = 0; i < atoms.size(); ++i) {
			atom_pos.push_back(pos[atoms[i]]);
		}
		for (i = 0; i < neighbor_atoms.size(); ++i) {
			neighbor_pos.push_back(pos[neighbor_atoms[i]]);
		}
		vector<vector<double> > dists = cdist(atom_pos, neighbor_pos);
		int atom;
		vector<double> atom_i_dists;
		for (i = 0; i < atoms.size(); ++i) {
			atom = atoms[i];
			atom_i_dists = dists[i];
			vector<int> atom_i_neighbors;
			for (j = 0; j < neighbor_atoms.size(); ++j) {
				if (atom_i_dists[j] < cutoff && atom != neighbor_atoms[j]) {
					atom_i_neighbors.push_back(neighbor_atoms[j]);
				}
			}
			neighbors_map[atom] = atom_i_neighbors;
		}
	}

	// Convert neighbors map into expected vector
	//	This should be complete
	vector<vector<int> > neighbors;
	for (j = 0; j < N; ++j) {
		neighbors.push_back(neighbors_map[j]);
	}
	
	return neighbors;
}

int main() {
    // Load test file of coordinates
    ifstream infile("C:/Users/xpspectre/workspace/neighbor-list/test_atoms.txt");
    string s;
    vector<vector<double> > pos;
    int i = 0;
    while (getline(infile, s)) {
        stringstream ss(s);
        vector<double> row;
        double num;
        while (ss >> num) {
            row.push_back(num);
        }
        pos.push_back(row);
        i++;
    }

    // Print out loaded coords
    /*cout << "Loaded atom coordinates:\n";
    for (vector<vector<double> >::size_type i = 0; i != pos.size(); ++i) {
        cout << pos[i][0] << "\t" << pos[i][1] << "\t" << pos[i][2] << endl;
    }*/

    
    double cutoff = 0.5;

	clock_t start;
	start = clock();
    //vector<vector<int> > neighbors = simple_neighbors(pos, cutoff);
	vector<vector<int> > neighbors = cell_list_neighbors(pos, cutoff);
	cout << "Time: " << (clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << endl;

    // Print out neighbor list
    cout << "calculated neighbor list:\n";
    for (i = 0; i != neighbors.size(); ++i) {
        vector<int> neighbors_i = neighbors[i];
        string line;
        line += to_string(i) + ":\t";
        for (int j = 0; j < neighbors_i.size(); ++j) {
            line += to_string(neighbors_i[j]) + " ";
        }
        line += "\n";
        cout << line;
    }

    return 0;
}
