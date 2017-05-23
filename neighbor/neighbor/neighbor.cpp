#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>

using namespace std;

vector<vector<int> > simple_neighbors(vector<vector<double> > pos, double cutoff) {
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

int main() {
    // Load test file of coordinates
    ifstream infile("../../test_atoms.txt");
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
    cout << "Loaded atom coordinates:\n";
    for (vector<vector<double> >::size_type i = 0; i != pos.size(); ++i) {
        cout << pos[i][0] << "\t" << pos[i][1] << "\t" << pos[i][2] << endl;
    }

    
    double cutoff = 0.5;
    vector<vector<int> > neighbors = simple_neighbors(pos, cutoff);

    // Print out neighbor list
    cout << "Calculated neighbor list:\n";
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
