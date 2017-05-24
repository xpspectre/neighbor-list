#pragma once

#include <vector>

std::vector<double> arange(double, double, double);
int assign_bin(double, const std::vector<double>&);
std::vector<std::vector<int> > get_neighbor_cells(const std::vector<int>&);
std::vector<std::vector<double> > cdist(const std::vector<std::vector<double> >&, const std::vector<std::vector<double> >&);
std::vector<std::vector<int> > simple_neighbors(const std::vector<std::vector<double> >&, double);
std::vector<std::vector<int> > cell_list_neighbors(const std::vector<std::vector<double> >&, double);