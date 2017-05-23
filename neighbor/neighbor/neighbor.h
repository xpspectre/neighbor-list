#pragma once

#include <vector>

std::vector<double> arange(double, double, double);
int assign_bin(double, std::vector<double>);
std::vector<std::vector<int> > get_neighbor_cells(std::vector<int>, std::vector<std::vector<int> >);
std::vector<std::vector<double> > cdist(std::vector<std::vector<double> >, std::vector<std::vector<double> >);
std::vector<std::vector<int> > simple_neighbors(std::vector<std::vector<double> >, double);
std::vector<std::vector<int> > cell_list_neighbors(std::vector<std::vector<double> >, double);