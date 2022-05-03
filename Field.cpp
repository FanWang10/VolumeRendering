#include <vector>
#include "Field.hpp"
#include "Solution.hpp"
#include "Grid.hpp"
#include <iostream>

Field::Field(Solution s, Grid g)
{
    solution = s;
    grid = g;
}

void Field::assign_grid(Grid g)
{
    grid = g;
}

void Field::assign_solution(Solution s)
{
    solution = s;
}

std::vector<double> Field::project_on_solution_coord(std::vector<double> v)
{
    std::vector<double> result = v;
    std::vector<double> min = grid.get_min_coord();
    std::vector<double> max = grid.get_max_coord();
    std::vector<double> delta = grid.get_delta();
    for(int i = 0; i < v.size(); i++)
    {
        v.at(i) = (v.at(i) - min.at(i)) / delta.at(i);
    }
    return v;
}

double Field::get_value(double x, double y, double z)
{
    /* Do Some Thing Here So That x, y, z: [-1, 1] -> [0, 255] */
    std::vector<double> pt_old {x, y, z};
    //std::cout << "pt_old " << x << " " << y << " " << z << " \n";
    std::vector<double> pt = project_on_solution_coord(pt_old);
    //std::cout << "pt " << pt.at(0) << " " << pt.at(1) << " " << pt.at(2) << "\n";
    std::vector<int> cell_index = grid.get_cell(pt);
    //std::cout << "cell_index " << cell_index.at(0) << " " << cell_index.at(1) << " " << cell_index.at(2) << " \n";
    std::vector<double> delta = grid.get_delta();
    return solution.read_at(pt.at(0), pt.at(1), pt.at(2), cell_index[0], cell_index[1], cell_index[2], delta);
}

std::vector<double> Field::get_gradient(double x, double y, double z)
{
    std::vector<double> pt {x, y, z};
    std::vector<double> pt_new = project_on_solution_coord(pt);
    std::vector<double> delta = grid.get_delta();
    std::vector<int> cell_index = grid.get_cell(pt_new);
    return solution.read_gradient_at(pt_new.at(0), pt_new.at(1), pt_new.at(2), cell_index[0], cell_index[1], cell_index[2], delta);
}

bool Field::is_point_in_bound(double x, double y, double z)
{
    return false;
}

