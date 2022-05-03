#ifndef FIELD_HPP
#define FIELD_HPP

#include "Solution.hpp"
#include "Grid.hpp"
#include <vector>

class Field
{
    private:
        Solution solution;
        Grid grid;
        std::vector<double> project_on_solution_coord(std::vector<double> v);
    public:
        Field(Solution s, Grid g);
        void assign_grid(Grid g);
        void assign_solution(Solution s);
        double get_value(double x, double y, double z);
        std::vector<double> get_gradient(double x, double y, double z);
        bool is_point_in_bound(double x, double y, double z);
        std::vector<int> get_length_per_dim_grid();
        std::vector<int> get_length_per_dim_solution();
};


#endif