#include <iostream>
#include "Grid.hpp"
#include "Solution.hpp"
#include "Field.hpp"
#include "png.h"
#include <string>
#include <vector>
#include <math.h>

using namespace std;

double calcultae_function(double x, double y, double z, double a, double b)
{

    return (x * x + ((1 + b) * y) * ((1 + b) * y) + z * z -1) * (x * x + ((1 + b) * y) * ((1 + b) * y) + z * z -1) * (x * x + ((1 + b) * y) * ((1 + b) * y) + z * z - 1) - 
                            (x * x) * (z * z * z) - a * (y * y) * (z * z * z);

}

void set_color(vector<double> & color, double red, double green, double blue)
{
    color[0] = red;
    color[1] = green;
    color[2] = blue;
}

vector<double> color_lookup(double scalar_value)
{
    vector<double> color {0, 0, 0};
    double k, r, g, b;

    // Color LookUp For Heart Data
    // if(scalar_value <= 0)
    // {

    //     set_color(color, 0.9, 0.1, 0.1);
    // }

    // Color Lookup For Earth Mantle
    if(scalar_value >= 5000)
    {
        r = 0.1;
        g = 0;
        b = 0;
        set_color(color, 0.1, 0, 0);
    }else if(scalar_value >= 3000)
    {
        k = (scalar_value - 3000.0) / 2000.0;
        r = (1 - k) * 0.9 + 0.1 * k;
        g = (1 - k) * 0.1 + 0 * k;
        b = (1 - k) * 0.2 + 0 * k;
        set_color(color, r, g, b);
    }else if(scalar_value >= 2500)
    {
        k = (scalar_value - 2500.0) / 500.0;
        r = (1 - k) * 0.2 + 0.9 * k;
        g = (1 - k) * 0.1 + 0.1 * k;
        b = (1 - k) * 0.9 + 0.2 * k;
        set_color(color, r, g, b);
    }else if(scalar_value >= 2000)
    {
        k = (scalar_value - 2000.0) / 500.0;
        r = (1 - k) * 0.1 + 0.2 * k;
        g = (1 - k) * 0.1 + 0.1 * k;
        b = (1 - k) * 0.1 + 0.9 * k;
        set_color(color, r, g, b);
    }else{
        k = (scalar_value - 0) / 2000.0;
        r = (1 - k) * 0 + 0.1 * k;
        g = (1 - k) * 0 + 0.1 * k;
        b = (1 - k) * 0 + 0.1 * k;
        set_color(color, r, g, b);
    }


    return color;
} 

double alpha_lookup(double scalar_value)
{
    double alpha = 0;
    double k = 0;

    // Alpha lookup for Heart Data
    // if(scalar_value < 0)
    // {
    //     alpha = 0.8;
    // }

    // Alpha Lookup For Earth    
    if(scalar_value >= 5000)
    {
        alpha = 0;
    }else if(scalar_value >= 3000)
    {
        k = (scalar_value - 3000.0) / 2000.0;
        alpha = (1 - k) * 0.15 + k * 0;
    }else if(scalar_value >= 2500)
    {
        k = (scalar_value - 2500.0) / 500.0;
        alpha = (1 - k) * 0.01 + k * 0.15;
    }else if(scalar_value >= 2000)
    {
        k = (scalar_value - 2000.0) / 500.0;
        alpha = (1 - k) * 0.01 + k * 0.01;
    }else{
        k = (scalar_value - 0) / 2000.0;
        alpha = (1 - k) * 0 + k * 0.01;
    }

    return alpha;
}

void add_two_vector(vector<double> & v1, vector<double> & v2)
{
    v1.at(0) = v1.at(0) + v2.at(0);
    v1.at(1) = v1.at(1) + v2.at(1);
    v1.at(2) = v1.at(2) + v2.at(2);
}

vector<double> blend_cout(vector<double> & cin, vector<double> & color, double alpha, double alpha_in)
{
    vector<double> cout {0, 0, 0};
    cout[0] = cin[0] + color[0] * alpha * (1 - alpha_in);
    cout[1] = cin[1] + color[1] * alpha * (1 - alpha_in);
    cout[2] = cin[2] + color[2] * alpha * (1 - alpha_in);
    return cout;
}

double blend_aout(double alpha_in, double alpha)
{
    return alpha_in + alpha * (1 - alpha_in);
}

void normalize(vector<double> & vector)
{
    double length = vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2];
    double length_sqrt = sqrt(length);
    vector[0] = vector[0] / length_sqrt;
    vector[1] = vector[1] / length_sqrt;
    vector[2] = vector[2] / length_sqrt;
}


vector<double> calculate_L_new(vector<double> & light_position, vector<double> & location)
{
    vector<double> L  = light_position;
    L[0] = L[0] - location[0];
    L[1] = L[1] - location[1];
    L[2] = L[2] - location[2];
    normalize(L);
    return L;
}

double calculate_N_L(vector<double> & L, vector<double> & N)
{
    double dot_product = L[0] * N[0] + L[1] * N[1] + L[2] * N[2];
    double result = dot_product;
    return result;
}

void phong_shading(vector<double> & color, vector<double> & L, vector<double> & N)
{
    double N_L = calculate_N_L(L, N);
    if(N_L > 0)
    {
        color[0] = color[0] * N_L;
        color[1] = color[1] * N_L;
        color[2] = color[2] * N_L;
    }
}



int main(int argc, char *argv[])
{
    /*Task 1 Grid Memeber Vectors*/
    // Task 1 Grid Min Vector
    vector<double> min_grid_task_1 {-2, -2, -2};
    // Task 1 Grid Max Vector
    vector<double> max_grid_task_1 {2, 2, 2};
    // Task 1 Num of pts Vector
    vector<int> num_of_pts_grid_task_1 {512, 512, 512};
    // Task 1 Delta Vector
    vector<double> delta_grid_task_1 {4.0 / 511.0, 4.0 / 511.0, 4.0 / 511.0}; 

    /*Task 1 Grid*/
    Grid grid_task_1(min_grid_task_1, max_grid_task_1, num_of_pts_grid_task_1, delta_grid_task_1);

    /*Task 1 Solution Member*/
    vector<double> solution_member_task_1(512 * 512 * 512, 0);

    /*Task 1 Solution*/
    Solution solution_task_1(solution_member_task_1, 512, 512, 512);

    /*Pop Solution With Data*/
    double a, b;
    a = 1;
    b = 1;

    double x, y, z;
    x = min_grid_task_1[0];
    y = min_grid_task_1[1];
    z = min_grid_task_1[2];

    for(int i = 0; i <= num_of_pts_grid_task_1[2] - 1; i++)
    {
        for(int j = 0; j <= num_of_pts_grid_task_1[1] - 1; j++)
        {
            for(int k = 0; k <= num_of_pts_grid_task_1[0] - 1; k++)
            {
                solution_task_1.assign(k, j, i , calcultae_function(x, y, z, a, b));
                x += delta_grid_task_1[0];
            }
            x = min_grid_task_1[0];
            y += delta_grid_task_1[1];
        }
        y = min_grid_task_1[1];
        z += delta_grid_task_1[2];
    }
    /*Pop Solution With Data*/

    /*Create Task 1 Field*/
    Field field_task_1(solution_task_1, grid_task_1);

    /*Get Data From The Field*/
    std::vector<double> task_1_a_data_vector = solution_task_1.get_solution_vector();


    /*
    
        Task Data 2
    
    */
    /*Task 2 Grid Memeber Vectors*/
    // Task 2 Grid Min Vector
    vector<double> min_grid_task_2 {0, 0, 0};
    // Task 2 Grid Max Vector
    vector<double> max_grid_task_2 {255, 255, 255};
    // Task 2 Num of pts Vector
    vector<int> num_of_pts_grid_task_2 {256, 256, 256};
    // Task 2 Delta Vector
    vector<double> delta_grid_task_2 {1, 1, 1}; 

    /*Task 2 Grid*/
    Grid grid_task_2(min_grid_task_2, max_grid_task_2, num_of_pts_grid_task_2, delta_grid_task_2);

    /*Task 2 Solution Member*/
    vector<double> solution_member_task_2(256 * 256 * 256, 0);

    /*Task 2 Solution*/
    Solution solution_task_2(solution_member_task_2, 256, 256, 256);

    /*Read From raw And Pop Data Into Solution*/
    string task2_raw_name = "task2.raw";
    solution_task_2.read_from_raw(task2_raw_name);

    /*Create Task 2 Field*/
    Field field_task_2(solution_task_2, grid_task_2);



    /*
    
        Task Data 2

    */

    
    /* 
    
    
        Lab 2 
    
    
    */
    /*
        Set Resolution
    */
    int num_of_pixel_x = 512;
    int num_of_pixel_y = 512;
    int max_depth = 512;

    /* 
        Step Direction For Earth Mantle 
    */
    // YZ Plane
    //vector<double> step_ray_direction {255.0 / max_depth, 0, 0};
    // XZ Plane
    vector<double> step_ray_direction {0, 255.0 / max_depth, 0};
    // XY Plane
    //vector<double> step_ray_direction {0, 0, 255.0 / max_depth};
    /* 
        Step Direction For Earth Mantle 
    */

    /* 
        Step Direction For Heart Data 
    */
    // YZ Plane
    //vector<double> step_ray_direction {4.0 / max_depth, 0, 0};
    // XZ Plane
    //vector<double> step_ray_direction {0, 4.0 / max_depth, 0};
    // XY Plane
    //vector<double> step_ray_direction {0, 0, 4.0 / max_depth};
    /* 
        Step Direction For Earth Mantle 
    */


    // Phong Flag
    // true for using Phong, false for not
    bool phong_flag = true;


    /*

            LIBPNG

    */
    /* LIBPNG */
	FILE *fp = NULL;
	png_structp png_ptr = NULL;
	png_infop info_ptr = NULL;
	png_bytep row = NULL;

    // Change File Name
    fp = fopen("task2_b_XZ_512.png", "wb");

    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    info_ptr = png_create_info_struct(png_ptr);
    png_init_io(png_ptr, fp);
    png_set_IHDR(png_ptr, info_ptr, num_of_pixel_x, num_of_pixel_y,
			8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
			PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    png_write_info(png_ptr, info_ptr);
    row = (png_bytep) malloc(3 * num_of_pixel_x * sizeof(png_byte));
    /*
    
            LIBPNG
    
    */



    for(int j = 0; j < num_of_pixel_y; j++)
    {
        for(int i = 0; i < num_of_pixel_x; i++)
        {
            double length_of_each_pixel = (max_grid_task_1[0] - min_grid_task_1[0]) / num_of_pixel_x;
            double start_position0 = min_grid_task_1[0];
            double start_position1 = min_grid_task_1[1];
            double start_position2 = min_grid_task_1[2];

            // Initial Location In The Volume
            /*
                Initial Localtion For Earth Mantle
            */
            // XY Plane
            //vector<double> location_in_volume {(i + 0.5) * (255.0 / num_of_pixel_x), (j + 0.5) * (255.0 / num_of_pixel_y), 0};
            // XZ Plane
            vector<double> location_in_volume {(i + 0.5) * (255.0 / num_of_pixel_x), 0, (j + 0.5) * (255.0 / num_of_pixel_y)};
            // YZ Plane
            //vector<double> location_in_volume {0, (i + 0.5) * (255.0 / num_of_pixel_x), (j + 0.5) * (255.0 / num_of_pixel_y)};

            /*
                Initial Location For Heart Data
            */
            // XY Plane
            //vector<double> location_in_volume {(i + 0.5) * length_of_each_pixel + start_position0,  -2, (j + 0.5) * length_of_each_pixel + start_position2};
            // YZ Plane
            //vector<double> location_in_volume {-2,  (i + 0.5) * length_of_each_pixel + start_position0, (j + 0.5) * length_of_each_pixel + start_position2};
            // XY Plane
            //vector<double> location_in_volume {(i + 0.5) * length_of_each_pixel + start_position0, (j + 0.5) * length_of_each_pixel + start_position2, -2};

            vector<double> color {0, 0, 0};
            double alpha = 0;

            // Initial Color Black
            vector<double> cout{0, 0, 0};
            vector<double> cin = cout;

            // Initial Alpha 0
            double alpha_out = 0;
            double alpha_in = 0;

            for(int k = 0; k < max_depth; k++)
            {
                // p = p + step_ray_direction
                add_two_vector(location_in_volume, step_ray_direction);
                std::cout << "The location_in_volume is " << location_in_volume[0] << " , " << location_in_volume[1] << " , " << location_in_volume[2] << " \n";

                // Get Value For Heart
                double scalar_value = field_task_2.get_value(location_in_volume[0], location_in_volume[1], location_in_volume[2]);

                // Get Value For Earth Mantle
                //double scalar_value = field_task_1.get_value(location_in_volume[0], location_in_volume[1], location_in_volume[2]);
                

                // Color Look Ups
                color = color_lookup(scalar_value);

                // Alpha Look Ups
                alpha = alpha_lookup(scalar_value);

                if(phong_flag)
                {
                    /*
                        Use Phong
                    */
                    // Get Gradient For Task 1
                    //vector<double> gradient_N_from_data = field_task_1.get_gradient(location_in_volume[0], location_in_volume[1], location_in_volume[2]);


                    // Get Gradient For Task 2
                    vector<double> gradient_N_from_data = field_task_2.get_gradient(location_in_volume[0], location_in_volume[1], location_in_volume[2]);
                 
                    //normalize(gradient_N);
                    normalize(gradient_N_from_data);



                    /*
                        Light Position For Earth Mantle
                    */
                    // YZ Plane
                    //vector<double> light_position {-5, 255.0 / 2.0, 255.0 / 2.0};
                    // XY Plane
                    //vector<double> light_position {255.0 / 2.0, 255.0 / 2.0, -5};
                    // XZ Plane
                    vector<double> light_position {255.0 / 2.0, -5, 255.0 / 2.0};

                    /*
                        Light Position For Heart Data
                    */
                    // YZ Plane
                    //vector<double> light_position {-4, 0, 0};
                    // XY Plane
                    //vector<double> light_position {0, 0, -4};
                    // XZ Plane
                    //vector<double> light_position {0, -4, 0};


                    vector<double> L  = calculate_L_new(light_position, location_in_volume);
                    phong_shading(color, gradient_N_from_data, L);
                }

                cout = blend_cout (cin, color, alpha, alpha_in);
                alpha_out = blend_aout(alpha_in, alpha);
                cin = cout;
                alpha_in = alpha_out;
            }
            row[0 + i * 3] = 255.0 * cout[0];
            row[1 + i * 3] = 255.0 * cout[1];
            row[2 + i * 3] = 255.0 * cout[2];

        }
        png_write_row(png_ptr, row);
    }

    png_write_end(png_ptr, NULL);
    if (fp != NULL) fclose(fp);
	if (info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
	if (png_ptr != NULL) png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
	if (row != NULL) free(row);
    return 0;
}

