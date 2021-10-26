#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include "bbs_eigen.h"

#define VNAME(x) #x
#define show(x) std::cout << #x << " = " << x << std::endl

void file_to_variable(std::vector<double>& array, string file_name);
void show_array(std::vector<double>& array);
void show_array_2D(std::vector<std::vector<double>>& array1);
void make2D(const std::vector<double>& array, std::vector<std::vector<double>>& array_2d, int rows, int cols);
std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& m);
void read_file_1D(std::vector<double>& array, string file_name);
void read_file_2D(std::vector<std::vector<double>>& array, string file_name, const int row, const int col);
void read_file_3D(std::vector<std::vector<std::vector<double>>>& array, string file_name, const int third_dim, const int row, const int col);
void write_to_file(ofstream& myfile_out, std::vector<std::vector<double>> particles_pos, const int nx, const int ny, const double t, string zone);
void write_to_file_irregular(ofstream& myfile_out, std::vector<std::vector<double>> particles_pos, DoubleVector2D tri, const double t, string zone);
void write_to_file_points(ofstream& myfile_out, std::vector<std::vector<double>> particles_pos, const double t, string zone);
void write_to_file_data(ofstream& myfile_out, std::vector<std::vector<double>> data, int dimension);
void write_to_file_data(ofstream& myfile_out, std::vector<double> data);
void store_data(ofstream& myfile_out, std::vector<std::vector<double>> vector);
void save_template(DoubleVector2D mesh_points, DoubleVector2D tri, Mat Template_image);
std::vector<std::vector<double>> sort_2D(const std::vector<std::vector<double>>& unsorted);
// void find_boundary(DoubleVector2D pixels_total, DoubleVector2D& boundary, DoubleVector1D& indices);
DoubleVector2D pixel_to_ray(DoubleVector2D pixels, DoubleVector2D K);
DoubleVector2D matrix_multiplication(DoubleVector2D first, DoubleVector2D second);
DoubleVector1D matrix_multiplication(DoubleVector2D first, DoubleVector1D second);
DoubleVector1D unique_indices(DoubleVector1D v);
DoubleVector2D normalize_boundary(vector<cv::Point> contour, int num_points);
DoubleVector2D normalize_boundary2(DoubleVector2D contour, int num_points);
Mat plot_points(Mat image_init, DoubleVector2D points, string window, bool line_show, bool text_show, Scalar point_clr, Scalar line_clr, Scalar string_clr, int size_image);
Mat plot_Grid(Mat image_init, DoubleVector2D points, int nx, int ny, string window, Scalar point_clr, Scalar line_clr, int size_image);
Mat plot_Triangulation(Mat image_init, DoubleVector2D shape, DoubleVector2D tri, string window, Scalar point_clr, Scalar line_clr, int size_image);
void draw_camera_head(Mat& blank_img, DoubleVector2D camera_points_Camera);
void draw_camera_body(Mat& blank_img, DoubleVector2D camera_points_Camera);
void plot_3D(DoubleVector2D Shape_3d, DoubleVector2D Tri, string Name_WINDOW, string focus_type, DoubleVector1D Camera_location, DoubleVector1D Focus_point);
DoubleVector2D projection(DoubleVector2D K, DoubleVector2D pixels);
int find_previous(DoubleVector1D vector, double value);
int find_in_vector(DoubleVector1D vector, int value);
int find_in_vector(vector<int> vector, int value);
int find_in_vector(DoubleVector1D vector, double value);
int find_in_vector(Eigen::VectorXi vector, int value);
int find_index_max(DoubleVector1D vector);
// void random_gen(std::vector<double>& visible_points, int num_points, int max, DoubleVector1D tindex);
// double link_error(Template init_template, DoubleVector2D particles_pose);
// double link_error2(Template init_template, DoubleVector2D particles_pose);
DoubleVector2D rot_matrix(string axis, double angle, string mode); 
void translate_linear(DoubleVector2D& vector, double t_x, double t_y, double t_z);
void transformation_linear(DoubleVector2D& vector, double rotx, double roty, double rotz, double t_x, double t_y, double t_z);
double random_num(int min, int max);
double random_num(double min, double max);
double initial_pose(DoubleVector2D& initial_shape_vertices_image, DoubleVector2D  K, DoubleVector1D visible_points, DoubleVector2D pixels_visible_points, int nx, int ny);
double calc_median(DoubleVector1D vector);
DoubleVector2D triangulation_simple(DoubleVector2D matches);
DoubleVector1D difference_of_vectors(DoubleVector2D v1, DoubleVector2D v2);
DoubleVector2D irregular_template_nodes(Mat image1,DoubleVector1D& bounding_box, vector<Point>& contour, int silhouette_template_num_temp = 30, int num_inside_gridPoints = 20);
DoubleVector2D removing_irregular_outsideBoundary_triangles(DoubleVector2D template_mesh_tri, vector<Point> contour, DoubleVector2D template_grid_points_2D);
VectorXd vector_to_eigen_1D(DoubleVector1D Li);
MatrixXd vector_to_eigen(DoubleVector2D Li);
DoubleVector2D eigen_to_vector(MatrixXd Ms_eigen);
double norm_vector(DoubleVector1D vector);
DoubleVector1D crossProduct(DoubleVector1D vect_A, DoubleVector1D vect_B);
DoubleVector2D warp(bbs_t bbs, int nC, double er, DoubleVector2D K, DoubleVector2D matches_template, DoubleVector2D matches_image, DoubleVector2D mesh_2D_initial);
std::vector<std::vector<double>> unique_2d(std::vector<std::vector<double>> v);
std::vector<std::vector<std::vector<double>>> groupEdges(std::vector<std::vector<double>> edgesL);
void gpu_to_cpu_features_descriptors(popsift::Features* feature_list, Mat &descriptors, DoubleVector2D &features);
void matching_dataFromGPU(Mat descriptors1, Mat descriptors2, DoubleVector2D &matches_template_project, DoubleVector2D &matches_image_project, DoubleVector2D keypoints1, DoubleVector2D keypoints2, DoubleVector1D &goodmatches_template_ids);

#endif // FUNCTIONS_H_INCLUDED
