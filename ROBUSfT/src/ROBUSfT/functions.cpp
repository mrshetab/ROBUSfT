#include "headers.h"
#include "functions.h"
#include "delaunator_h.hpp"

void file_to_variable(std::vector<double>& array, string file_name)
{
    //Opening file------------------------------------------
    ifstream inFile;
    double value;
    string address{file_name};
    inFile.open(address);
    //Counting lines------------------------------------------
//    int line_counts{0};
//    std::string line;
//    while (getline(inFile, line))
//    {
//        line_counts++;
//    }
//    cout << line_counts << endl;

//    inFile.clear();
//    inFile.seekg(0, ios::beg);
    //Reading file------------------------------------------
    if (!inFile)
    {
        cout << "\nError opening file.\n";
        //        return 13.0;
    }
    vector<double> file_reading{};
    int i{ 0 };
    while (inFile >> value)
    {
        file_reading.push_back(value);
    }
    inFile.close();
    //Output preparation--------------------------------------
    auto begin{ &file_reading[0] };
    auto end{ begin + file_reading.size() };
    array.resize(file_reading.size());
    for (auto ptr{ begin }; ptr != end; ++ptr)
    {
        array[i] = *ptr;
        i++;
    }
}
//================================================================================
//================================================================================
void show_array(std::vector<double>& array)
{
    auto begin{ &array[0] };
    auto end{ begin + array.size() };
    cout << "\n1D Array (" << array.size() << ") ============================" << endl;
    int i{ 0 };
    for (auto ptr{ begin }; ptr != end; ++ptr)
    {
        cout << *ptr << ' ';
    }
    cout << "\n";
}
//================================================================================
//================================================================================
void show_array_2D(std::vector<std::vector<double>>& array1)
{

    int rows{ static_cast<int>(array1.size()) };
    int cols{ static_cast<int>(array1[0].size()) };
    cout << "\n2D Array (" << rows << "," << cols << ") ============================" << endl;
    if (cols > 0)
    {
        for (int j = 0; j < rows; j++)
        {
            for (int k = 0; k < cols; k++)
            {
                cout << array1[j][k] << '\t';
            }
            cout << '\n';
        }
    }
    cout << "\n";
}
//================================================================================
//================================================================================
void make2D(const std::vector<double>& array, std::vector<std::vector<double>>& array_2d, int rows, int cols)
{
    array_2d.resize(rows);

    for (int i = 0; i < rows; i++)
    {
        array_2d[i].resize(cols);
    }
    int i{ 0 };
    for (int j = 0; j < rows; j++)
    {
        for (int k = 0; k < cols; k++)
        {
            array_2d[j][k] = array[i++];
        }
    }
}
//================================================================================
//================================================================================
std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& m)
{
    vector<vector<double>> result(m[0].size(), vector<double>(m.size()));

    for (vector<double>::size_type i(0); i < m[0].size(); ++i)
        for (vector<double>::size_type j(0); j < m.size(); ++j)
            result[i][j] = m[j][i];

    return result;
}
//================================================================================
//================================================================================
void read_file_1D(std::vector<double>& array, string file_name)
{
    file_to_variable(array, file_name);
}
//================================================================================
//================================================================================
void read_file_2D(std::vector<std::vector<double>>& array, string file_name, const int row, const int col)
{
    if (row < 1 || col < 1)
        cout << "Please check the row and the columns" << endl;
    else
    {
        std::vector<double> array_1d;
        file_to_variable(array_1d, file_name);
        make2D(array_1d, array, row, col);
    }
}
//================================================================================
//================================================================================
void read_file_3D(std::vector<std::vector<std::vector<double>>>& array, string file_name, const int third_dim, const int row, const int col)
{
    if (row < 1 || col < 1 || third_dim < 1)
        cout << "Please check the dimensions" << endl;
    else
    {
        array.resize(third_dim);
        std::vector<double> array_1d;
        std::vector<double> array_1d_temp(row * col);
        file_to_variable(array_1d, file_name);

        for (int i{ 0 }; i < third_dim; i++)
        {
            for (int j{ 0 }; j < row * col; j++)
            {
                array_1d_temp[j] = array_1d[i * (row * col) + j];
            }
            make2D(array_1d_temp, array[i], row, col);
        }
    }
}
//================================================================================
//================================================================================
void write_to_file(ofstream& myfile_out, std::vector<std::vector<double>> particles_pos, const int nx, const int ny, const double t, string zone)
{
    if (myfile_out.is_open())
    {
        string name = "solid";
        myfile_out << "VARIABLES = \"X\", \"Y\", \"Z\" \n";
        if (zone != "") name = zone;
        myfile_out << "ZONE T=\"" << name << "\" \n";
        myfile_out << "I=" << nx << "\n";
        myfile_out << "J=" << ny << "\n";
        myfile_out << "K=1\n";
        if (t != -100) myfile_out << "SolutionTime=" << t << "\n";
        //        int n_particles{static_cast<int>(size(particles_pos[0]))};
        for (int i{ 0 }; i < (nx * ny); i++)
        {
            myfile_out << particles_pos[0][i] << " " << particles_pos[1][i] << " " << particles_pos[2][i] << " \n";
        }
    }
}
//================================================================================
//================================================================================
void write_to_file_irregular(ofstream& myfile_out, std::vector<std::vector<double>> particles_pos, DoubleVector2D tri, const double t, string zone)
{
    if (myfile_out.is_open())
    {
        int num_points(particles_pos[0].size());
        int num_triangles(tri.size());

        string name = "solid";
        myfile_out << "VARIABLES = \"X\", \"Y\", \"Z\" \n";
        if (zone != "") name = zone;
        myfile_out << "ZONE T=\"" << name << "\" \n";
        myfile_out << " N=" << num_points << " F=FEBLOCK, ET=TRIANGLE \n";

        if (t != -100) myfile_out << "SolutionTime=" << t << "\n";
        //        int n_particles{static_cast<int>(size(particles_pos[0]))};
        for (int i{ 0 }; i < num_points; i++)
        {
            myfile_out << particles_pos[0][i] << " ";
        }
        myfile_out << "\n";
        for (int i{ 0 }; i < num_points; i++)
        {
            myfile_out << particles_pos[1][i] << " ";
        }
        myfile_out << "\n";
        for (int i{ 0 }; i < num_points; i++)
        {
            myfile_out << particles_pos[2][i] << " ";
        }
        myfile_out << "\n";
        for (int i{ 0 }; i < num_triangles; i++)
        {
            myfile_out << tri[i][0]+1 << " " << tri[i][1]+1 << " " << tri[i][2]+1 << " \n";
        }
        myfile_out << "\n";
    }
}
//================================================================================
//================================================================================
void write_to_file_points(ofstream& myfile_out, std::vector<std::vector<double>> particles_pos, const double t, string zone)
{
    if (myfile_out.is_open())
    {
        string name = "solid";
        myfile_out << "VARIABLES = \"X\", \"Y\", \"Z\" \n";
        if (zone != "") name = zone;
        myfile_out << "ZONE T=\"" << name << "\" \n";
        if (t != -100) myfile_out << "SolutionTime=" << t << "\n";
        //        int n_particles{static_cast<int>(size(particles_pos[0]))};
        for (int i{ 0 }; i < particles_pos[0].size(); i++)
        {
            myfile_out << particles_pos[0][i] << " " << particles_pos[1][i] << " " << particles_pos[2][i] << " \n";
        }
    }
}
//================================================================================
//================================================================================
void write_to_file_data(ofstream& myfile_out, std::vector<std::vector<double>> data, int dimension)
{
    if (myfile_out.is_open())
    {
        for (int j{ 0 }; j < dimension; j++)
        {
            for (int i{ 0 }; i < data[0].size(); i++)
            {
            
                myfile_out << data[j][i] << " ";
            }
            myfile_out << "\n";
        }
    }
}
//================================================================================
//================================================================================
void write_to_file_data(ofstream& myfile_out, std::vector<double> data)
{
    if (myfile_out.is_open())
    {
        for (int i{ 0 }; i < data.size(); i++)
        {

            myfile_out << data[i] << " ";
        }
    }
}
//================================================================================
//================================================================================
void store_data(ofstream& myfile_out, std::vector<std::vector<double>> vector)
{
    if (myfile_out.is_open())
    {
        int rows(vector.size());
        int cols(vector[0].size());
        for (int i{ 0 }; i < rows; i++)
        {
            for (int j{ 0 }; j < cols; j++)
            {
                myfile_out << " " << vector[i][j];
            }
            myfile_out << " \n";
        }
    }
}
//================================================================================
//================================================================================
void save_template(DoubleVector2D mesh_points, DoubleVector2D tri, Mat Template_image)
{
    ofstream file_meshpoints ("./mesh_points.txt");
    ofstream file_tri ("./tri.txt");
   
    store_data(file_meshpoints, mesh_points);
    store_data(file_tri, tri);
    imwrite("./Template_image.jpg", Template_image);

    file_meshpoints.close();
    file_tri.close();
}
//================================================================================
//================================================================================
std::vector<std::vector<double>> sort_2D(const std::vector<std::vector<double>>& unsorted)
{
    vector<vector<double>> result(unsorted.size(), vector<double>(unsorted[0].size()));
    vector<vector<double>> unsorted_trans(transpose(unsorted));
    sort(unsorted_trans[0].begin(), unsorted_trans[0].end());

    for (int i = 0; i < unsorted.size(); i++)
    {
        result[i][0] = unsorted_trans[0][i];
        for (int j = 0; j < unsorted.size(); j++)
        {
            if (unsorted_trans[0][i] == unsorted[j][0])
            {
                result[i][1] = unsorted[j][1];
                break;
            }
        }
    }
    return result;
}
//================================================================================
//================================================================================
DoubleVector2D pixel_to_ray(DoubleVector2D pixels, DoubleVector2D K)
{
    double fx{ K[0][0] };
    double fy{ K[1][1] };
    double u0{ K[0][2] };
    double v0{ K[1][2] };

    int n_visible_particles{ static_cast<int>(pixels[0].size()) };

    DoubleVector2D xy(3, vector<double>(n_visible_particles));
    DoubleVector1D xy_norm(n_visible_particles);
    DoubleVector2D particles_unit(3, vector<double>(n_visible_particles));
    for (int i = 0; i < n_visible_particles; i++)
    {
        xy[0][i] = (pixels[0][i] - u0) / fx;
        xy[1][i] = (pixels[1][i] - v0) / fy;
        xy[2][i] = 1;
        xy_norm[i] = sqrt(pow(xy[0][i], 2) + pow(xy[1][i], 2) + pow(xy[2][i], 2));
        particles_unit[0][i] = xy[0][i] / xy_norm[i];
        particles_unit[1][i] = xy[1][i] / xy_norm[i];
        particles_unit[2][i] = xy[2][i] / xy_norm[i];
    }

    return particles_unit;
}
//================================================================================
//================================================================================
DoubleVector2D matrix_multiplication(DoubleVector2D first, DoubleVector2D second)
{

    int row_1{ static_cast<int>(first.size()) };
    int col_1{ static_cast<int>(first[0].size()) };
    int row_2{ static_cast<int>(second.size()) };
    int col_2{ static_cast<int>(second[0].size()) };
    // if(col_1 != row_2) 
    //     cout << "Col_1 = " << col_1 << " and Row_2 = " << row_2 <<", Dimensions do not match!" << endl;
    DoubleVector2D result(row_1, vector<double>(col_2));
    for (int i = 0; i < row_1; i++)
    {
        for (int j = 0; j < col_2; j++)
        {
            result[i][j] = 0;
            for (int k = 0; k < col_1; k++)
            {
                result[i][j] += first[i][k] * second[k][j];
            }
        }
    }
    return result;
}
//================================================================================
//================================================================================
DoubleVector1D matrix_multiplication(DoubleVector2D first, DoubleVector1D second)
{

    int row_1{ static_cast<int>(first.size()) };
    int col_1{ static_cast<int>(first[0].size()) };
    int row_2{ static_cast<int>(second.size()) };
    // if(col_1 != row_2) 
    //     cout << "Col_1 = " << col_1 << " and Row_2 = " << row_2 <<", Dimensions do not match!" << endl;
    DoubleVector1D result(row_1);
    for (int i = 0; i < row_1; i++)
    {
        result[i] = 0;
        for (int k = 0; k < col_1; k++)
        {
            result[i] += first[i][k] * second[k];
        }
    }
    return result;
}
//================================================================================
//================================================================================
DoubleVector1D unique_indices(DoubleVector1D v)
{
    DoubleVector1D result;
    result.push_back(0);
    bool valid(true);
    for (int i = 1; i < v.size(); i++)
    {
        valid = true;
        for (int j = 0; j < i; j++)
        {
            if (v[i] == v[j]) {
                valid = false;
                break;
            }
        }
        if (valid) result.push_back(i);
    }
    return result;
}
//================================================================================
//================================================================================
DoubleVector2D normalize_boundary(vector<cv::Point> contour, int num_points)
{
    int size_initial(contour.size());
    DoubleVector2D ctr(2, vector<double>(size_initial));
    for (int i = 0; i < size_initial; i++)
    {
        ctr[0][i] = contour[i].x;
        ctr[1][i] = contour[i].y;
    }
    double contour_length(arcLength(contour, true));
    double interval(contour_length / num_points);
    DoubleVector2D contour_interpolated(2, vector<double>(2 * num_points));
    contour_interpolated[0][0] = ctr[0][0];
    contour_interpolated[1][0] = ctr[1][0];
    int index(0);
    int index_ctr(1);
    double distance_with_first_node(100000);
    while (index < 2 * num_points && index_ctr < size_initial)
    {
        double distance(sqrt(pow(ctr[0][index_ctr] - contour_interpolated[0][index], 2) + pow(ctr[1][index_ctr] - contour_interpolated[1][index], 2)));
        if (distance > interval)
        {
            double cos((ctr[0][index_ctr] - contour_interpolated[0][index]) / distance);
            double sin((ctr[1][index_ctr] - contour_interpolated[1][index]) / distance);
            contour_interpolated[0][index + 1] = contour_interpolated[0][index] + cos * interval;
            contour_interpolated[1][index + 1] = contour_interpolated[1][index] + sin * interval;
            distance_with_first_node = sqrt(pow(contour_interpolated[0][index + 1] - contour_interpolated[0][0], 2) + pow(contour_interpolated[1][index + 1] - contour_interpolated[1][0], 2));
            index++;
            if (distance_with_first_node < interval && index > 1) break;
        }
        else
        {
            index_ctr++;
        }
    }

    DoubleVector2D contour_interpolated_final(2, vector<double>(index));
    for (int i = 0; i < index; i++)
    {
        contour_interpolated_final[0][i] = contour_interpolated[0][i];
        contour_interpolated_final[1][i] = contour_interpolated[1][i];
    }
    return contour_interpolated_final;
}
//================================================================================
//================================================================================
DoubleVector2D normalize_boundary2(DoubleVector2D contour_v, int num_points)
{
    vector<cv::Point> contour(contour_v[0].size() + 1);
    for (int i = 0; i < contour_v[0].size(); i++)
    {
        contour[i].x = contour_v[0][i];
        contour[i].y = contour_v[1][i];
    }
    contour[contour_v[0].size()].x = contour_v[0][0];
    contour[contour_v[0].size()].y = contour_v[1][0];

    int size_initial(contour.size());
    DoubleVector2D ctr(2, vector<double>(size_initial));
    for (int i = 0; i < size_initial; i++)
    {
        ctr[0][i] = contour[i].x;
        ctr[1][i] = contour[i].y;
    }
    double contour_length(arcLength(contour, true));
    double interval(contour_length / num_points);
    DoubleVector2D contour_interpolated(2, vector<double>(2 * num_points));
    contour_interpolated[0][0] = ctr[0][0];
    contour_interpolated[1][0] = ctr[1][0];
    int index(0);
    int index_ctr(1);
    double distance_with_first_node(100000);
    while (index < 2 * num_points && index_ctr < size_initial)
    {
        double distance(sqrt(pow(ctr[0][index_ctr] - contour_interpolated[0][index], 2) + pow(ctr[1][index_ctr] - contour_interpolated[1][index], 2)));
        if (distance > interval)
        {
            double cos((ctr[0][index_ctr] - contour_interpolated[0][index]) / distance);
            double sin((ctr[1][index_ctr] - contour_interpolated[1][index]) / distance);
            contour_interpolated[0][index + 1] = contour_interpolated[0][index] + cos * interval;
            contour_interpolated[1][index + 1] = contour_interpolated[1][index] + sin * interval;
            distance_with_first_node = sqrt(pow(contour_interpolated[0][index + 1] - contour_interpolated[0][0], 2) + pow(contour_interpolated[1][index + 1] - contour_interpolated[1][0], 2));
            index++;
            if (distance_with_first_node < interval && index > 1) break;
        }
        else
        {
            index_ctr++;
        }
    }

    DoubleVector2D contour_interpolated_final(2, vector<double>(index));
    for (int i = 0; i < index; i++)
    {
        contour_interpolated_final[0][i] = contour_interpolated[0][i];
        contour_interpolated_final[1][i] = contour_interpolated[1][i];
    }

    return contour_interpolated_final;
}
//================================================================================
//================================================================================
Mat plot_points(Mat image_init, DoubleVector2D points, string window, bool line_show, bool text_show, Scalar point_clr, Scalar line_clr, Scalar string_clr, int size_image)
{
    if (!getWindowProperty(window, WND_PROP_VISIBLE)) // 0 if window is closed, 1 when it is open
    {
        namedWindow(window, cv::WINDOW_NORMAL);
        resizeWindow(window, 600, 600);
    }
    Vec3b white(255, 255, 255);
    Mat image;

    int rows = image_init.rows;
    if (rows != 0) {
        image = image_init;
    }
    else
    {
        image = Mat::ones(size_image, size_image, CV_8UC3);
        image.setTo(white);
    }

    /*double maxx;
    double minn;
    double max_points(0);
    double min_points(10000000);
    for (int i = 0; i < size(points[0]); i++)
    {
        for (int j = 0; j < 2; j++)
        {
            if (max_points > points[j][i]) max_points = points[j][i];
            if (min_points < points[j][i]) min_points = points[j][i];
        }
    }

    double image_min, image_max;
    minMaxLoc(image, &image_min, &image_max);

    if (1.1 * abs(max_points - min_points) > abs(image_max - image_min))
    {

    }*/

    if(points.size() > 0)
    {
        if(points[0].size() > 0)
        {
            int n_points(static_cast<int>(points[0].size()));
            int size_circle(size_image/250);
            int size_line(size_image / 500);
            double size_charac(size_image / 1000.0);

            /*double max_x(0);
            double max_y(0);
            double min_x(100000000);
            double min_y(100000000);
            for (int i = 0; i < n_points; i++)
            {
                if (points[0][i] > max_x) max_x = points[0][i];
                if (points[1][i] > max_y) max_y = points[1][i];
                if (points[0][i] < min_x) min_x = points[0][i];
                if (points[1][i] < min_y) min_y = points[1][i];
            }

            cout << min_x << " " << max_x << " " << min_y << " " << max_y << endl;
            if((max_x - min_x) > 500 || (max_y - min_y) > 500 )
                resize(image, image, Size(1.2*(max_x - min_x), 1.2 * (max_y - min_y)));*/

            for (int i = 0; i < n_points; i++) 
            {
                int index{ i };
                int index_next{ i + 1 };
                if (i == n_points - 1) { 
                    index_next = 0;
                }

                circle(image, Point2d(points[0][index], points[1][index]), size_circle, point_clr, size_circle, LINE_8, 0);
                if (line_show) line(image, Point2d(points[0][index], points[1][index]), Point2d(points[0][index_next], points[1][index_next]), line_clr, size_line, LINE_8, 0);
                if (text_show)  putText(image, to_string(i), cv::Point(points[0][index] + 5, points[1][index]), FONT_HERSHEY_COMPLEX_SMALL, size_charac, string_clr, 2* size_charac, LINE_AA);
            }    
        }               
    }
    imshow(window, image);
    return image;
}
//================================================================================
//================================================================================
Mat plot_Grid(Mat image_init, DoubleVector2D points, int nx, int ny, string window, Scalar point_clr, Scalar line_clr, int size_image)
{
    if (!getWindowProperty(window, WND_PROP_VISIBLE)) // 0 if window is closed, 1 when it is open
    {
        namedWindow(window, cv::WINDOW_NORMAL);
        resizeWindow(window, 600, 600);
    }
    Vec3b white(255, 255, 255);
    Mat image;

    int rows = image_init.rows;
    if (rows != 0) {
        image = image_init;
    }
    else
    {
        image = Mat::ones(size_image, size_image, CV_8UC3);
        image.setTo(white);
    }

    int n_points(static_cast<int>(points[0].size()));
    int size_circle(size_image / 250);
    int size_line(size_image / 200);
    double size_charac(size_image / 1000.0);

    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx - 1; i++)
        {
            //circle(image, Point2d(points[0][index], points[1][index]), size_circle, point_clr, size_circle, LINE_8, 0);
            //line(image, Point2d(points[0][i * ny + j], points[1][i * ny + j]), Point2d(points[0][i * ny + j +1], points[1][i * ny + j +1]), line_clr, size_line, LINE_8, 0);
            line(image, Point2d(points[0][j * nx + i], points[1][j * nx + i]), Point2d(points[0][j * nx + i + 1], points[1][j * nx + i +1]), line_clr, size_line, LINE_8, 0);
        }
    }

    for (int j = 0; j < ny-1; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            //circle(image, Point2d(points[0][index], points[1][index]), size_circle, point_clr, size_circle, LINE_8, 0);
            //line(image, Point2d(points[0][i * ny + j], points[1][i * ny + j]), Point2d(points[0][i * ny + j +1], points[1][i * ny + j +1]), line_clr, size_line, LINE_8, 0);
            line(image, Point2d(points[0][j * nx + i], points[1][j * nx + i]), Point2d(points[0][(j + 1) * nx + i], points[1][(j + 1) * nx + i]), line_clr, size_line, LINE_8, 0);
        }
    }
    imshow(window, image);
    return image;
}
//================================================================================
//================================================================================
void draw_camera_head(Mat& blank_img, DoubleVector2D camera_points_Camera)
{
	Point rook_points[1][3];
	rook_points[0][0]  = Point( camera_points_Camera[0][8], camera_points_Camera[1][8] );
	rook_points[0][1]  = Point( camera_points_Camera[0][9], camera_points_Camera[1][9] );
	rook_points[0][2]  = Point( camera_points_Camera[0][10], camera_points_Camera[1][10] );
	const Point* ppt[1] = { rook_points[0] };
	int npt[] = { 3 };
	fillPoly( blank_img, ppt,npt, 1,Scalar( 255, 255, 255 ),LINE_8);
	rook_points[0][0]  = Point( camera_points_Camera[0][8], camera_points_Camera[1][8] );
	rook_points[0][1]  = Point( camera_points_Camera[0][10], camera_points_Camera[1][10] );
	rook_points[0][2]  = Point( camera_points_Camera[0][11], camera_points_Camera[1][11] );
	fillPoly( blank_img, ppt,npt, 1,Scalar( 255, 255, 255 ),LINE_8);
	rook_points[0][0]  = Point( camera_points_Camera[0][8], camera_points_Camera[1][8] );
	rook_points[0][1]  = Point( camera_points_Camera[0][11], camera_points_Camera[1][11] );
	rook_points[0][2]  = Point( camera_points_Camera[0][12], camera_points_Camera[1][12] );
	fillPoly( blank_img, ppt,npt, 1,Scalar( 255, 255, 255 ),LINE_8);
	rook_points[0][0]  = Point( camera_points_Camera[0][8], camera_points_Camera[1][8] );
	rook_points[0][1]  = Point( camera_points_Camera[0][12], camera_points_Camera[1][12] );
	rook_points[0][2]  = Point( camera_points_Camera[0][9], camera_points_Camera[1][9] );
	fillPoly( blank_img, ppt,npt, 1,Scalar( 255, 255, 255 ),LINE_8);

	line(blank_img, Point2d(camera_points_Camera[0][9], camera_points_Camera[1][9]), Point2d(camera_points_Camera[0][10], camera_points_Camera[1][10]), Scalar(0, 0, 0), 1, LINE_8, 0);
	line(blank_img, Point2d(camera_points_Camera[0][10], camera_points_Camera[1][10]), Point2d(camera_points_Camera[0][11], camera_points_Camera[1][11]), Scalar(0, 0, 0), 1, LINE_8, 0);
	line(blank_img, Point2d(camera_points_Camera[0][11], camera_points_Camera[1][11]), Point2d(camera_points_Camera[0][12], camera_points_Camera[1][12]), Scalar(0, 0, 0), 1, LINE_8, 0);
	line(blank_img, Point2d(camera_points_Camera[0][12], camera_points_Camera[1][12]), Point2d(camera_points_Camera[0][9], camera_points_Camera[1][9]), Scalar(0, 0, 0), 1, LINE_8, 0);
	line(blank_img, Point2d(camera_points_Camera[0][9], camera_points_Camera[1][9]), Point2d(camera_points_Camera[0][8], camera_points_Camera[1][8]), Scalar(0, 0, 0), 1, LINE_8, 0);
	line(blank_img, Point2d(camera_points_Camera[0][10], camera_points_Camera[1][10]), Point2d(camera_points_Camera[0][8], camera_points_Camera[1][8]), Scalar(0, 0, 0), 1, LINE_8, 0);
	line(blank_img, Point2d(camera_points_Camera[0][11], camera_points_Camera[1][11]), Point2d(camera_points_Camera[0][8], camera_points_Camera[1][8]), Scalar(0, 0, 0), 1, LINE_8, 0);
	line(blank_img, Point2d(camera_points_Camera[0][12], camera_points_Camera[1][12]), Point2d(camera_points_Camera[0][8], camera_points_Camera[1][8]), Scalar(0, 0, 0), 1, LINE_8, 0);
}
//================================================================================
//================================================================================
void draw_camera_body(Mat& blank_img, DoubleVector2D camera_points_Camera)
{
	Point rook_points[1][4];
	rook_points[0][0]  = Point( camera_points_Camera[0][0], camera_points_Camera[1][0] );
	rook_points[0][1]  = Point( camera_points_Camera[0][1], camera_points_Camera[1][1] );
	rook_points[0][2]  = Point( camera_points_Camera[0][5], camera_points_Camera[1][5] );
	rook_points[0][3]  = Point( camera_points_Camera[0][4], camera_points_Camera[1][4] );
	const Point* ppt[1] = { rook_points[0] };
	int npt[] = { 4 };
	fillPoly( blank_img, ppt,npt, 1,Scalar( 0, 0, 0 ),LINE_8);
	rook_points[0][0]  = Point( camera_points_Camera[0][2], camera_points_Camera[1][2] );
	rook_points[0][1]  = Point( camera_points_Camera[0][3], camera_points_Camera[1][3] );
	rook_points[0][2]  = Point( camera_points_Camera[0][7], camera_points_Camera[1][7] );
	rook_points[0][3]  = Point( camera_points_Camera[0][6], camera_points_Camera[1][6] );
	fillPoly( blank_img, ppt,npt, 1,Scalar( 0, 0, 0 ),LINE_8);
	rook_points[0][0]  = Point( camera_points_Camera[0][0], camera_points_Camera[1][0] );
	rook_points[0][1]  = Point( camera_points_Camera[0][3], camera_points_Camera[1][3] );
	rook_points[0][2]  = Point( camera_points_Camera[0][7], camera_points_Camera[1][7] );
	rook_points[0][3]  = Point( camera_points_Camera[0][4], camera_points_Camera[1][4] );
	fillPoly( blank_img, ppt,npt, 1,Scalar( 0, 0, 0 ),LINE_8);
	rook_points[0][0]  = Point( camera_points_Camera[0][1], camera_points_Camera[1][1] );
	rook_points[0][1]  = Point( camera_points_Camera[0][2], camera_points_Camera[1][2] );
	rook_points[0][2]  = Point( camera_points_Camera[0][6], camera_points_Camera[1][6] );
	rook_points[0][3]  = Point( camera_points_Camera[0][5], camera_points_Camera[1][5] );
	fillPoly( blank_img, ppt,npt, 1,Scalar( 0, 0, 0 ),LINE_8);
}
//================================================================================
//================================================================================
Mat plot_Triangulation(Mat image_init, DoubleVector2D shape, DoubleVector2D tri, string window, Scalar point_clr, Scalar line_clr, int size_image)
{
    if (!getWindowProperty(window, WND_PROP_VISIBLE)) // 0 if window is closed, 1 when it is open
    {
        namedWindow(window, cv::WINDOW_NORMAL);
        resizeWindow(window, 600, 600);
    }
    Vec3b white(255, 255, 255);
    Mat image;

    int rows = image_init.rows;
    if (rows != 0) {
        image = image_init;
    }
    else
    {
        image = Mat::ones(size_image, size_image, CV_8UC3);
        image.setTo(white);
    }

    int n_points(static_cast<int>(shape[0].size()));
    int n_triangles(static_cast<int>(tri.size()));
    int size_circle(size_image / 250);
    int size_line(size_image / 300);
    double size_charac(size_image / 1000.0);
    
    for (int i = 0; i < n_triangles; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int index(j);
            int index_next(j+1);
            if (index == 2) index_next = 0;
            line(image, Point2d(shape[0][tri[i][index]], shape[1][tri[i][index]]), Point2d(shape[0][tri[i][index_next]], shape[1][tri[i][index_next]]), line_clr, size_line, LINE_8, 0);
        }
    }

    imshow(window, image);
    return image;
}
//================================================================================
//================================================================================
void plot_3D(DoubleVector2D Shape_3d, DoubleVector2D Tri, string Name_WINDOW = "3D", string focus_type = "point", DoubleVector1D Camera_location = {0.8, -0.5, -1}, DoubleVector1D Focus_point = {0, 0, 1})
{
	if(Shape_3d.size() == 0)
	{
		cout << "Shape 3D is empty!" << endl;
		return;
	}
	int nParticles = Shape_3d[0].size();

	DoubleVector1D z_new;
	if(focus_type == "point")
	{
		z_new = {Focus_point[0] - Camera_location[0], Focus_point[1] - Camera_location[1], Focus_point[2] - Camera_location[2]};
	}
	else if(focus_type == "object")
	{
		DoubleVector1D Object_location(3);
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < nParticles; j++)
			{
				Object_location[i] += Shape_3d[i][j];
			}
			Object_location[i] = Object_location[i]/nParticles;
		}
		z_new = {Object_location[0] - Camera_location[0], Object_location[1] - Camera_location[1], Object_location[2] - Camera_location[2]};
	}

	double norm_z = norm_vector(z_new);
	z_new[0] = z_new[0]/norm_z; 
	z_new[1] = z_new[1]/norm_z; 
	z_new[2] = z_new[2]/norm_z;
	DoubleVector1D y_temp = {0,1,0}; 
	// DoubleVector1D x_new = {0,0,-1}; 
	DoubleVector1D x_new = crossProduct(y_temp, z_new); 
	double norm_x = norm_vector(x_new);
	x_new[0] = x_new[0]/norm_x; 
	x_new[1] = x_new[1]/norm_x; 
	x_new[2] = x_new[2]/norm_x;
	DoubleVector1D y_new = crossProduct(z_new, x_new); 

	DoubleVector2D T_Camera(4, vector<double>(4));
	DoubleVector2D Rotation_matrix(3, vector<double>(3));		

	for (int i = 0; i < 3; i++)
	{
		Rotation_matrix[i][0] = x_new[i];
		Rotation_matrix[i][1] = y_new[i];
		Rotation_matrix[i][2] = z_new[i];
	}

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			T_Camera[i][j] = Rotation_matrix[i][j];
		}
	}

	T_Camera[0][3] = Camera_location[0];
	T_Camera[1][3] = Camera_location[1];
	T_Camera[2][3] = Camera_location[2];
	T_Camera[3][3] = 1;

	MatrixXd T_Camera_eigen =  vector_to_eigen(T_Camera);
	T_Camera =  eigen_to_vector(T_Camera_eigen.inverse());

	DoubleVector1D last_row(nParticles, 1);
	DoubleVector2D shape_3d_extended = Shape_3d;
	shape_3d_extended.push_back(last_row);
	DoubleVector2D Camera_Matrix{ 	{800, 	0, 		330},
									{0, 	800, 	230},
									{0,		0,		1} };
	DoubleVector2D shape_2d_Camera = projection(Camera_Matrix, matrix_multiplication(T_Camera, shape_3d_extended));

	DoubleVector2D camera_points = {{0.02, 0.02, -0.02, -0.02, 0.02, 0.02, -0.02, -0.02,  	0, 0.02, 0.02, -0.02, -0.02 }, 
									{0.02, -0.02, -0.02, 0.02, 0.02, -0.02, -0.02, 0.02,  	0, 0.02, -0.02, -0.02, 0.02 },
									{0, 0, 0, 0, -0.05, -0.05, -0.05, -0.05,		  			0, 0.02, 0.02, 0.02, 0.02  }};

	DoubleVector1D last_row_camera(camera_points[0].size(), 1);
	DoubleVector2D camera_points_extended = camera_points;
	camera_points_extended.push_back(last_row_camera);
	DoubleVector2D camera_points_Camera = projection(Camera_Matrix, matrix_multiplication(T_Camera, camera_points_extended));							

	Mat blank_img;				
	blank_img = plot_Triangulation(blank_img, shape_2d_Camera, Tri, Name_WINDOW, Scalar(255, 0, 0), Scalar(255, 0, 0), 500);
	
	if(Camera_location[2] > 0) 
	{
		draw_camera_body(blank_img, camera_points_Camera);
		draw_camera_head(blank_img, camera_points_Camera);
	}
	else
	{
		draw_camera_head(blank_img, camera_points_Camera);
		draw_camera_body(blank_img, camera_points_Camera);
	}
	
	imshow(Name_WINDOW, blank_img);
	waitKey(1);
}
//================================================================================
//================================================================================
DoubleVector2D projection(DoubleVector2D K, DoubleVector2D pixels)
{
    int num_points(pixels[0].size());

    DoubleVector2D pixels_2D_temp(matrix_multiplication(K, pixels));
    DoubleVector2D pixels_2D(2, vector<double>(num_points));
    for (int i = 0; i < num_points; i++)
    {
        pixels_2D[0][i] = pixels_2D_temp[0][i] / pixels_2D_temp[2][i];
        pixels_2D[1][i] = pixels_2D_temp[1][i] / pixels_2D_temp[2][i];
    }

    return pixels_2D;
}
//================================================================================
//================================================================================
int find_previous(DoubleVector1D vector, double value)
{
    double final;
    for (int i = 0; i < vector.size(); i++)
    {
        int index = i;
        int index_previous = i - 1;
        if (index == 0) index_previous = vector.size() - 1;
        if (vector[index] == value)
        {
            final = vector[index - 1];
        }
    }
    return final;
}
//================================================================================
//================================================================================
int find_in_vector(DoubleVector1D vector, int value)
{
    int result = -2;
    for (int i = 0; i < vector.size(); i++)
    {
        if (vector[i] == value)
        {
            result = i;
            break;
        }
    }
    return result;
}
//================================================================================
//================================================================================
int find_in_vector(DoubleVector1D vector, double value)
{
    int result = -2;
    for (int i = 0; i < vector.size(); i++)
    {
        if ( abs(vector[i] - value) < 1e-5)
        {
            result = i;
            break;
        }
    }
    return result;
}
//================================================================================
//================================================================================
int find_in_vector(vector<int> vector, int value)
{
    int result = -2;
    for (int i = 0; i < vector.size(); i++)
    {
        if (vector[i] == value)
        {
            result = i;
            break;
        }
    }
    return result;
}
//================================================================================
//================================================================================
int find_in_vector(Eigen::VectorXi vector, int value)
{
    int result = -2;
    for (int i = 0; i < vector.size(); i++)
    {
        if (vector(i) == value)
        {
            result = i;
            break;
        }
    }
    return result;
}
//================================================================================
//================================================================================
int find_index_max(DoubleVector1D vector)
{
    double max(0);
    int index;
    for (int i = 0; i < vector.size(); i++)
    {
        if (vector[i] > max)
        {
            max = vector[i];
            index = i;
        }
    }
    return index;
}
//================================================================================
//================================================================================
DoubleVector2D rot_matrix(string rot_axis, double angle, string mode)
{
    if (mode == "deg") angle = angle * M_PI / 180;
    DoubleVector2D rotaton(3, vector<double>(3));

    if (rot_axis == "x")
    {
        rotaton[0] = {1,0,0};
        rotaton[1] = {0, cos(angle), -1*sin(angle) };
        rotaton[2] = {0, sin(angle),    cos(angle) };
    }
    else if (rot_axis == "y")
    {
        rotaton[0] = { cos(angle),    0, sin(angle) };
        rotaton[1] = { 0, 1, 0 };
        rotaton[2] = { -1*sin(angle), 0, cos(angle) };
    }
    else 
    {
        rotaton[0] = { cos(angle), -1*sin(angle), 0};
        rotaton[1] = { sin(angle),    cos(angle), 0};
        rotaton[2] = { 0,0,1 };
    }

    return rotaton;
}
//================================================================================
//================================================================================
void translate_linear(DoubleVector2D& vector,double t_x, double t_y, double t_z)
{
    int dimension(static_cast<int>(vector.size()));
    int n_particles(static_cast<int>(vector[0].size()));
    for (int i = 0; i < n_particles; i++)
    {
        vector[0][i] += t_x;
        vector[1][i] += t_y;
        if(dimension==3) vector[2][i] += t_z;
    } 
}
//================================================================================
//================================================================================
void transformation_linear(DoubleVector2D& vector, double rot_x,double rot_y,double rot_z, double t_x, double t_y, double t_z)
{
    DoubleVector2D rotx(rot_matrix("x", rot_x, "deg"));
    DoubleVector2D roty(rot_matrix("y", rot_y, "deg"));
    DoubleVector2D rotz(rot_matrix("z", rot_z, "deg"));
    vector = matrix_multiplication(rotx, vector);
    vector = matrix_multiplication(roty, vector);
    vector = matrix_multiplication(rotz, vector);
    translate_linear(vector, t_x, t_y, t_z);
}
//================================================================================
//================================================================================
double random_num(int min, int max)
{
    double random_number(rand()%(max - min + 1) + min);
    return random_number;
}
//================================================================================
//================================================================================
double random_num(double min, double max)
{
    double random_number((max - min) * ((double)rand() / (double)RAND_MAX) + min);
    return random_number;
}
//================================================================================
//================================================================================
double initial_pose(DoubleVector2D& initial_shape_vertices_image, DoubleVector2D  K, DoubleVector1D visible_points, DoubleVector2D pixels_visible_points, int nx, int ny)
{
    int n_visible_particles(visible_points.size());
    double fx{ K[0][0] };
    double fy{ K[1][1] };
    double u0{ K[0][2] };
    double v0{ K[1][2] };
    double Z(0.0);
    double Z_average(0.0);
    //hamsoon kardan bayad inja anjam shavad
    // with Knowing point we can put the direction between the known points as the main direction
    double max_distance(0);
    double distance(0);
    int first_index(0);
    int sec_index(0);
    for (int i = 0; i < n_visible_particles - 1; i++)
    {
        for (int j = i + 1; j < n_visible_particles; j++)
        {
            distance = sqrt(pow(initial_shape_vertices_image[0][visible_points[i]], 2) + pow(initial_shape_vertices_image[1][visible_points[j]], 2));
            if (distance > max_distance)
            {
                max_distance = distance;
                first_index = i;
                sec_index = j;
            }
        }
    }

    DoubleVector1D vector1({ initial_shape_vertices_image[0][visible_points[first_index]] - initial_shape_vertices_image[0][visible_points[sec_index]], initial_shape_vertices_image[1][visible_points[first_index]] - initial_shape_vertices_image[1][visible_points[sec_index]] });
    DoubleVector1D vector2({ pixels_visible_points[0][first_index] - pixels_visible_points[0][sec_index], pixels_visible_points[1][first_index] - pixels_visible_points[1][sec_index] });
    //double angle(acos( (vector1[0]*vector2[0]+ vector1[1]*vector2[1])/(sqrt(pow(vector1[0],2)+ pow(vector1[1], 2))*sqrt(pow(vector2[0], 2) + pow(vector2[1], 2))) ));
    double angle = atan2(vector2[1], vector2[0]) - atan2(vector1[1], vector1[0]);

    //angle = -1 * angle;
    /*show(angle);
    cout << angle * 180 / 3.1415 << endl;*/

    // finding Z distance
    double initial_shape_total(0);
    double visible_points_total(0);
    for (int i = 0; i < n_visible_particles - 1; i++)
    {
        for (int j = i + 1; j < n_visible_particles; j++)
        {
            Z = sqrt(pow(fx * (initial_shape_vertices_image[0][visible_points[i]] - initial_shape_vertices_image[0][visible_points[j]]), 2) +
                pow(fy * (initial_shape_vertices_image[1][visible_points[i]] - initial_shape_vertices_image[1][visible_points[j]]), 2)) /
                sqrt(pow(pixels_visible_points[0][i] - pixels_visible_points[0][j], 2) +
                    pow(pixels_visible_points[1][i] - pixels_visible_points[1][j], 2));
            Z_average += Z;

            initial_shape_total += sqrt(pow((initial_shape_vertices_image[0][visible_points[i]] - initial_shape_vertices_image[0][visible_points[j]]), 2) +
                pow((initial_shape_vertices_image[1][visible_points[i]] - initial_shape_vertices_image[1][visible_points[j]]), 2));

            visible_points_total += sqrt(pow(pixels_visible_points[0][i] - pixels_visible_points[0][j], 2) +
                pow(pixels_visible_points[1][i] - pixels_visible_points[1][j], 2));
        }
    }

    double initial_shape_ave(0);
    double visible_points_ave(0);

    initial_shape_ave = sqrt(pow((initial_shape_vertices_image[0][visible_points[first_index]] - initial_shape_vertices_image[0][visible_points[sec_index]]), 2) +
        pow((initial_shape_vertices_image[1][visible_points[first_index]] - initial_shape_vertices_image[1][visible_points[sec_index]]), 2)) / initial_shape_total;

    visible_points_ave = sqrt(pow(pixels_visible_points[0][first_index] - pixels_visible_points[0][sec_index], 2) +
        pow(pixels_visible_points[1][first_index] - pixels_visible_points[1][sec_index], 2)) / visible_points_total;

    double Z_coef(visible_points_ave / initial_shape_ave);
    if (Z_coef > 1) Z_coef = 1;
    Z_average = Z_coef * abs(Z_average / ((n_visible_particles) * (n_visible_particles - 1) / 2));
    //Z_average = abs(Z_average / ((n_visible_particles) * (n_visible_particles - 1) / 2));

    //Setting values of initial surface
    DoubleVector2D temp(initial_shape_vertices_image);
    for (int i = 0; i < nx * ny; i++)
    {
        initial_shape_vertices_image[0][i] = temp[0][i] * cos(angle) - temp[1][i] * sin(angle);
        initial_shape_vertices_image[1][i] = temp[0][i] * sin(angle) + temp[1][i] * cos(angle);
        initial_shape_vertices_image[2][i] = Z_average;
    }

    DoubleVector2D xy(3, vector<double>(n_visible_particles));
    DoubleVector1D xy_norm(n_visible_particles);
    DoubleVector2D particles_unit(3, vector<double>(n_visible_particles));
    for (int i = 0; i < n_visible_particles; i++)
    {
        xy[0][i] = (pixels_visible_points[0][i] - u0) / fx;
        xy[1][i] = (pixels_visible_points[1][i] - v0) / fy;
        xy[2][i] = 1;
        xy_norm[i] = sqrt(pow(xy[0][i], 2) + pow(xy[1][i], 2) + pow(xy[2][i], 2));
        particles_unit[0][i] = xy[0][i] / xy_norm[i];
        particles_unit[1][i] = xy[1][i] / xy_norm[i];
        particles_unit[2][i] = xy[2][i] / xy_norm[i];
    }

    double min_distance(10000000);
    int target_index(0);
    DoubleVector1D target_vector;
    for (int i = 0; i < n_visible_particles; i++)
    {
        DoubleVector2D temp2(initial_shape_vertices_image);
        DoubleVector2D temp3(initial_shape_vertices_image);
        temp2[0][visible_points[i]] = Z_average * particles_unit[0][i] / particles_unit[2][i];
        temp2[1][visible_points[i]] = Z_average * particles_unit[1][i] / particles_unit[2][i];

        DoubleVector1D vect{ {temp2[0][visible_points[i]] - initial_shape_vertices_image[0][visible_points[i]],
                            temp2[1][visible_points[i]] - initial_shape_vertices_image[1][visible_points[i]],
                            temp2[2][visible_points[i]] - initial_shape_vertices_image[2][visible_points[i]] } };

        double distance(0);
        for (int j = 0; j < n_visible_particles; j++)
        {
            temp3[0][visible_points[j]] += vect[0];
            temp3[1][visible_points[j]] += vect[1];
            temp3[2][visible_points[j]] += vect[2];
            double x_element(fx * temp3[0][visible_points[j]] / temp3[2][visible_points[j]] + u0);
            double y_element(fy * temp3[1][visible_points[j]] / temp3[2][visible_points[j]] + v0);
            distance += sqrt(pow(x_element - pixels_visible_points[0][j], 2) + pow(y_element - pixels_visible_points[1][j], 2));
        }
        if (distance < min_distance)
        {
            min_distance = distance;
            target_index = i;
            target_vector = vect;
        }
    }


    for (int i = 0; i < nx * ny; i++)
    {
        initial_shape_vertices_image[0][i] += target_vector[0];
        initial_shape_vertices_image[1][i] += target_vector[1];
        initial_shape_vertices_image[2][i] += target_vector[2];
    }
}
//================================================================================
//================================================================================
double calc_median(DoubleVector1D vector)
{   
    int length(vector.size());
    double median(0);
    if (length > 0)
    {
        sort(vector.begin(), vector.end());

        if (length % 2 == 1)
        {
            median = vector[length / 2];
        }
        else
        {
            median = (vector[length / 2] + vector[(length / 2) - 1]) / 2.0;
        }
    }
    return median;
}
//================================================================================
//================================================================================
DoubleVector2D triangulation_simple(DoubleVector2D matches)
{
    std::vector<double> coords(2 * matches[0].size());
    for (int i = 0; i < matches[0].size(); i++)
    {
        coords[2 * i] = matches[0][i];
        coords[2 * i + 1] = matches[1][i];
    }

    //triangulation happens here
    delaunator::Delaunator template_d(coords);
    DoubleVector2D template_tri(template_d.triangles.size() / 3, vector<double>(3));

    for (std::size_t i = 0; i < template_d.triangles.size() / 3; i++)
    {
        template_tri[i][0] = template_d.triangles[i * 3];
        template_tri[i][1] = template_d.triangles[i * 3 + 1];
        template_tri[i][2] = template_d.triangles[i * 3 + 2];
    }

    return template_tri;
}
//================================================================================
//================================================================================
DoubleVector1D difference_of_vectors(DoubleVector2D v1, DoubleVector2D v2)
{
    int n(v1.size());
    DoubleVector2D total(n);
    DoubleVector2D common(n);
    DoubleVector1D percentage(n);
    for (int i = 0; i < n; i++)
    {
        vector<double> v(v1[i].size() + v2[i].size());
        vector<double> v_prime(v1[i].size() + v2[i].size());
        vector<double>::iterator it, st;
        vector<double>::iterator it2, st2;
        it = set_union(v1[i].begin(),
            v1[i].end(),
            v2[i].begin(),
            v2[i].end(),
            v.begin());

        for (st = v.begin(); st != it; ++st)
        {
            total[i].push_back(*st);
        }

        it2 = set_intersection(v1[i].begin(),
            v1[i].end(),
            v2[i].begin(),
            v2[i].end(),
            v_prime.begin());

        for (st2 = v_prime.begin(); st2 != it2; ++st2)
        {
            common[i].push_back(*st2);
        }

        if (total[i].size() > 0)
            percentage[i] = 100 * (total[i].size() - common[i].size()) / total[i].size();
        else
            percentage[i] = 101.0;

    }

    return percentage;
}
//================================================================================
//================================================================================
DoubleVector2D irregular_template_nodes(Mat image1,DoubleVector1D& bounding_box, vector<Point>& contour, int  silhouette_template_num_temp, int num_inside_gridPoints)
{

    Mat image1_temp1;
    bounding_box.resize(4);
	
	GaussianBlur(image1, image1_temp1, Size( 9, 9 ), 1.0);
    cv::threshold(image1_temp1, image1_temp1, 250, 255, THRESH_BINARY_INV);

    //Find the contours. Use the contourOutput Mat so the original image doesn't get overwritten
    std::vector<std::vector<cv::Point> > contours;
    cv::findContours( image1_temp1, contours, RETR_LIST, CHAIN_APPROX_NONE );

    //Extracting the largest contour
	int contour_index;
	double max_area;
    cv::Mat contourImage(image1.size(), CV_8UC3, cv::Scalar(0,0,0));
    for (int idx = 0; idx < contours.size(); idx++) 
	{
		double area = contourArea(contours[idx]);
		if(area > max_area)
		{
			max_area = area;
			contour_index= idx;
		}
    }
	
    // boudnaries of the contour
	double min_x(1000000);
	double max_x(0);
	double min_y(1000000);
	double max_y(0);
	for (int i = 0; i < contours[contour_index].size(); i++) 
	{
		if(contours[contour_index][i].x < min_x) min_x = contours[contour_index][i].x;
		if(contours[contour_index][i].x > max_x) max_x = contours[contour_index][i].x;
		if(contours[contour_index][i].y < min_y) min_y = contours[contour_index][i].y;
		if(contours[contour_index][i].y > max_y) max_y = contours[contour_index][i].y;
	}

    bounding_box[0] = min_x;
    bounding_box[1] = max_x;
    bounding_box[2] = min_y;
    bounding_box[3] = max_y;

    contour = contours[contour_index];

	srand((unsigned)time(NULL));
	
	DoubleVector2D random_point(2);
    double dist_crit = 0.5*sqrt(pow(max_x-min_x,2) + pow(max_y-min_y,2)) / num_inside_gridPoints;
	while(random_point[0].size() < num_inside_gridPoints)
	{
        bool far = true;
		double rand_x(min_x +  rand() % static_cast<int>(max_x - min_x) );
		double rand_y(min_y +  rand() % static_cast<int>(max_y - min_y) );
		double inside_outside = pointPolygonTest(contours[contour_index], Point(rand_x,rand_y), false);
		if(inside_outside == 1)	
		{
			if(pointPolygonTest(contours[contour_index], Point(rand_x,rand_y), true) > dist_crit)
			{
                for (int i = 0; i < random_point[0].size(); i++) 
	            {
                    if(sqrt(pow(rand_x - random_point[0][i],2) + pow(rand_y - random_point[1][i],2)) < dist_crit )
                    {
                        far = false;
                        break;
                    }
                }
                if(far)
                {
                    random_point[0].push_back(rand_x);
				    random_point[1].push_back(rand_y);
                }
			}
		}
	}

	DoubleVector2D contour_silhouette(normalize_boundary(contours[contour_index], silhouette_template_num_temp));
	
	DoubleVector2D template_grid_points_2D(2, vector<double>(num_inside_gridPoints + contour_silhouette[0].size()));
	for (int i = 0; i < num_inside_gridPoints; i++) 
	{
		template_grid_points_2D[0][i] = random_point[0][i];
		template_grid_points_2D[1][i] = random_point[1][i];
	}
	for (int i = 0; i < contour_silhouette[0].size(); i++) 
	{
		template_grid_points_2D[0][i + num_inside_gridPoints] = contour_silhouette[0][i];
		template_grid_points_2D[1][i + num_inside_gridPoints] = contour_silhouette[1][i];
	}

    return template_grid_points_2D;

}
//================================================================================
//================================================================================
DoubleVector2D removing_irregular_outsideBoundary_triangles(DoubleVector2D template_mesh_tri, vector<Point> contour, DoubleVector2D template_grid_points_2D)
{
    DoubleVector2D temp_tri;
	for (int i = 0; i < template_mesh_tri.size(); i++) 
	{
		double middle_x = ( template_grid_points_2D[0][template_mesh_tri[i][0]] + template_grid_points_2D[0][template_mesh_tri[i][1]] + template_grid_points_2D[0][template_mesh_tri[i][2]] ) /3.;
		double middle_y = ( template_grid_points_2D[1][template_mesh_tri[i][0]] + template_grid_points_2D[1][template_mesh_tri[i][1]] + template_grid_points_2D[1][template_mesh_tri[i][2]] ) /3.;
		double inside_outside = pointPolygonTest(contour, Point(middle_x,middle_y), false);
		if(inside_outside == 1)	
		{
			temp_tri.push_back(template_mesh_tri[i]);
		}
	}

	return temp_tri;
}
//================================================================================
//================================================================================
VectorXd vector_to_eigen_1D(DoubleVector1D Li)
{

    VectorXd mat(Li.size());
    mat = VectorXd::Map(&Li[0], Li.size());      

    return mat;
}
//================================================================================
//================================================================================
MatrixXd vector_to_eigen(DoubleVector2D Li)
{

    MatrixXd mat(Li.size(), Li[0].size());
    for (int i = 0; i < Li.size(); i++)
        mat.row(i) = VectorXd::Map(&Li[i][0], Li[i].size());

    return mat;
}
//================================================================================
//================================================================================
DoubleVector2D eigen_to_vector(MatrixXd Ms_eigen)
{
    DoubleVector2D Ms;

	for (int i=0; i<Ms_eigen.cols(); ++i)
	{
        vector<double> row_vec(Ms_eigen.col(i).data(), Ms_eigen.col(i).data() + Ms_eigen.rows());
        Ms.push_back(row_vec);
		// const double* begin = Ms_eigen.col(i).data();
		// Ms.push_back( std::vector<double>(begin, begin + Ms_eigen.rows()) );
		// Ms[i] = std::vector<double>(Ms_eigen.row(i).data(), Ms_eigen.row(i).data() + Ms_eigen.cols());
        // show_array(Ms[i]);
        // cout << i << " " << Ms_eigen.cols() << " " << begin << endl;
        // getchar();
	}

    Ms = transpose(Ms);

    return Ms;
}
//================================================================================
//================================================================================
double norm_vector(DoubleVector1D vector)
{
    int size_v(vector.size());
    double norm(0);
    for (int i=0; i<size_v; i++)
	{
        norm += pow(vector[i],2);
    }
    norm = sqrt(norm);
    return norm;
}
//================================================================================
//================================================================================
DoubleVector1D crossProduct(DoubleVector1D vect_A, DoubleVector1D vect_B)
{
    if(vect_A.size() != 3 || vect_B.size() != 3) cout << "Error Size...." << endl;
    DoubleVector1D cross_P(3);
    cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
    cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
    cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];

    return cross_P;
}
//================================================================================
//================================================================================
DoubleVector2D warp(bbs_t bbs, int nC, double er, DoubleVector2D K, DoubleVector2D matches_template, DoubleVector2D matches_image, DoubleVector2D mesh_2D_initial)
{
    
    DoubleVector2D ctrlpts_final(BBS_Function(bbs, matches_template, matches_image, K, nC, er));
    DoubleVector2D mesh_2D_transformed = BBS_Evaluation(bbs, ctrlpts_final, mesh_2D_initial, K, 0, 0);
    return mesh_2D_transformed;
}
//================================================================================
//================================================================================
std::vector<std::vector<double>> unique_2d(std::vector<std::vector<double>> v)
{
	//show_array_2D(v);
	int rows_original{ static_cast<int>(v.size()) };
	std::vector<std::vector<double>> unique_vector;
	unique_vector.push_back(v[0]);
	for (int i{ 1 }; i < rows_original; i++)
	{
		bool match = true;
		int rows_new{ static_cast<int>(unique_vector.size()) };
		for (int j{ 0 }; j < rows_new; j++)
		{
			if ((v[i][0] == unique_vector[j][0] && v[i][1] == unique_vector[j][1]) || (v[i][1] == unique_vector[j][0] && v[i][0] == unique_vector[j][1]))
			{
				match = false;
				break;
			}
		}
		if (match) unique_vector.push_back(v[i]);
	}
	return unique_vector;
}
//================================================================================
//================================================================================
std::vector<std::vector<std::vector<double>>> groupEdges(const std::vector<std::vector<double>> edgesL)
{
	std::vector<std::vector<double>> edgesL_copy(edgesL);
	std::vector<std::vector<std::vector<double>>> Groups;
	//show_array_2D(edgesL_copy);

	while (static_cast<int>(edgesL_copy.size()) > 1)
	{

		std::vector<std::vector<double>> G;
		for (int i{ 0 }; i < edgesL_copy.size(); i++)
		{
			int nx(G.size());
			std::vector<double> anedge(edgesL_copy[i]);
			std::vector<std::vector<double>> anedgeRep(2 * nx, vector<double>(2));
			for (int j{ 0 }; j < anedgeRep.size(); j++)
			{
				anedgeRep[j][0] = anedge[0];
				anedgeRep[j][1] = anedge[1];
			}
			std::vector<double> anedgeRep_linearized(4 * nx);
			for (int j{ 0 }; j < anedgeRep.size(); j++)
			{
				anedgeRep_linearized[j] = anedgeRep[j][0];
				anedgeRep_linearized[2 * nx + j] = anedgeRep[j][1];
			}

			//std::vector<std::vector<double>> G_repmat(2 * nx, vector<double>(2));
			std::vector<double> G_repmat(4 * nx);
			for (int j{ 0 }; j < G.size(); j++)
			{
				G_repmat[j] = G[j][0];
				G_repmat[nx + j] = G[j][1];
				G_repmat[2 * nx + j] = G[j][0];
				G_repmat[3 * nx + j] = G[j][1];
			}
			//if (G.size() > 0) show_array_2D(G_repmat);
			//std::vector<std::vector<double>> diffv(2 * nx, vector<double>(2));
			std::vector<double> diffv(4 * nx);
			for (int j{ 0 }; j < diffv.size(); j++)
			{
				diffv[j] = anedgeRep_linearized[j] - G_repmat[j];
			}
			/*if (G.size() > 0) show_array(anedgeRep_linearized);
			cout << "----------------------------------" << endl;
			if (G.size() > 0) show_array(G_repmat);
			cout << "----------------------------------" << endl;
			if (G.size() > 0) show_array(diffv);
			cout << "==================================" << endl;
			getchar();*/
			//if (G.size() > 0) show_array_2D(diffv);
			bool check_all_nonzero{ true };
			for (int j{ 0 }; j < diffv.size(); j++)
			{
				if (diffv[j] == 0)
				{
					check_all_nonzero = false;
					break;
				}
			}
			if (check_all_nonzero)
			{
				G.insert(G.begin(), anedge);
				//G.push_back(anedge);
			}
		}

		/*if (edgesL_copy.size() > 0) show_array_2D(edgesL_copy);
		getchar();*/

		if (G.size() > 0)
		{
			Groups.push_back(G);
			std::vector<std::vector<double>> setdiff;
			for (int k{ 0 }; k < edgesL_copy.size(); k++)
			{
				bool available = true;
				for (int j{ 0 }; j < G.size(); j++)
				{
					if ((edgesL_copy[k][0] == G[j][0] && edgesL_copy[k][1] == G[j][1]) )
					{
						available = false;
						break;
					}
				}
				if (available) setdiff.push_back(edgesL_copy[k]);
			}
			edgesL_copy = setdiff;
			sort(edgesL_copy.begin(), edgesL_copy.end());
			//if (setdiff.size() > 0) show_array_2D(setdiff);
		}
	}

	if (static_cast<int>(edgesL_copy.size()) == 1)
	{
		Groups.push_back(edgesL_copy);
	}

	return Groups;
}
//================================================================================
//================================================================================
void gpu_to_cpu_features_descriptors(popsift::Features* feature_list, Mat &descriptors, DoubleVector2D &features)
{
	int num_desriptors(feature_list->getFeatureCount());
	popsift::Feature* features1 = feature_list->getFeatures();
	for (int i = 0; i < num_desriptors; i++)
	{
		features[0][i] = features1[i].xpos;
		features[1][i] = features1[i].ypos;

		for (int j = 0; j < 128; j++)
		{
			popsift::Descriptor* descc = *features1[i].desc;
			descriptors.at<float>(Point(j,i)) = descc->features[j];
		}
	}
}
//================================================================================
//================================================================================
void matching_dataFromGPU(Mat descriptors1, Mat descriptors2, DoubleVector2D &matches_template_project, DoubleVector2D &matches_image_project, DoubleVector2D keypoints1, DoubleVector2D keypoints2, DoubleVector1D &goodmatches_template_ids)
{
	// Brute Force Matching
		BFMatcher matcher(NORM_L2, false);
		//std::vector< DMatch > matches;
		std::vector< std::vector<DMatch> > matches;
		//matcher.match(descriptors1, descriptors2, matches);

		matcher.knnMatch(descriptors1, descriptors2, matches, 2);

		//-- Filter matches using the Lowe's ratio test
		//std::vector< DMatch > good_matches(matches);
		const float ratio_thresh = 0.8f;
		std::vector<DMatch> good_matches;
		for (int i = 0; i < matches.size(); i++)
		{
			if (matches[i][0].distance < ratio_thresh * matches[i][1].distance)
			{
				good_matches.push_back(matches[i][0]);
			}
		}


		//-- Localize the object
		// std::vector<Point2f> obj(good_matches.size());
		// std::vector<Point2f> scene(good_matches.size());
		// for (size_t i = 0; i < good_matches.size(); i++)
		// {
		// 	//-- Get the keypoints from the good matches
		// 	obj[i] = keypoints1[0][good_matches[i].queryIdx].pt;
		// 	scene[i] = keypoints2[good_matches[i].trainIdx].pt;
		// }

		DoubleVector2D matches_template_project_0(2, vector<double>(good_matches.size()));
		DoubleVector2D matches_image_project_0(2, vector<double>(good_matches.size()));
		goodmatches_template_ids.resize(good_matches.size());

		for (size_t i = 0; i < good_matches.size(); i++)
		{
			goodmatches_template_ids[i] = good_matches[i].queryIdx;
			matches_template_project_0[0][i] = keypoints1[0][good_matches[i].queryIdx];
			matches_template_project_0[1][i] = keypoints1[1][good_matches[i].queryIdx];
			matches_image_project_0[0][i] = keypoints2[0][good_matches[i].trainIdx];
			matches_image_project_0[1][i] = keypoints2[1][good_matches[i].trainIdx];
		}

		matches_template_project = matches_template_project_0;
		matches_image_project = matches_image_project_0;
}
