#include <iostream>
#include "headers.h"
#include "delaunator_h.hpp"
#include "functions.h"
#include "bbs_eigen.h"
#include "ROBUSfT.h"
#include "thread"
#include "mutex"
#include "condition_variable"
#include <omp.h>

mutex my_mutex;
bool job_done = false;
condition_variable con_var;

#if POPSIFT_IS_DEFINED(POPSIFT_USE_NVTX)
#include <nvToolsExtCuda.h>
#else
#define nvtxRangePushA(a)
#define nvtxRangePop()
#endif

void Robusft::set_template_image(string address)
	{
		Template_image = imread(address);
		cvtColor(Template_image, Template_image_BW, COLOR_BGR2GRAY);
	}

void Robusft::setCalibrationMat(DoubleVector2D calibrationMat)
	{
		K = calibrationMat;
	}

void Robusft::build_template(bool load_template)
	{	
		cout << "####################################################" << endl;
		cout << "Creating Template ......." << endl;
		// Estimating the dimension of the object 
		int rows = Template_image.rows;
		int cols = Template_image.cols;

		Template.Mesh_project_min_x = 0;
		Template.Mesh_project_max_x = cols;
		Template.Mesh_project_min_y = 0;
		Template.Mesh_project_max_y = rows;

		Template.Mesh_z = K[0][0]*Template.Width/(Template.Mesh_project_max_x - Template.Mesh_project_min_x);
		Template.Mesh_min_x = Template.Mesh_z*(Template.Mesh_project_min_x - K[0][2]) / K[0][0];
		Template.Mesh_max_x = Template.Mesh_z*(Template.Mesh_project_max_x - K[0][2]) / K[0][0];
		Template.Mesh_min_y = Template.Mesh_z*(Template.Mesh_project_min_y - K[1][2]) / K[1][1];
		Template.Mesh_max_y = Template.Mesh_z*(Template.Mesh_project_max_y - K[1][2]) / K[1][1];

		Template.Sx = (Template.Mesh_max_x - Template.Mesh_min_x);
		Template.Sy = (Template.Mesh_max_y - Template.Mesh_min_y);

		// Initialization
		DoubleVector1D coords;
		DoubleVector2D mesh;
		DoubleVector2D template_grid_points_2D;	
		DoubleVector1D bounding_box;  
		vector<Point> contour;	

		if(Template.Mesh_type == "regular")
		{
			cout << "Template type is regular" << endl;
			if(Template.Nx == 0 || Template.Nx == 0) cout << "Template.Nx or Template.Ny is 0. Please check their value again. Are you sure your template is regular?" << endl;
			Template.nParticles = Template.Nx * Template.Ny;
			template_grid_points_2D.resize(2, vector<double>(Template.nParticles));
			mesh.resize(Template.nParticles, vector<double>(3));
			DoubleVector2D x(Template.Nx, vector<double>(Template.Ny));
			DoubleVector2D y(Template.Nx, vector<double>(Template.Ny));
			coords.resize(2 * Template.nParticles);
			
			int k = 0;
			int s = 0;
			for (int j = 0; j < Template.Ny; j++) {
				for (int i = 0; i < Template.Nx; i++) {
					x[i][j] = Template.Mesh_min_x + i * Template.Sx / (Template.Nx - 1);
					y[i][j] = Template.Mesh_min_y + j * Template.Sy / (Template.Ny - 1);
					template_grid_points_2D[0][s] = 0 + i * (Template.Mesh_project_max_x - Template.Mesh_project_min_x) / (Template.Nx - 1);
					template_grid_points_2D[1][s] = 0 + j * (Template.Mesh_project_max_y - Template.Mesh_project_min_y) / (Template.Ny - 1);
					mesh[s][0] = x[i][j];
					mesh[s][1] = y[i][j];
					mesh[s][2] = Template.Mesh_z;
					coords[k] = x[i][j];
					coords[k + 1] = y[i][j];
					s++;
					k = k + 2;
				}
			}

			Template.nGrid_on_boundary = 2*Template.Nx + 2*(Template.Ny-2);
			Template.nGrid_on_boundary = Template.nParticles - Template.nGrid_on_boundary;
		}
		else if(Template.Mesh_type == "irregular")
		{
			cout << "Template type is irregular" << endl;
			template_grid_points_2D = irregular_template_nodes(Template_image_BW, bounding_box, contour, Template.nGrid_on_boundary, Template.nGrid_inside_boundary);
			Template.nParticles = template_grid_points_2D[0].size();
			mesh.resize(Template.nParticles, vector<double>(3));
			coords.resize(2 * Template.nParticles);
			for (int i = 0; i < Template.nParticles; i++) 
			{
				mesh[i][0] = Template.Mesh_z*(template_grid_points_2D[0][i] - K[0][2]) / K[0][0];
				mesh[i][1] = Template.Mesh_z*(template_grid_points_2D[1][i] - K[1][2]) / K[1][1];
				mesh[i][2] = Template.Mesh_z;
			}
			for (int i = 0; i < Template.nParticles; i++)
			{
					coords[2*i] = mesh[i][0];
					coords[2*i + 1] = mesh[i][1];
			}
			Template.Nx = Template.nParticles;
			Template.Ny = 1;
		}

		//triangulation happens here
		delaunator::Delaunator d(coords);
		DoubleVector2D tri(d.triangles.size() / 3, vector<double>(3));

		for (std::size_t i = 0; i < d.triangles.size() / 3; i++)
		{
			tri[i][0] = d.triangles[i * 3];
			tri[i][1] = d.triangles[i * 3 + 1];
			tri[i][2] = d.triangles[i * 3 + 2];
		}

		if(Template.Mesh_type == "irregular") tri = removing_irregular_outsideBoundary_triangles(tri, contour, template_grid_points_2D);
		cout << "Triangulation is done." << endl;

		DoubleVector2D combined(3 * tri.size(), vector<double>(2));
		for (int i = 0; i < tri.size(); i++)
		{
			combined[i][0] = tri[i][1];
			combined[i][1] = tri[i][0];

			combined[tri.size() + i][0] = tri[i][2];
			combined[tri.size() + i][1] = tri[i][1];

			combined[2 * tri.size() + i][0] = tri[i][0];
			combined[2 * tri.size() + i][1] = tri[i][2];
		}

		for (int i = 0; i < combined.size(); i++)
		{
			sort(combined[i].begin(), combined[i].end());
		}

		DoubleVector2D stretching_edge_list(unique_2d(combined));
		sort(stretching_edge_list.begin(), stretching_edge_list.end());

		DoubleVector3D SGroups(groupEdges(stretching_edge_list));
		DoubleVector2D  Sgroup_length0(SGroups.size());
		DoubleVector2D  Sgroup_k(SGroups.size());
		for (int i = 0; i < SGroups.size(); i++)
		{

			DoubleVector2D group_edges(SGroups[i]);
			int n_edges(group_edges.size());
			DoubleVector2D p1(3, vector<double>(n_edges));
			DoubleVector2D p2(3, vector<double>(n_edges));
			for (int j = 0; j < n_edges; j++)
			{
				p1[0][j] = mesh[group_edges[j][0]][0];
				p1[1][j] = mesh[group_edges[j][0]][1];
				p1[2][j] = mesh[group_edges[j][0]][2];
				p2[0][j] = mesh[group_edges[j][1]][0];
				p2[1][j] = mesh[group_edges[j][1]][1];
				p2[2][j] = mesh[group_edges[j][1]][2];
			}
			Sgroup_length0[i].resize(n_edges);
			Sgroup_k[i].resize(n_edges);

			for (int j = 0; j < n_edges; j++)
			{
				Sgroup_length0[i][j] = sqrt(pow(p2[0][j] - p1[0][j], 2) + pow(p2[1][j] - p1[1][j], 2) + pow(p2[2][j] - p1[2][j], 2));
				Sgroup_k[i][j] = Template.K_stretch; 
			}
		}
		cout << "Stretching Groups are created." << endl;

		// Bending Constraints
		int nt(tri.size());
		for (int i = 0; i < nt; i++)
		{
			sort(tri[i].begin(), tri[i].end());
		}

		std::vector<double> trinum(3 * nt);
		for (int i = 0; i < nt; i++)
		{
			trinum[i] = i;
			trinum[nt + i] = i;
			trinum[2 * nt + i] = i;
		}

		DoubleVector2D edges(3 * nt, vector<double>(2));
		for (int i = 0; i < nt; i++)
		{
			edges[i][0] = tri[i][0];
			edges[i][1] = tri[i][1];

			edges[nt + i][0] = tri[i][0];
			edges[nt + i][1] = tri[i][2];

			edges[2 * nt + i][0] = tri[i][1];
			edges[2 * nt + i][1] = tri[i][2];
		}

		DoubleVector2D edges_unsorted{ edges };
		sort(edges.begin(), edges.end());
		DoubleVector1D tags(3 * nt);
		DoubleVector1D temp_zero;
		DoubleVector1D help_index;
		for (int i = 0; i < 3 * nt; i++)
		{
			if (i > 0)
			{
				if (edges[i][0] != edges[i - 1][0] || edges[i][1] != edges[i - 1][1])
				{
					help_index = temp_zero;
				}
			}
			for (int j = 0; j < 3 * nt; j++)
			{
				if (edges_unsorted[j][0] == edges[i][0] && edges_unsorted[j][1] == edges[i][1])
				{
					int indexx = find_in_vector(help_index, j);
					if (indexx == -2)
					{
						help_index.push_back(j);
						tags[i] = j;
						break;
					}
				}
			}
		}

		DoubleVector1D trinum_new{ trinum };
		for (int i = 0; i < trinum.size(); i++)
		{
			trinum_new[i] = trinum[tags[i]];
		}
		trinum = trinum_new;

		// adjacent triangles list
		DoubleVector2D diff(3 * nt - 1, vector<double>(2));
		DoubleVector1D all_diff(3 * nt - 1);
		DoubleVector1D k_all_diff;
		for (int i = 0; i < (3 * nt - 1); i++)
		{
			diff[i][0] = edges[i + 1][0] - edges[i][0];
			diff[i][1] = edges[i + 1][1] - edges[i][1];
			if (diff[i][0] == 0 && diff[i][1] == 0)
			{
				all_diff[i] = 1;
				k_all_diff.push_back(i);
			}
		}

		DoubleVector2D tri_list(k_all_diff.size(), vector<double>(2));
		for (int i = 0; i < tri_list.size(); i++)
		{
			tri_list[i][0] = trinum[k_all_diff[i]];
			tri_list[i][1] = trinum[k_all_diff[i] + 1];
		}

		sort(tri_list.begin(), tri_list.end());

		// create bending constraints
		DoubleVector2D bending_edge_list(tri_list.size(), vector<double>(2));
		for (int i = 0; i < tri_list.size(); i++)
		{
			DoubleVector1D t1(3);
			DoubleVector1D t2(3);
			t1[0] = tri[tri_list[i][0]][0];
			t1[1] = tri[tri_list[i][0]][1];
			t1[2] = tri[tri_list[i][0]][2];
			t2[0] = tri[tri_list[i][1]][0];
			t2[1] = tri[tri_list[i][1]][1];
			t2[2] = tri[tri_list[i][1]][2];

			double ip3{ 0 };
			double ip4{ 0 };
			for (int j = 0; j < 3; j++)
			{
				ip3 = t1[j];
				bool no_match{ true };
				for (int k = 0; k < 3; k++)
				{
					if (t1[j] == t2[k])
					{
						no_match = false;
					}
				}
				if (no_match) break;
			}
			for (int j = 0; j < 3; j++)
			{
				ip4 = t2[j];
				bool no_match{ true };
				for (int k = 0; k < 3; k++)
				{
					if (t2[j] == t1[k])
					{
						no_match = false;
					}
				}
				if (no_match) break;
			}
			bending_edge_list[i][0] = ip3;
			bending_edge_list[i][1] = ip4;
		}

		DoubleVector3D BGroups(groupEdges(bending_edge_list));


		DoubleVector2D  Bgroup_length0(BGroups.size());
		DoubleVector2D  Bgroup_k(BGroups.size());

		for (int i = 0; i < BGroups.size(); i++)
		{
			DoubleVector2D group_edges(BGroups[i]);
			int n_edges(group_edges.size());

			DoubleVector2D p1(3, vector<double>(n_edges));
			DoubleVector2D p2(3, vector<double>(n_edges));

			for (int j = 0; j < n_edges; j++)
			{
				p1[0][j] = mesh[group_edges[j][0]][0];
				p1[1][j] = mesh[group_edges[j][0]][1];
				p1[2][j] = mesh[group_edges[j][0]][2];
				p2[0][j] = mesh[group_edges[j][1]][0];
				p2[1][j] = mesh[group_edges[j][1]][1];
				p2[2][j] = mesh[group_edges[j][1]][2];
			}

			Bgroup_length0[i].resize(n_edges);
			Bgroup_k[i].resize(n_edges);

			for (int j = 0; j < n_edges; j++)
			{
				Bgroup_length0[i][j] = sqrt(pow(p2[0][j] - p1[0][j], 2) + pow(p2[1][j] - p1[1][j], 2) + pow(p2[2][j] - p1[2][j], 2));
				Bgroup_k[i][j] = Template.K_bend; 
			}
		}
		cout << "Bending Groups are created." << endl;

		Template.SGroups = SGroups;
		Template.Sgroups_length0 = Sgroup_length0;
		Template.Sgroups_k = Sgroup_k;
		Template.nSGroups = SGroups.size();

		Template.BGroups = BGroups;
		Template.Bgroups_length0 = Bgroup_length0;
		Template.Bgroups_k = Bgroup_k;
		Template.nBGroups = BGroups.size();
		Template.Vertices = mesh;
		Template.Tri = tri;

		// Moving the inial mesh to a new location in 3D space (optional)
		Template.Mesh = transpose(Template.Vertices);
		// translate_linear(Template.Mesh, Template.Mesh_min_x, Template.Mesh_min_y, Template.Mesh_z - 0.01);
		Template.Mesh_project = projection(K, Template.Mesh);

		Shape_3d = Template.Mesh;
		Shape_2d = Template.Mesh_project;


		// Extracting Keypoints on the Template Texture-map ========
		extract_keypoints_GPU("Template");

		// Extracting Barycenteric Coordinates of Keypoints ========
		get_BaryCentricCoordsOfKeypoints();

		// Saving Template =========================================
		save_template(Template.Vertices, Template.Tri, Template_image);

		cout << "Number of mesh points = " << Template.nParticles << endl;
		cout << "Number of mesh points in X direction = " << Template.Nx << endl;
		cout << "Number of mesh points in Y direction = " << Template.Ny << endl;
		cout << "Number of Stretching Groups = " << Template.nSGroups << endl;
		cout << "Number of Bending Groups = " << Template.nBGroups << endl;

		if(ROBUSfT_plot)
		{
			namedWindow("Template", cv::WINDOW_NORMAL);
			resizeWindow("Template", 400, 400);
			Mat temp_img = Template_image.clone();
			plot_Triangulation(temp_img, Template.Mesh_project, tri, "Template", Scalar(0, 255, 0), Scalar(0, 255, 0), 600);
			plot_points(temp_img, Template.Mesh_project, "Template", false, false, Scalar(0, 0, 255), Scalar(255, 0, 0), Scalar(0, 255, 0), 600);
			waitKey(1);
		}
		
		cout << "Creating Template is finished" << endl;
		cout << "####################################################" << endl;
	}

void Robusft::extract_keypoints_GPU(string type, VideoCapture cap)
{
	popsift::Config config;
	PopSift PopSift( config, popsift::Config::ExtractingMode);

	if(type == "Template")
	{
		const auto w = Template_image_BW.cols;
		const auto h = Template_image_BW.rows;
		uchar *template_data = Template_image_BW.data;
		SiftJob* job = PopSift.enqueue( w, h, template_data );
		popsift::Features* feature_list = job->get();
		int nKeys = feature_list->getFeatureCount();
		Mat desc(Size(128, nKeys), CV_32F);
		DoubleVector2D features(2, vector<double>(nKeys));
		gpu_to_cpu_features_descriptors(feature_list, desc, features);
		Template.nKeypoints = nKeys;
		Template.Descriptors = desc;
		Template.Features = features;
		delete feature_list;
		delete job;
		PopSift.uninit();
	}
	else if(type == "Image")
	{
		while(true)
		{
			capture_image(cap);

			job_done = false;
			lock_guard<mutex> lg(my_mutex);
			
			// Extracting Keypoints in the Image ===========================
			auto start_Keypoint = chrono::high_resolution_clock::now();

			const auto w = Captured_image_BW.cols;
			const auto h = Captured_image_BW.rows;
			uchar *template_data = Captured_image_BW.data;
			SiftJob* job = PopSift.enqueue( w, h, template_data );
			popsift::Features* feature_list = job->get();
			int nKeys = feature_list->getFeatureCount();
			Mat desc(Size(128, nKeys), CV_32F);
			DoubleVector2D features(2, vector<double>(nKeys));
			gpu_to_cpu_features_descriptors(feature_list, desc, features);
			nKeypoints = nKeys;
			Descriptors = desc;
			Features = features;

			delete feature_list;
			delete job;

			auto stop_Keypoint = chrono::high_resolution_clock::now();
			if(ROBUSfT_report) cout << "Keypoint Extraction GPU Thread = " << std::chrono::duration<double>(stop_Keypoint - start_Keypoint).count() << " (s)" << endl;

			job_done = true;
			con_var.notify_one();
		}
	}

	PopSift.uninit();		
}

void Robusft::get_BaryCentricCoordsOfKeypoints()
	{
		Template.nKeypoints = Template.Features[0].size();
		Template.Keypoints_BC_related_triangle.resize(Template.nKeypoints); // triangle that surrounds the match in templates... related_particles[number of matches_in_template] = row number of tri
		Template.Keypoints_BC_Coefficients.resize(3, vector<double>(Template.nKeypoints)); //coefficients of barycentric coordinates
		Template.Keypoints_BC_related_meshpoints.resize(3, vector<double>(Template.nKeypoints));
		for (int i = 0; i < Template.nKeypoints; i++)
		{
			double x(Template.Features[0][i]);
			double y(Template.Features[1][i]);

			for (int j = 0; j < Template.Tri.size(); j++)
			{
				double x1(Template.Mesh_project[0][Template.Tri[j][0]]);
				double x2(Template.Mesh_project[0][Template.Tri[j][1]]);
				double x3(Template.Mesh_project[0][Template.Tri[j][2]]);
				double y1(Template.Mesh_project[1][Template.Tri[j][0]]);
				double y2(Template.Mesh_project[1][Template.Tri[j][1]]);
				double y3(Template.Mesh_project[1][Template.Tri[j][2]]);
				double landa_1(((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)));
				double landa_2(((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)));
				double landa_3(1 - landa_1 - landa_2);
				if (landa_1 > 0 && landa_2 > 0 && landa_3 > 0)
				{
					Template.Keypoints_BC_related_triangle[i] = j;
					Template.Keypoints_BC_Coefficients[0][i] = landa_1;
					Template.Keypoints_BC_Coefficients[1][i] = landa_2;
					Template.Keypoints_BC_Coefficients[2][i] = landa_3;
					Template.Keypoints_BC_related_meshpoints[0][i] = Template.Tri[j][0];
					Template.Keypoints_BC_related_meshpoints[1][i] = Template.Tri[j][1];
					Template.Keypoints_BC_related_meshpoints[2][i] = Template.Tri[j][2];
					break;
				}
			}
		}
	}

void Robusft::capture_image(VideoCapture cap)
	{
		bool bSuccess = cap.read(Captured_image); // read a new frame from video
		cvtColor(Captured_image, Captured_image_BW, COLOR_BGR2GRAY);
		if (!bSuccess) //if not success, break loop
		{
			cout << "Cannot read a frame from video stream" << endl;
		}
	}

void Robusft::match(float ratio_thresh)
	{
		// Brute Force Matching
		BFMatcher matcher(NORM_L2, false);
		std::vector< std::vector<DMatch> > matches;
		matcher.knnMatch(Template.Descriptors, Descriptors, matches, 2);

		//-- Filter matches using the Lowe's ratio test
		// const float ratio_thresh = 0.8f;
		
		std::vector<DMatch> good_matches;
		for (int i = 0; i < matches.size(); i++)
		{
			if (matches[i][0].distance < ratio_thresh * matches[i][1].distance)
			{
				good_matches.push_back(matches[i][0]);
			}
		}

		nMatches = good_matches.size();

		Matches_template = std::vector<std::vector<double>>(2, vector<double>(nMatches));
		Matches_image = std::vector<std::vector<double>>(2, vector<double>(nMatches));
		Goodmatches_ids = vector<double>(nMatches);

		for (size_t i = 0; i < nMatches; i++)
		{
			Goodmatches_ids[i] = good_matches[i].queryIdx;
			Matches_template[0][i] = Template.Features[0][good_matches[i].queryIdx];
			Matches_template[1][i] = Template.Features[1][good_matches[i].queryIdx];
			Matches_image[0][i] = Features[0][good_matches[i].trainIdx];
			Matches_image[1][i] = Features[1][good_matches[i].trainIdx];
		}
	}

void Robusft::outlier_removal_algorithm(int crit_nMatches_stepI, double crit_MAD, double alpha_s)
	{
		if(ROBUSfT_plot)
		{
			namedWindow("Outlier removal, Step I", cv::WINDOW_NORMAL);
			resizeWindow("Outlier removal, Step I", 400, 400);
			namedWindow("Outlier removal, Step II", cv::WINDOW_NORMAL);
			resizeWindow("Outlier removal, Step II", 400, 400);
			namedWindow("Outlier removal, Step III", cv::WINDOW_NORMAL);
			resizeWindow("Outlier removal, Step III", 400, 400);
		}
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// Step I
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		auto start_stepI = chrono::high_resolution_clock::now();

		// Generating triangulations
		DoubleVector2D template_tri(triangulation_simple(Matches_template));
		DoubleVector2D image_tri(triangulation_simple(Matches_image));

		//finding neighbors in the template and in the image
		DoubleVector2D template_neighbors(nMatches);
		int other_1(0);
		int other_2(0);
		for (int i = 0; i < template_tri.size(); i++)
		{
			for (int j = 0; j < 3; j++)
			{
				switch (j) {
				case 0:
					other_1 = 1;
					other_2 = 2;
					break;
				case 1:
					other_1 = 0;
					other_2 = 2;
					break;
				case 2:
					other_1 = 0;
					other_2 = 1;
					break;
				}

				if (find_in_vector(template_neighbors[template_tri[i][j]], template_tri[i][other_1]) == -2)
					template_neighbors[template_tri[i][j]].push_back(template_tri[i][other_1]);
				if (find_in_vector(template_neighbors[template_tri[i][j]], template_tri[i][other_2]) == -2)
					template_neighbors[template_tri[i][j]].push_back(template_tri[i][other_2]);
			}
		}

		DoubleVector2D image_neighbors(nMatches);
		for (int i = 0; i < image_tri.size(); i++)
		{
			for (int j = 0; j < 3; j++)
			{
				switch (j) {
				case 0:
					other_1 = 1;
					other_2 = 2;
					break;
				case 1:
					other_1 = 0;
					other_2 = 2;
					break;
				case 2:
					other_1 = 0;
					other_2 = 1;
					break;
				}

				if (find_in_vector(image_neighbors[image_tri[i][j]], image_tri[i][other_1]) == -2)
					image_neighbors[image_tri[i][j]].push_back(image_tri[i][other_1]);
				if (find_in_vector(image_neighbors[image_tri[i][j]], image_tri[i][other_2]) == -2)
					image_neighbors[image_tri[i][j]].push_back(image_tri[i][other_2]);
			}
		}

		DoubleVector1D OF(difference_of_vectors(template_neighbors, image_neighbors));

		//finding average
		double OF_th;
		double sumTotal = 0;
		for (int k = 0; k < OF.size(); ++k) 
		{
			sumTotal += OF[k];
		}
		OF_th = sumTotal / OF.size();

		//selecting inliers based on the percentage_difference
		DoubleVector2D Matches_template_stepI(2);
		DoubleVector2D Matches_image_stepI(2);
		DoubleVector1D id_selectedMatchPoint;
		int nMatches_stepI(0);
		for (int i = 0; i < nMatches; i++)
		{
			if (OF[i] < OF_th)
			{
				Matches_template_stepI[0].push_back(Matches_template[0][i]);
				Matches_template_stepI[1].push_back(Matches_template[1][i]);

				Matches_image_stepI[0].push_back(Matches_image[0][i]);
				Matches_image_stepI[1].push_back(Matches_image[1][i]);
				id_selectedMatchPoint.push_back(i);
				nMatches_stepI++;
			}
		}

		if (nMatches_stepI > crit_nMatches_stepI)
		{
			// Calculating Warp W1
			Mesh_transformed_stepI = warp(bbs, Warp_nC, Warp_er, K, Matches_template_stepI, Matches_image_stepI, Template.Mesh_project);

			auto stop_stepI = chrono::high_resolution_clock::now();
			
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			// Step II
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			auto start_stepII = chrono::high_resolution_clock::now();

			// Evaluating Distances using barycentric coordinates 
			DoubleVector2D matches_image_stepII_barrycentric(2, vector<double>(nMatches_stepI));
			for (int i = 0; i < nMatches_stepI; i++)
			{
				int match_index(Goodmatches_ids[id_selectedMatchPoint[i]]);
				matches_image_stepII_barrycentric[0][i] = Template.Keypoints_BC_Coefficients[0][match_index] * Mesh_transformed_stepI[0][Template.Tri[Template.Keypoints_BC_related_triangle[match_index]][0]] +
					Template.Keypoints_BC_Coefficients[1][match_index] * Mesh_transformed_stepI[0][Template.Tri[Template.Keypoints_BC_related_triangle[match_index]][1]] +
					Template.Keypoints_BC_Coefficients[2][match_index] * Mesh_transformed_stepI[0][Template.Tri[Template.Keypoints_BC_related_triangle[match_index]][2]];

				matches_image_stepII_barrycentric[1][i] = Template.Keypoints_BC_Coefficients[0][match_index] * Mesh_transformed_stepI[1][Template.Tri[Template.Keypoints_BC_related_triangle[match_index]][0]] +
					Template.Keypoints_BC_Coefficients[1][match_index] * Mesh_transformed_stepI[1][Template.Tri[Template.Keypoints_BC_related_triangle[match_index]][1]] +
					Template.Keypoints_BC_Coefficients[2][match_index] * Mesh_transformed_stepI[1][Template.Tri[Template.Keypoints_BC_related_triangle[match_index]][2]];
			}

			DoubleVector1D pixel_difference(nMatches_stepI);
			double ave_pixel_difference(0);
			for (int i = 0; i < nMatches_stepI; i++)
			{
				pixel_difference[i] = sqrt(pow(matches_image_stepII_barrycentric[0][i] - Matches_image_stepI[0][i], 2) +
					pow(matches_image_stepII_barrycentric[1][i] - Matches_image_stepI[1][i], 2));
			}

			// Removing potential outliers using Median Absolute Deviation
			DoubleVector1D temp_rho(nMatches_stepI);
			double rho_median(calc_median(pixel_difference));
			for (int i = 0; i < nMatches_stepI; i++)
			{
				temp_rho[i] = abs(pixel_difference[i] - rho_median);
			}
			double mad(1.4826 * calc_median(temp_rho));

			DoubleVector1D deviation_MAD_ratios(nMatches_stepI);
			for (int i = 0; i < nMatches_stepI; i++)
			{
				deviation_MAD_ratios[i] = abs(pixel_difference[i] - rho_median) / mad;
			}

			DoubleVector2D Matches_template_stepII(2);
			DoubleVector2D Matches_image_stepII(2);
			for (int i = 0; i < nMatches_stepI; i++)
			{
				if (deviation_MAD_ratios[i] < crit_MAD)
				{
					Matches_template_stepII[0].push_back(Matches_template_stepI[0][i]);
					Matches_template_stepII[1].push_back(Matches_template_stepI[1][i]);
					Matches_image_stepII[0].push_back(Matches_image_stepI[0][i]);
					Matches_image_stepII[1].push_back(Matches_image_stepI[1][i]);
				}
			}
			int nMatches_stepII(Matches_image_stepII[0].size());

			// Calculating Average distance between Matches in Step II
			int s(0);
			double ave_dist_matches(0);
			for (int i = 1; i < nMatches_stepII; i++)
			{
				for (int j = 0; j < nMatches_stepII; j++)
				{
					if (i != j)
					{
						ave_dist_matches += sqrt(pow(Matches_image_stepII[0][i] - Matches_image_stepII[0][j], 2) +
							pow(Matches_image_stepII[1][i] - Matches_image_stepII[1][j], 2));
						s++;
					}
				}
			}
			ave_dist_matches = ave_dist_matches / s;
			
			if (nMatches_stepII > 0)
			{
				// Calculating Warp W2
				Mesh_transformed_stepII = warp(bbs, Warp_nC, Warp_er, K, Matches_template_stepII, Matches_image_stepII, Template.Mesh_project);
				
				auto stop_stepII = chrono::high_resolution_clock::now();
				
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				// Step III
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				auto start_stepIII = chrono::high_resolution_clock::now();

				// Checking All matches ===========================================================
				// Evaluating Distances using barycentric coordinates 
				DoubleVector2D matches_image_stepIII_barrycentric(2, vector<double>(nMatches));

				for (int i = 0; i < nMatches; i++)
				{
					int match_index(Goodmatches_ids[i]);
					matches_image_stepIII_barrycentric[0][i] = Template.Keypoints_BC_Coefficients[0][match_index] * Mesh_transformed_stepII[0][Template.Tri[Template.Keypoints_BC_related_triangle[match_index]][0]] +
						Template.Keypoints_BC_Coefficients[1][match_index] * Mesh_transformed_stepII[0][Template.Tri[Template.Keypoints_BC_related_triangle[match_index]][1]] +
						Template.Keypoints_BC_Coefficients[2][match_index] * Mesh_transformed_stepII[0][Template.Tri[Template.Keypoints_BC_related_triangle[match_index]][2]];

					matches_image_stepIII_barrycentric[1][i] = Template.Keypoints_BC_Coefficients[0][match_index] * Mesh_transformed_stepII[1][Template.Tri[Template.Keypoints_BC_related_triangle[match_index]][0]] +
						Template.Keypoints_BC_Coefficients[1][match_index] * Mesh_transformed_stepII[1][Template.Tri[Template.Keypoints_BC_related_triangle[match_index]][1]] +
						Template.Keypoints_BC_Coefficients[2][match_index] * Mesh_transformed_stepII[1][Template.Tri[Template.Keypoints_BC_related_triangle[match_index]][2]];
				}

				DoubleVector1D pixel_difference_stepIII(nMatches);
				for (int i = 0; i < nMatches; i++)
				{
					pixel_difference_stepIII[i] = sqrt(pow(matches_image_stepIII_barrycentric[0][i] - Matches_image[0][i], 2) +
						pow(matches_image_stepIII_barrycentric[1][i] - Matches_image[1][i], 2));
				}

				
				int num_selected_inliers_all(0);
				int num_selected_outliers_all(0);
				for (int i = 0; i < nMatches; i++)
				{
					if (pixel_difference_stepIII[i] < alpha_s*ave_dist_matches)
					{
						Matches_template_final[0].push_back(Matches_template[0][i]);
						Matches_template_final[1].push_back(Matches_template[1][i]);
						Matches_image_final[0].push_back(Matches_image[0][i]);
						Matches_image_final[1].push_back(Matches_image[1][i]);
						Matches_template_final_index.push_back(Goodmatches_ids[i]);
					}
				}

				nMatches_final = Matches_template_final_index.size();

				auto stop_stepIII = chrono::high_resolution_clock::now();
				
				if(nMatches_final > 0)
				{
					if(ROBUSfT_plot)
					{
						Mat img_stepI = plot_Grid(Captured_image.clone(), Mesh_transformed_stepI, Template.Nx, Template.Ny, "Outlier removal, Step I", Scalar(0, 135, 255), Scalar(0, 135, 255), 400);
						plot_points(img_stepI, Matches_image_stepI, "Outlier removal, Step I", false, false, Scalar(0, 255, 0), Scalar(255, 0, 0), Scalar(0, 255, 0), 600);	

						Mat img_stepII = plot_Grid(Captured_image.clone(), Mesh_transformed_stepII, Template.Nx, Template.Ny, "Outlier removal, Step II", Scalar(0, 255, 255), Scalar(0, 255, 255), 400);
						plot_points(img_stepII, Matches_image_stepII, "Outlier removal, Step II", false, false, Scalar(0, 255, 0), Scalar(255, 0, 0), Scalar(0, 255, 0), 600);	

						if(Mesh_transformed_stepIII.size() > 0) 
						{
							Mat img_stepIII = plot_Grid(Captured_image.clone(), Mesh_transformed_stepIII, Template.Nx, Template.Ny, "Outlier removal, Step III", Scalar(255, 255, 0), Scalar(255, 255, 0), 400);
							plot_points(img_stepIII, Matches_image_final, "Outlier removal, Step III", false, false, Scalar(0, 255, 0), Scalar(255, 0, 0), Scalar(0, 255, 0), 600);			
						}
					}
				}
				else
				{
					if(ROBUSfT_plot)
					{
						imshow("Outlier removal, Step I", Captured_image);
						imshow("Outlier removal, Step II", Captured_image);
						imshow("Outlier removal, Step III", Captured_image);
					}
				}				
				
				if(ROBUSfT_report) 
				{
					cout << "=======================================================" << endl;
					cout << "Outlier Removal Results" << endl;
					cout << "=======================================================" << endl;
					cout << "Number of Total Matches = " << nMatches << endl; 
					cout << "Number of Final Selected Correct Matches = " << nMatches_final << endl; 
					cout << "Number of Final Selected Wrong Matches = " << nMatches - nMatches_final << endl; 
					cout << "Number of Selected Matches Step I = " << nMatches_stepI << endl; 
					cout << "Number of Selected Matches Step II = " << nMatches_stepII << endl;
					cout << "Execution Times ------------------------" << endl;
					cout << "Step I = " 	<< std::chrono::duration<double>(stop_stepI - start_stepI).count() << " (s)" << endl;
					cout << "Step II = " 	<< std::chrono::duration<double>(stop_stepII - start_stepII).count() << " (s)" << endl;
					cout << "Step III = " 	<< std::chrono::duration<double>(stop_stepIII - start_stepIII).count() << " (s)" << endl;
					cout << "Total Outlier Removal = " 	<< std::chrono::duration<double>(stop_stepIII - start_stepI).count() << " (s)" << endl;
					cout << "=======================================================" << endl;
				}
			
			}
			else
			{
				if(ROBUSfT_plot)
				{
					imshow("Outlier removal, Step I", Captured_image);
					imshow("Outlier removal, Step II", Captured_image);
					imshow("Outlier removal, Step III", Captured_image);
				}
			}		

		}
		else
		{
			if(ROBUSfT_plot)
			{
				imshow("Outlier removal, Step I", Captured_image);
				imshow("Outlier removal, Step II", Captured_image);
				imshow("Outlier removal, Step III", Captured_image);
			}
		}

		if(ROBUSfT_plot)
		{
			moveWindow("Outlier removal, Step I", 100, 530);
			moveWindow("Outlier removal, Step II", 500, 530);
			moveWindow("Outlier removal, Step III", 900, 530);
			waitKey(1);
		}
	}

void Robusft::set_sightlines()
	{
		Sightlines.clear();
		Sightlines.resize(2);
		Sightlines_ids.clear();	

		for (int i = 0; i < nMatches_final; i++)
		{
			// Finding related mesh grid points surrounding the match points to be used as visible points
			int relevant_index(Matches_template_final_index[i]);
			for (int k = 0; k < 3; k++)
			{
				int mesh_num(Template.Keypoints_BC_related_meshpoints[k][relevant_index]);
				if (find_in_vector(Sightlines_ids, mesh_num) == -2)
				{
					Sightlines_ids.push_back(mesh_num);
					Sightlines[0].push_back(Mesh_transformed_stepIII[0][mesh_num]);
					Sightlines[1].push_back(Mesh_transformed_stepIII[1][mesh_num]);
				}
			}
		}

		nSightlines = Sightlines_ids.size();
	}

void Robusft::shapeInference()
	{
		// Initialization (physics & position based model) ---------
		DoubleVector2D particles_pose(Shape_3d);
		DoubleVector2D particles_pos_prev(Shape_3d);
		DoubleVector2D particles_vel(3, vector<double>(Template.nParticles));
		DoubleVector2D particles_force(3, vector<double>(Template.nParticles));
		DoubleVector1D particles_mass(Template.nParticles);
		fill(particles_mass.begin(), particles_mass.end(), Template.Mass / Template.nParticles);

		// Getting Sightlines unit vectors --------------------------
		double fx{K[0][0]};
		double fy{K[1][1]};
		double u0{K[0][2]};
		double v0{K[1][2]};

		DoubleVector2D xy(3, vector<double>(nSightlines));
		DoubleVector1D xy_norm(nSightlines);
		DoubleVector2D Sightlines_unit(3, vector<double>(nSightlines));

		for (int i = 0; i < nSightlines; i++)
		{
			xy[0][i] = (Sightlines[0][i] - u0) / fx;
			xy[1][i] = (Sightlines[1][i] - v0) / fy;
			xy[2][i] = 1;
			xy_norm[i] = sqrt(pow(xy[0][i], 2) + pow(xy[1][i], 2) + pow(xy[2][i], 2));
			Sightlines_unit[0][i] = xy[0][i] / xy_norm[i];
			Sightlines_unit[1][i] = xy[1][i] / xy_norm[i];
			Sightlines_unit[2][i] = xy[2][i] / xy_norm[i];
		}

		// Loop ---------------------------------------------------------------------------------
		double ts{1};
		for (int iter = 0; iter < SfT_maxIter; iter++)
		{

			// Damp velocities & estimate new positions -----------------------------------------
			for (int j = 0; j < Template.nParticles; j++)
			{
				for (int i = 0; i < 3; i++)
				{
					particles_vel[i][j] = particles_vel[i][j]*(1 - Template.K_damping);
					particles_pose[i][j] = Shape_3d[i][j] + ts*particles_vel[i][j];
				}
			}

			// Apply Sight-line constraints -----------------------------------------------------
			for (int i = 0; i < nSightlines; i++)
			{
				int j{ static_cast<int>(Sightlines_ids[i]) };
				double sum_temp{ particles_pose[0][j] * Sightlines_unit[0][i] + particles_pose[1][j] * Sightlines_unit[1][i] + particles_pose[2][j] * Sightlines_unit[2][i] };
				if (sum_temp < 0) sum_temp = -1 * sum_temp;
				particles_pose[0][j] = sum_temp * Sightlines_unit[0][i];
				particles_pose[1][j] = sum_temp * Sightlines_unit[1][i];
				particles_pose[2][j] = sum_temp * Sightlines_unit[2][i];
			}

			// Solver ---------------------------------------------------------------------------
			for (int solverIter = 0; solverIter < SfT_nSolverIterations; solverIter++)
			{
				// apply stretching constraints -------------------------------------------------
				for (int j = 0; j < Template.nSGroups; j++)
				{

					DoubleVector2D group_edges(Template.SGroups[j]);
					DoubleVector1D length0(Template.Sgroups_length0[j]);
					DoubleVector1D k_stiffness(Template.Sgroups_k[j]);

					DoubleVector1D i1(transpose(group_edges)[0]);
					DoubleVector1D i2(transpose(group_edges)[1]);
					int i1_size{static_cast<int>(i1.size())};
					DoubleVector1D d(3);
					double d_norm;
					DoubleVector1D d_unit(3);
					DoubleVector1D e(3);
					double coef_temp{0.0};
					
					for (int i = 0; i < i1_size; i++)
					{
						d[0] = particles_pose[0][i2[i]] - particles_pose[0][i1[i]];
						d[1] = particles_pose[1][i2[i]] - particles_pose[1][i1[i]];
						d[2] = particles_pose[2][i2[i]] - particles_pose[2][i1[i]];

						d_norm = sqrt(pow(d[0],2)+pow(d[1],2)+pow(d[2],2));

						d_unit[0] = d[0]/d_norm;
						d_unit[1] = d[1]/d_norm;
						d_unit[2] = d[2]/d_norm;

						coef_temp = 0.5*k_stiffness[i]*(d_norm - length0[i]);
						e[0] = coef_temp*d_unit[0];
						e[1] = coef_temp*d_unit[1];
						e[2] = coef_temp*d_unit[2];

						particles_pose[0][i1[i]] += e[0];
						particles_pose[1][i1[i]] += e[1];
						particles_pose[2][i1[i]] += e[2];

						particles_pose[0][i2[i]] -= e[0];
						particles_pose[1][i2[i]] -= e[1];
						particles_pose[2][i2[i]] -= e[2];
					}
				}

				// apply bending constraints -------------------------------------------------
				for (int j = 0; j < Template.nBGroups; j++)
				{

					DoubleVector2D group_edges(Template.BGroups[j]);
					DoubleVector1D length0(Template.Bgroups_length0[j]);
					DoubleVector1D k_stiffness(Template.Bgroups_k[j]);

					DoubleVector1D i1(transpose(group_edges)[0]);
					DoubleVector1D i2(transpose(group_edges)[1]);
					int i1_size{static_cast<int>(i1.size())};
					DoubleVector1D d(3);
					double d_norm;
					DoubleVector1D d_unit(3);
					DoubleVector1D e(3);
					double coef_temp{0.0};
					
					for (int i = 0; i < i1_size; i++)
					{
						d[0] = particles_pose[0][i2[i]] - particles_pose[0][i1[i]];
						d[1] = particles_pose[1][i2[i]] - particles_pose[1][i1[i]];
						d[2] = particles_pose[2][i2[i]] - particles_pose[2][i1[i]];

						d_norm = sqrt(pow(d[0],2)+pow(d[1],2)+pow(d[2],2));

						d_unit[0] = d[0]/d_norm;
						d_unit[1] = d[1]/d_norm;
						d_unit[2] = d[2]/d_norm;

						coef_temp = 0.5*k_stiffness[i]*(d_norm - length0[i]);
						e[0] = coef_temp*d_unit[0];
						e[1] = coef_temp*d_unit[1];
						e[2] = coef_temp*d_unit[2];

						particles_pose[0][i1[i]] += e[0];
						particles_pose[1][i1[i]] += e[1];
						particles_pose[2][i1[i]] += e[2];

						particles_pose[0][i2[i]] -= e[0];
						particles_pose[1][i2[i]] -= e[1];
						particles_pose[2][i2[i]] -= e[2];
					}
				}

				// Applying known Particles ----------------------------------------------------
				if (KnownParticles_ids.size() > 0)
				{
					for (int j = 0; j < KnownParticles_ids.size(); j++)
					{
						particles_pose[0][KnownParticles_ids[j]] = KnownParticles_3D[0][j];
						particles_pose[1][KnownParticles_ids[j]] = KnownParticles_3D[1][j];
						particles_pose[2][KnownParticles_ids[j]] = KnownParticles_3D[2][j];
					}
				}
			}

			// Update Particle States ----------------------------------------------------------
			particles_pos_prev = Shape_3d;
			for (int j = 0; j < Template.nParticles; j++)
			{
				for (int i = 0; i < 3; i++)
				{
					particles_vel[i][j] = (particles_pose[i][j] - Shape_3d[i][j])/ts;
				}
			}
			Shape_3d = particles_pose;

			// break loop condition (while the mesh moves no more than the threshold RMSE value)
			double criterion(0);
			for (int j = 0; j < Template.nParticles; j++)
			{
				double power(0);
				for (int i = 0; i < 3; i++)
				{
					power += pow((particles_pos_prev[i][j] - Shape_3d[i][j]),2);
				}
				criterion += power;
			}
			criterion = sqrt(criterion/Template.nParticles);
			if(criterion < SfT_threshold)
			{
				Shape_2d = projection(K, Shape_3d);
				break;
			} 
		}

		// Positioning the 3D shape in front of the camera if it is formed behind the camera -----
		if(Shape_3d[2][0] < 0)
		{
			for (int j = 0; j < Template.nParticles; j++)
			{
				for (int i = 0; i < 3; i++)
				{
					Shape_3d[i][j] = -1*Shape_3d[i][j];
				}
			}
		}

	}

void Robusft::resetVariables()
	{
		Matches_template_final.clear();
		Matches_template_final.resize(2);
		Matches_image_final.clear();
		Matches_image_final.resize(2);
		Matches_template_final_index.clear();
		nMatches_final = 0;
		nMatches = 0;
	}

void Robusft::pipeline_CPU()
{
	if(ROBUSfT_plot)
	{
		namedWindow("Original", cv::WINDOW_NORMAL);
		resizeWindow("Original", 400, 400);
		namedWindow("3D Projection", cv::WINDOW_NORMAL);
		resizeWindow("3D Projection", 400, 400);
		namedWindow("3D", cv::WINDOW_NORMAL);
		resizeWindow("3D", 400, 400);
	}

	while(true)
	{
		
		resetVariables();

		auto start_total = chrono::high_resolution_clock::now();

		unique_lock<mutex> ul(my_mutex);	
		con_var.wait(ul, [] { return (job_done) ? true : false; });


		if(ROBUSfT_report) cout << "#######################################################" << endl;

		// ======================================================================================================
		// Matching
		// ======================================================================================================
		auto start_matching = chrono::high_resolution_clock::now();

		match();

		my_mutex.unlock();	

		auto stop_matching = chrono::high_resolution_clock::now();
		if(ROBUSfT_report) cout << "Matching CPU Thread = " << std::chrono::duration<double>(stop_matching - start_matching).count() << " (s)" << endl;
		
		if(nMatches > 3)
		{
			//=====================================================================================
			//	Outlier Removal Algorithm
			//=====================================================================================
			auto start_outlier_removal = chrono::high_resolution_clock::now();
			
			outlier_removal_algorithm();

			auto stop_outlier_removal = chrono::high_resolution_clock::now();
			if(ROBUSfT_report) cout << "Outlier Removal CPU Thread = " << std::chrono::duration<double>(stop_outlier_removal - start_outlier_removal).count() << " (s)" << endl;

			if(nMatches_final > 0)
			{
				//================================================================================
				// Warp
				//================================================================================
				auto start_Warp = chrono::high_resolution_clock::now();

				// Calculating Warp W3
				Mesh_transformed_stepIII = warp(bbs, Warp_nC, Warp_er, K, Matches_template_final, Matches_image_final, Template.Mesh_project);
				
				auto stop_Warp = chrono::high_resolution_clock::now();
				if(ROBUSfT_report) cout << "Final Warp CPU Thread = " << std::chrono::duration<double>(stop_Warp - start_Warp).count() << " (s)" << endl;

				//================================================================================
				// Shape Inference
				//================================================================================
				auto start_Inference = chrono::high_resolution_clock::now();
				
				// Setting Sightlines for Particle-SfT
				set_sightlines();

				// Setting 3D positions of Known particles for Particle-SfT, Here is an example
				// softObject.KnownParticles_ids = {3, 54};										// Id of mesh points that their 3D positions in the camera frame are known					
				// softObject.KnownParticles_3D = {{0.32, 0.25, 1.07}, {0.65, 0.29, 1.12}};		// 3D positions of those nodes in the camera frame	

				shapeInference();

				auto stop_Inference = chrono::high_resolution_clock::now();
				if(ROBUSfT_report) cout << "Shape Inference CPU Thread = " << std::chrono::duration<double>(stop_Inference - start_Inference).count() << " (s)" << endl;

				auto stop_total = chrono::high_resolution_clock::now();
				if(ROBUSfT_report) cout << "Total CPU Thread = " << std::chrono::duration<double>(stop_total - start_total).count() << " (s)" << endl;
				

				//================================================================================
				// Plotting
				//================================================================================
				if(ROBUSfT_plot)
				{
					// Original Captured Image
					imshow("Original", Captured_image);

					// Projection of the 3D shape and sightlines
					Mat img_SfT = plot_Grid(Captured_image.clone(), Shape_2d, Template.Nx, Template.Ny, "3D Projection", Scalar(0, 255, 0), Scalar(0, 255, 0), 400);
					plot_points(img_SfT, Sightlines, "3D Projection", false, false, Scalar(0, 0, 255), Scalar(255, 0, 0), Scalar(0, 255, 0), 600);			

					// 3D shape
					DoubleVector1D Camera_location = {0.8, -0.5, -1}; 
					DoubleVector1D Focus_point = {0, 0, 1};
					plot_3D(Shape_3d, Template.Tri, "3D", "point", Camera_location, Focus_point);
				}
				
			}
			else
			{
				if(ROBUSfT_plot)
				{
					imshow("Original", Captured_image);
					imshow("3D Projection", Captured_image);
				}
			}
		}
		else
		{
			if(ROBUSfT_plot)
			{
				imshow("Original", Captured_image);
				imshow("3D Projection", Captured_image);
			}
		}

		if(ROBUSfT_plot)
		{
			moveWindow("Original", 100, 100);
			moveWindow("3D Projection", 500, 100);
			moveWindow("3D", 900, 100);
			waitKey(1);
		}

	}
}

void Robusft::online_execution(VideoCapture cap)
{
	// Threads =================================================
	std::thread t1(&Robusft::extract_keypoints_GPU, this, "Image", cap);
	std::thread t2(&Robusft::pipeline_CPU, this);

	t1.join();
	t2.join();
}