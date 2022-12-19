
#ifndef __ROBUSFT_H__
#define __ROBUSFT_H__ 

#include <iostream>
#include "headers.h"
#include "delaunator_h.hpp"
#include "functions.h"
#include "bbs_eigen.h" 

extern mutex my_mutex;
extern bool job_done;
extern condition_variable con_var;
extern bool ROBUSfT_report;
extern bool ROBUSfT_plot;

class Robusft
{
public:
	
	struct TemplateVariables
	{
		DoubleVector3D SGroups;
		DoubleVector2D Sgroups_length0;
		DoubleVector2D Sgroups_k;
		int nSGroups;
		DoubleVector3D BGroups;
		DoubleVector2D Bgroups_length0;
		DoubleVector2D Bgroups_k;
		int nBGroups;
		DoubleVector2D Vertices;
		DoubleVector2D Tri;
		double Width;
		string Mesh_type;
		int Nx;
		int Ny;
		double Sx;
		double Sy;
		int nParticles;
		int nGrid_on_boundary;
		int nGrid_inside_boundary;
		double Mesh_min_x;
		double Mesh_max_x;
		double Mesh_min_y;
		double Mesh_max_y;
		double Mesh_project_min_x;
		double Mesh_project_max_x;
		double Mesh_project_min_y;
		double Mesh_project_max_y;
		double Mesh_z;
		double Mass;
		double K_damping;
		double K_stretch;
		double K_bend;
		DoubleVector2D Mesh;
		DoubleVector2D Mesh_project;
		DoubleVector1D Keypoints_BC_related_triangle;
		DoubleVector2D Keypoints_BC_Coefficients;
		DoubleVector2D Keypoints_BC_related_meshpoints;
		int nKeypoints;
		DoubleVector2D Features;
		Mat Descriptors;
		
	};

	TemplateVariables Template;
	Mat Template_image;
	Mat Template_image_BW;
	Mat Captured_image;
	Mat Captured_image_BW;
	DoubleVector2D K;

	// Matching ================================================
	int nKeypoints;
	DoubleVector2D Features;
	Mat Descriptors;
	DoubleVector1D Goodmatches_ids;
	DoubleVector2D Matches_template;
	DoubleVector2D Matches_image;
	int nMatches;

	// outlier removal =========================================
	DoubleVector2D Mesh_transformed_stepI;
	DoubleVector2D Mesh_transformed_stepII;
	DoubleVector2D Mesh_transformed_stepIII;
	DoubleVector2D Matches_template_final;
	DoubleVector2D Matches_image_final;
	DoubleVector1D Matches_template_final_index;
	int nMatches_final;

	// Warp ====================================================
	int Warp_nC;
	double Warp_er;
	bbs_t bbs;

	// SfT =====================================================
	int SfT_maxIter;
	int SfT_nSolverIterations;
	double SfT_threshold;
	DoubleVector2D Shape_2d;
	DoubleVector2D Shape_3d;
	int nSightlines;
	DoubleVector1D Sightlines_ids;
	DoubleVector2D Sightlines;
	DoubleVector1D KnownParticles_ids;
	DoubleVector2D KnownParticles_3D;

	void set_template_image(string address);
	void setCalibrationMat(DoubleVector2D calibrationMat);
	void build_template(bool load_template = false);
	void extract_keypoints_GPU(string type, VideoCapture cap = VideoCapture(0));
	void get_BaryCentricCoordsOfKeypoints();
	void capture_image(VideoCapture cap);
	void match(float ratio_thresh = 0.8f);
	void outlier_removal_algorithm(int crit_nMatches_stepI = 50, double crit_MAD = 2.5, double alpha_s = 0.07);
	void calculate_warp();
	void set_sightlines();
	void shapeInference();
	void resetVariables();
	void pipeline_CPU();
	void online_execution(VideoCapture cap = VideoCapture(0));
};

#endif