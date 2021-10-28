#include <iostream>
#include "ROBUSfT/headers.h"
#include "ROBUSfT/delaunator_h.hpp"
#include "ROBUSfT/functions.h"
#include "ROBUSfT/bbs_eigen.h" 
#include "ROBUSfT/ROBUSfT.h"

#if POPSIFT_IS_DEFINED(POPSIFT_USE_NVTX)
#include <nvToolsExtCuda.h>
#else
#define nvtxRangePushA(a)
#define nvtxRangePop()
#endif

bool ROBUSfT_report = false;
bool ROBUSfT_plot = true; 

int main()
{
	//#############################################################################
	// OFFLINE SECTION
	//#############################################################################

	cudaDeviceReset();
	popsift::Config config;
	PopSift PopSift(config, popsift::Config::ExtractingMode);

	Robusft softObject;

	// Calibration Matrix ======================================
	DoubleVector2D Calibration_Matrix{ 	{817.6721, 	0, 			330.7831},
										{0, 		821.8128, 	235.4349},
										{0,			0,			1} };
	softObject.setCalibrationMat(Calibration_Matrix);

	// Template ================================================
	softObject.Template.Mesh_type = "irregular";	// Type of mesh, regular for rectangular shapes and irregular for nonrectangular shapes
	string address = string(ROOT_FOLDER) + string("/Texturemaps/shoe.jpg");
	softObject.set_template_image(address);
	
	// Input for Creating an Irregular Mesh --------------------
	// softObject.Template.Nx & Ny will be set automatically equal to nParticles and 1, respectively.
	softObject.Template.nGrid_on_boundary = 20;
	softObject.Template.nGrid_inside_boundary = 20;

	softObject.Template.Width = 0.30;		// Width of the object in X direction, It is 300mm for the shoe sole
	softObject.Template.Mass = 0.01;		// Will be used if the effect of the gravity is considered
	softObject.Template.K_damping = 0.008;	// Will be used for relaxing velocities id Particle-SfT
	softObject.Template.K_stretch = 1;		// Should be 1 for isometric deformation 
	softObject.Template.K_bend = 1;			// Should be 1 for isometric deformation 
	bool load_template = false;
	softObject.build_template(load_template);

	// Initial Pose of the Deformable Object ===================
	softObject.Shape_3d = softObject.Template.Mesh;				
	softObject.Shape_2d = softObject.Template.Mesh_project;

	// Video Configuration =====================================
	VideoCapture cap(0); //capture the video from web cam
	cap.set(cv::CAP_PROP_FRAME_WIDTH, 640); 
	cap.set(cv::CAP_PROP_FRAME_HEIGHT, 480);
	
	// Warp Parameters =========================================
	softObject.Warp_nC = 12;
	softObject.Warp_er = 500;
	softObject.bbs = {softObject.Template.Mesh_project_min_x-20, softObject.Template.Mesh_project_max_x+20, softObject.Warp_nC, softObject.Template.Mesh_project_min_y-20, softObject.Template.Mesh_project_max_y+20, softObject.Warp_nC, 2 };

	// SfT Parameters ==========================================
	softObject.SfT_maxIter = 2000;
	softObject.SfT_nSolverIterations = 2;
	softObject.SfT_threshold = 0.00001;

	//#############################################################################
	// ONLINE SECTION
	//#############################################################################

	softObject.online_execution(cap);
}