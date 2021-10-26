/* bbs.h
 *
 * General routines for cubic bidimensional b-splines.
 * 
 * History
 *  2009/??/??: First version
 *  2010/12/15: Translation to CPP
 *              Adding OpenMP support
 *
 * (c)2009-2010, Florent Brunet.
 */
 
 /*
 * BBS is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * BBS is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
 
#ifndef __BBS_EIGEN_H__
#define __BBS_EIGEN_H__ 1

//#include "Engine.h"
// #include "matrix.h"
#include "headers.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/SVD>
//#include <Eigen/SparseExtra>
#include<Eigen/SparseCholesky>
#include<Eigen/SparseQR>
#include <type_traits>
#include <limits>
#include <Eigen/IterativeLinearSolvers> 

using namespace Eigen;
typedef SparseMatrix<double, ColMajor> MatlabSparse;

/* The B-Spline structure */
typedef struct _bbs_t {
    double umin;
    double umax;
    int nptsu;
    double vmin;
    double vmax;
    int nptsv;
    int valdim;
} bbs_t;

void normalize(double xmin, double xmax, int npts, double *x, int nb_x, double* nx);
void normalize_with_inter(double xmin, double xmax, int npts, double *x, int nb_x, double* nx, int *inter);

void eval_basis(double nx, double *val_basis);
void eval_basis_d(double nx, double *val_basis);
void eval_basis_dd(double nx, double *val_basis);

// void eval(bbs_t *bbs, double *ctrlpts, double *u, double *v, int nval, double *val, int du, int dv);
void eval_new(bbs_t* bbs, double* ctrlpts, double* u, double* v, int nval, double* val, int du, int dv);

// int coloc(bbs_t *bbs, double *u, double *v, int nsites, double *pr, mwIndex *ir, mwIndex *jc);
int coloc_new(bbs_t* bbs, double* u, double* v, int nsites, double* pr, MatlabSparse::StorageIndex* ir, MatlabSparse::StorageIndex* jc);
// int coloc_deriv(bbs_t *bbs, double *u, double *v, int nsites, int du, int dv, double *pr, mwIndex *ir, mwIndex *jc);

// void bending_ur(bbs_t *bbs, double* lambdas, double *pr, mwIndex *ir, mwIndex *jc);
void bending_ur_new(bbs_t* bbs, double* lambdas, double* pr, MatlabSparse::StorageIndex* ir, MatlabSparse::StorageIndex* jc);

// Map<MatlabSparse> matlab_to_eigen_sparse(const mxArray* mat);
// mxArray* eigen_to_matlab_sparse(const Ref<const MatlabSparse, StandardCompressedFormat>& mat);

// void my_coloc(int nlhs, mxArray* plhs[], int nrhs, bbs_t bbs, mxArray* prhs2, mxArray* prhs3);
void my_coloc_new(int nlhs, MatlabSparse &plhs, int nrhs, bbs_t bbs, VectorXd prhs2, VectorXd prhs3);
// mxArray* my_bbs_coloc(bbs_t bbs, DoubleVector1D p_0, DoubleVector1D p_1);
MatlabSparse my_bbs_coloc_new(bbs_t bbs, DoubleVector1D p_0, DoubleVector1D p_1);
// void my_bending(int nlhs, mxArray* plhs[], bbs_t bbs, mxArray* prhs);
void my_bending_new(int nlhs, MatlabSparse &plhs, bbs_t bbs, MatrixXd prhs);
// mxArray* my_bbs_bending(bbs_t bbs, DoubleVector2D lambdas);
MatlabSparse my_bbs_bending_new(bbs_t bbs, DoubleVector2D lambdas);
// void my_eval(int nlhs, mxArray* plhs[], int nrhs, bbs_t bbs, mxArray* prhs4, mxArray* prhs2, mxArray* prhs3, int du, int dv);
void my_eval_new(int nlhs, MatlabSparse &plhs, bbs_t bbs, MatrixXd prhs1, VectorXd prhs2, VectorXd prhs3, int du, int dv);
// mxArray* my_bbs_eval(bbs_t bbs, DoubleVector2D ctrlpts, DoubleVector1D p_0, DoubleVector1D p_1, int du, int dv);
MatlabSparse my_bbs_eval_new(bbs_t bbs, DoubleVector2D ctrlpts, DoubleVector1D p_0, DoubleVector1D p_1, int du, int dv);
DoubleVector2D BBS_Function(bbs_t bbs, DoubleVector2D pixles_template, DoubleVector2D pixles_image, DoubleVector2D K, int nC, double er);
DoubleVector2D BBS_Evaluation(bbs_t bbs, DoubleVector2D ctrlpts, DoubleVector2D pixles, DoubleVector2D K, int deriv_1, int deriv_2);


#endif
