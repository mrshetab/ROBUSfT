/* bbs.cpp
 *
 * General routines for cubic bidimensional b-splines.
 * 
 * History
 *  2009/??/??: First version
 *  2010/12/15: Translation to CPP
 *              Adding OpenMP support
 *  2011/02/03: Defining number of thread by compilation directive
 *
 * (c)2009-2011, Florent Brunet.
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

#include "bbs_eigen.h"
#include <stdlib.h>
#include <math.h>
#include "headers.h"
#include "functions.h"
#include <omp.h>

#ifndef NTHREADS
#define NTHREADS 1
#endif

#define maxx(a,b) (((a)>=(b)) ? (a) : (b))
#define min(a,b) (((a)<=(b)) ? (a) : (b))

 /* Normalize the values in "x". "nb_x" is the size of the "x" and "nx" arrays. */
void normalize(double xmin, double xmax, int npts, double* x, int nb_x, double* nx) {
    int ninter = npts - 3;
    double width_inter = (xmax - xmin) / ninter;

#pragma omp parallel for num_threads(NTHREADS), schedule(guided)
    for (int i = 0; i < nb_x; ++i) {
        if (x[i] == xmax)
            nx[i] = 1.0;
        else if (x[i] < xmin)
            nx[i] = (x[i] - xmin) / width_inter;
        else if (x[i] > xmax)
            nx[i] = (x[i] - xmin) / width_inter - ninter;
        else {
            double scaled = (x[i] - xmin) / width_inter;
            nx[i] = scaled - floor(scaled);
        }
    }
}

/* Same as "normalize" but computes also the interval number to which the x's belongs. */
void normalize_with_inter(double xmin, double xmax, int npts, double* x, int nb_x, double* nx, int* inter) {
    int ninter = npts - 3;
    double width_inter = (xmax - xmin) / ninter;

#pragma omp parallel for num_threads(NTHREADS), schedule(guided)
    for (int i = 0; i < nb_x; ++i) {
        if (x[i] == xmax) {
            nx[i] = 1.0;
            inter[i] = ninter - 1;
        }
        else if (x[i] < xmin) {
            nx[i] = (x[i] - xmin) / width_inter;
            inter[i] = -1;
        }
        else if (x[i] > xmax) {
            nx[i] = (x[i] - xmin) / width_inter - ninter;
            inter[i] = ninter;
        }
        else {
            double scaled = (x[i] - xmin) / width_inter;
            inter[i] = (int)floor(scaled);
            nx[i] = scaled - inter[i];
        }
    }
}

/* Evaluate the 4 basis at nx. "nx" must be normalized values. */
void eval_basis(double nx, double* val_basis) {
    double nx2 = nx * nx;
    double nx3 = nx2 * nx;
    val_basis[0] = (-nx3 + 3.0 * nx2 - 3.0 * nx + 1.0) / 6.0;
    val_basis[1] = (3.0 * nx3 - 6.0 * nx2 + 4.0) / 6.0;
    val_basis[2] = (-3.0 * nx3 + 3.0 * nx2 + 3.0 * nx + 1.0) / 6.0;
    val_basis[3] = nx3 / 6.0;
}

/* Evaluate the 4 basis first derivatives at nx. "nx" must be normalized values. */
void eval_basis_d(double nx, double* val_basis) {
    double nx2 = nx * nx;
    val_basis[0] = (-nx2 + 2 * nx - 1) / 2.0;
    val_basis[1] = (3.0 * nx2 - 4.0 * nx) / 2.0;
    val_basis[2] = (-3 * nx2 + 2 * nx + 1) / 2.0;
    val_basis[3] = nx2 / 2.0;
}

/* Evaluate the 4 basis second derivatives at nx. "nx" must be normalized values. */
void eval_basis_dd(double nx, double* val_basis) {
    val_basis[0] = -nx + 1.0;
    val_basis[1] = 3.0 * nx - 2.0;
    val_basis[2] = -3.0 * nx + 1.0;
    val_basis[3] = nx;
}

/* Pointer type for a basis function evaluation */
typedef void (*basis_func_t)(double, double*);

/* Return a pointer to the basis function that correspond to the "order" of derivation. */
basis_func_t get_basis_ptr(int order) {
    switch (order) {
    case 0: return eval_basis;
    case 1: return eval_basis_d;
    case 2: return eval_basis_dd;
    }
    return NULL;
}

double get_deriv_fact(bbs_t* bbs, int uorder, int vorder) {
    double su = (bbs->umax - bbs->umin) / (bbs->nptsu - 3);
    double sv = (bbs->vmax - bbs->vmin) / (bbs->nptsv - 3);

    return 1.0 / (pow(su, uorder) * pow(sv, vorder));
}

/* Evaluate a spline (or its derivatives, according to the basis functions utilized).
 * - bbs: bidimensional bspline structure.
 * - ctrlpts: control points (valdim x npts matrix).
 * - u, v: locations where the b-spline is evaluated.
 * - nval: number of values in u (and v).
 * - val: valdim x npts matrix containing the result.
 * - du, dv: order of derivation. */
void eval_new(bbs_t* bbs, double* ctrlpts, double* u, double* v, int nval, double* val, int du, int dv) {
    double* nu = (double*)malloc(nval * sizeof(double));
    double* nv = (double*)malloc(nval * sizeof(double));
    int* interu = (int*)malloc(nval * sizeof(int));
    int* interv = (int*)malloc(nval * sizeof(int));

    // Compute the normalized evaluation values and their interval numbers
    normalize_with_inter(bbs->umin, bbs->umax, bbs->nptsu, u, nval, nu, interu);
    normalize_with_inter(bbs->vmin, bbs->vmax, bbs->nptsv, v, nval, nv, interv);

#pragma omp parallel for num_threads(NTHREADS), schedule(guided)
    for (int k = 0; k < nval; ++k) {
        int iu, iv, d, ind;
        double basis_u[4];
        double basis_v[4];
        double bu, bas;
        basis_func_t b_func_u = get_basis_ptr(du);
        basis_func_t b_func_v = get_basis_ptr(dv);
        double fact = get_deriv_fact(bbs, du, dv);

        b_func_u(nu[k], basis_u);
        b_func_v(nv[k], basis_v);

        for (d = 0; d < bbs->valdim; ++d)
            val[bbs->valdim * k + d] = 0.0;    

        for (iu = 0; iu < 4; ++iu) {
            bu = basis_u[iu];
            for (iv = 0; iv < 4; ++iv) {
                bas = bu * basis_v[iv];
                ind = bbs->valdim * ((iu + interu[k]) * bbs->nptsv + iv + interv[k]);
                //mexPrintf("iu=%d iv=%d ind=%f\n", iu, iv, ind);
                for (d = 0; d < bbs->valdim; ++d)
                {
                    //mexPrintf("ind=%d ctrlpts=%f bas=%f index=%d \n", ind, ctrlpts[ind], bas, bbs->valdim * k + d);
                    //mexPrintf("iu=%d iv=%d d=%d index=%d val=%f\n", iu, iv, d, bbs->valdim * k + d, val[bbs->valdim * k + d]);
                    val[bbs->valdim * k + d] += ctrlpts[ind++] * bas;
                    //mexPrintf("iu=%d iv=%d d=%d index=%d val=%f ==== \n", iu, iv, d, bbs->valdim * k + d, val[bbs->valdim * k + d]);
                    //mexPrintf("iu=%d iv=%d d=%d val=%f ==== \n", iu, iv, d, val[bbs->valdim * k + d]);
                }
            }
        }

        for (d = 0; d < bbs->valdim; ++d)
            val[bbs->valdim * k + d] *= fact;
    }

    free(nu);
    free(nv);
    free(interu);
    free(interv);
}

/* Build a sparse colocation matrix.
 * This function does not handle sites that are outside of the domain (an error is issued in such cases).
 * ARGUMENTS:
 * - bbs: b-spline structure.
 * - u,v: sites used to build the colocation matrix.
 * - nsites: number of sites (number of elements in u and v).
 * - pr: array containing the values of the matrix (size: number of non-zero elements, i.e. 16*nsites).
 * - ir: ir(i) is the row index (0-based) of pr(i). The sizes of ir and pr are equal.
 * - jc: jc(j) is the shift in pr (and ir) for the j-th column. numel(jc) = 1+nptsx*nptsy (number of columns).
 * RETURN VALUE:
 * This function return 0 if successful, an error code otherwise:
 * 0: successful.
 * 1: a colocation site was ouside of the spline definition domain. */
int coloc_new(bbs_t* bbs, double* u, double* v, int nsites, double* pr, MatlabSparse::StorageIndex* ir, MatlabSparse::StorageIndex* jc) {
    int k, iu, iv, col, Iu, Iv;
    int ret_code = 0;
    double* nu = (double*)malloc(nsites * sizeof(double));
    double* nv = (double*)malloc(nsites * sizeof(double));
    int* interu = (int*)malloc(nsites * sizeof(int));
    int* interv = (int*)malloc(nsites * sizeof(int));
    double basis_u[4];
    double basis_v[4];
    int npts = bbs->nptsu * bbs->nptsv;
    int* nb_elem_col = (int*)calloc(npts, sizeof(int));
    int* cur_ind = (int*)malloc(npts * sizeof(int));

    // Compute the normalized evaluation values and their interval numbers
    normalize_with_inter(bbs->umin, bbs->umax, bbs->nptsu, u, nsites, nu, interu);
    normalize_with_inter(bbs->vmin, bbs->vmax, bbs->nptsv, v, nsites, nv, interv);

    // FILLING jc
    // First, compute the number of elements in each column
    // Check at the same time if the sites are inside the definition domain.
    for (k = 0; k < nsites; ++k) {
        Iu = interu[k];
        Iv = interv[k];

        if (Iu < 0 || Iu > bbs->nptsu - 4 || Iv < 0 || Iv > bbs->nptsv - 4) {
            ret_code = 1;
            goto error;
        }

        for (iu = 0; iu < 4; ++iu)
            for (iv = 0; iv < 4; ++iv)
                nb_elem_col[(iu + Iu) * bbs->nptsv + iv + Iv] += 1;
    }

    // Second, compute jc from nb_elem_col
    for (k = 1; k <= npts; ++k)
        jc[k] = jc[k - 1] + nb_elem_col[k - 1];

    // FILLING pr and ir
    // cur_ind contains the current shift in pr (and ir) for each column
    for (k = 0; k < npts; ++k)
        cur_ind[k] = jc[k];

    for (k = 0; k < nsites; ++k) {
        eval_basis(nu[k], basis_u);
        eval_basis(nv[k], basis_v);
        for (iu = 0; iu < 4; ++iu)
            for (iv = 0; iv < 4; ++iv) {
                col = (iu + interu[k]) * bbs->nptsv + iv + interv[k];
                pr[cur_ind[col]] = basis_u[iu] * basis_v[iv];
                ir[cur_ind[col]] = k;
                cur_ind[col] += 1;
            }
    }

error:

    free(nu);
    free(nv);
    free(interu);
    free(interv);
    free(nb_elem_col);
    free(cur_ind);

    return ret_code;
}

/* Build a sparse colocation matrix that accounts for derivatives.
 * This function does not handle sites that are outside of the domain (an error is issued in such cases).
 * ARGUMENTS:
 * - bbs: b-spline structure.
 * - u,v: sites used to build the colocation matrix.
 * - nsites: number of sites (number of elements in u and v).
 * - du,dv: order of derivation along the u- and the v- axis resp.
 * - pr: array containing the values of the matrix (size: number of non-zero elements, i.e. 16*nsites).
 * - ir: ir(i) is the row index (0-based) of pr(i). The sizes of ir and pr are equal.
 * - jc: jc(j) is the shift in pr (and ir) for the j-th column. numel(jc) = 1+nptsx*nptsy (number of columns).
 * RETURN VALUE:
 * This function return 0 if successful, an error code otherwise:
 * 0: successful.
 * 1: a colocation site was ouside of the spline definition domain. */
// int coloc_deriv(bbs_t* bbs, double* u, double* v, int nsites, int du, int dv, double* pr, mwIndex* ir, mwIndex* jc) {
//     int k, iu, iv, col, Iu, Iv;
//     int ret_code = 0;
//     double* nu = (double*)malloc(nsites * sizeof(double));
//     double* nv = (double*)malloc(nsites * sizeof(double));
//     int* interu = (int*)malloc(nsites * sizeof(int));
//     int* interv = (int*)malloc(nsites * sizeof(int));
//     basis_func_t b_func_u = get_basis_ptr(du);
//     basis_func_t b_func_v = get_basis_ptr(dv);
//     double fact = get_deriv_fact(bbs, du, dv);
//     double basis_u[4];
//     double basis_v[4];
//     int npts = bbs->nptsu * bbs->nptsv;
//     mwIndex* nb_elem_col = (mwIndex*)calloc(npts, sizeof(mwIndex));
//     mwIndex* cur_ind = (mwIndex*)malloc(npts * sizeof(mwIndex));
//     // Compute the normalized evaluation values and their interval numbers
//     normalize_with_inter(bbs->umin, bbs->umax, bbs->nptsu, u, nsites, nu, interu);
//     normalize_with_inter(bbs->vmin, bbs->vmax, bbs->nptsv, v, nsites, nv, interv);
//     // FILLING jc
//     // First, compute the number of elements in each column
//     // Check at the same time if the sites are inside the definition domain.
//     for (k = 0; k < nsites; ++k) {
//         Iu = interu[k];
//         Iv = interv[k];
//         if (Iu < 0 || Iu > bbs->nptsu - 4 || Iv < 0 || Iv > bbs->nptsv - 4) {
//             ret_code = 1;
//             goto error;
//         }
//         for (iu = 0; iu < 4; ++iu)
//             for (iv = 0; iv < 4; ++iv)
//                 nb_elem_col[(iu + Iu) * bbs->nptsv + iv + Iv] += 1;
//     }
//     // Second, compute jc from nb_elem_col
//     for (k = 1; k <= npts; ++k)
//         jc[k] = jc[k - 1] + nb_elem_col[k - 1];
//     // FILLING pr and ir
//     // cur_ind contains the current shift in pr (and ir) for each column
//     for (k = 0; k < npts; ++k)
//         cur_ind[k] = jc[k];
//     for (k = 0; k < nsites; ++k) {
//         //eval_basis(nu[k], basis_u);
//         //eval_basis(nv[k], basis_v);
//         b_func_u(nu[k], basis_u);
//         b_func_v(nv[k], basis_v);
//         for (iu = 0; iu < 4; ++iu)
//             for (iv = 0; iv < 4; ++iv) {
//                 col = (iu + interu[k]) * bbs->nptsv + iv + interv[k];
//                 pr[cur_ind[col]] = fact * basis_u[iu] * basis_v[iv];
//                 ir[cur_ind[col]] = k;
//                 cur_ind[col] += 1;
//             }
//     }
// error:
//     free(nu);
//     free(nv);
//     free(interu);
//     free(interv);
//     free(nb_elem_col);
//     free(cur_ind);
//     return ret_code;
// }

/* Precomputed coefficients used in the bending matrix.
 * These arrays are used by the "bending_ur" function and they are
 * not intended to be used directly by the user. */
double __bxx_coeff[] = { 1.0 / 756.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0 / 504.0, 1.0 / 252.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0 / 504.0, 1.0 / 252.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 / 1512.0, 0.0, -1.0 / 504.0, 1.0 / 756.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 43.0 / 5040.0, -43.0 / 3360.0, 0.0, 43.0 / 10080.0, 11.0 / 140.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -43.0 / 3360.0, 43.0 / 1680.0, -43.0 / 3360.0, 0.0, -33.0 / 280.0, 33.0 / 140.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -43.0 / 3360.0, 43.0 / 1680.0, -43.0 / 3360.0, 0.0, -33.0 / 280.0, 33.0 / 140.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 43.0 / 10080.0, 0.0, -43.0 / 3360.0, 43.0 / 5040.0, 11.0 / 280.0, 0.0, -33.0 / 280.0, 11.0 / 140.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 / 252.0, -1.0 / 168.0, 0.0, 1.0 / 504.0, 311.0 / 5040.0, -311.0 / 3360.0, 0.0, 311.0 / 10080.0, 11.0 / 140.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0 / 168.0, 1.0 / 84.0, -1.0 / 168.0, 0.0, -311.0 / 3360.0, 311.0 / 1680.0, -311.0 / 3360.0, 0.0, -33.0 / 280.0, 33.0 / 140.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0 / 168.0, 1.0 / 84.0, -1.0 / 168.0, 0.0, -311.0 / 3360.0, 311.0 / 1680.0, -311.0 / 3360.0, 0.0, -33.0 / 280.0, 33.0 / 140.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 / 504.0, 0.0, -1.0 / 168.0, 1.0 / 252.0, 311.0 / 10080.0, 0.0, -311.0 / 3360.0, 311.0 / 5040.0, 11.0 / 280.0, 0.0, -33.0 / 280.0, 11.0 / 140.0, 0.0, 0.0, 0.0, 0.0, 1.0 / 15120.0, -1.0 / 10080.0, 0.0, 1.0 / 30240.0, 1.0 / 252.0, -1.0 / 168.0, 0.0, 1.0 / 504.0, 43.0 / 5040.0, -43.0 / 3360.0, 0.0, 43.0 / 10080.0, 1.0 / 756.0, 0.0, 0.0, 0.0, -1.0 / 10080.0, 1.0 / 5040.0, -1.0 / 10080.0, 0.0, -1.0 / 168.0, 1.0 / 84.0, -1.0 / 168.0, 0.0, -43.0 / 3360.0, 43.0 / 1680.0, -43.0 / 3360.0, 0.0, -1.0 / 504.0, 1.0 / 252.0, 0.0, 0.0, 0.0, -1.0 / 10080.0, 1.0 / 5040.0, -1.0 / 10080.0, 0.0, -1.0 / 168.0, 1.0 / 84.0, -1.0 / 168.0, 0.0, -43.0 / 3360.0, 43.0 / 1680.0, -43.0 / 3360.0, 0.0, -1.0 / 504.0, 1.0 / 252.0, 0.0, 1.0 / 30240.0, 0.0, -1.0 / 10080.0, 1.0 / 15120.0, 1.0 / 504.0, 0.0, -1.0 / 168.0, 1.0 / 252.0, 43.0 / 10080.0, 0.0, -43.0 / 3360.0, 43.0 / 5040.0, 1.0 / 1512.0, 0.0, -1.0 / 504.0, 1.0 / 756.0 };
double __bxy_coeff[] = { 1.0 / 200.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.0 / 1200.0, 17.0 / 600.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0 / 100.0, -29.0 / 1200.0, 17.0 / 600.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0 / 1200.0, -1.0 / 100.0, 7.0 / 1200.0, 1.0 / 200.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.0 / 1200.0, 49.0 / 7200.0, -7.0 / 600.0, -7.0 / 7200.0, 17.0 / 600.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 49.0 / 7200.0, 119.0 / 3600.0, -203.0 / 7200.0, -7.0 / 600.0, 119.0 / 3600.0, 289.0 / 1800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -7.0 / 600.0, -203.0 / 7200.0, 119.0 / 3600.0, 49.0 / 7200.0, -17.0 / 300.0, -493.0 / 3600.0, 289.0 / 1800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -7.0 / 7200.0, -7.0 / 600.0, 49.0 / 7200.0, 7.0 / 1200.0, -17.0 / 3600.0, -17.0 / 300.0, 119.0 / 3600.0, 17.0 / 600.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0 / 100.0, -7.0 / 600.0, 1.0 / 50.0, 1.0 / 600.0, -29.0 / 1200.0, -203.0 / 7200.0, 29.0 / 600.0, 29.0 / 7200.0, 17.0 / 600.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -7.0 / 600.0, -17.0 / 300.0, 29.0 / 600.0, 1.0 / 50.0, -203.0 / 7200.0, -493.0 / 3600.0, 841.0 / 7200.0, 29.0 / 600.0, 119.0 / 3600.0, 289.0 / 1800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 / 50.0, 29.0 / 600.0, -17.0 / 300.0, -7.0 / 600.0, 29.0 / 600.0, 841.0 / 7200.0, -493.0 / 3600.0, -203.0 / 7200.0, -17.0 / 300.0, -493.0 / 3600.0, 289.0 / 1800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 / 600.0, 1.0 / 50.0, -7.0 / 600.0, -1.0 / 100.0, 29.0 / 7200.0, 29.0 / 600.0, -203.0 / 7200.0, -29.0 / 1200.0, -17.0 / 3600.0, -17.0 / 300.0, 119.0 / 3600.0, 17.0 / 600.0, 0.0, 0.0, 0.0, 0.0, -1.0 / 1200.0, -7.0 / 7200.0, 1.0 / 600.0, 1.0 / 7200.0, -1.0 / 100.0, -7.0 / 600.0, 1.0 / 50.0, 1.0 / 600.0, 7.0 / 1200.0, 49.0 / 7200.0, -7.0 / 600.0, -7.0 / 7200.0, 1.0 / 200.0, 0.0, 0.0, 0.0, -7.0 / 7200.0, -17.0 / 3600.0, 29.0 / 7200.0, 1.0 / 600.0, -7.0 / 600.0, -17.0 / 300.0, 29.0 / 600.0, 1.0 / 50.0, 49.0 / 7200.0, 119.0 / 3600.0, -203.0 / 7200.0, -7.0 / 600.0, 7.0 / 1200.0, 17.0 / 600.0, 0.0, 0.0, 1.0 / 600.0, 29.0 / 7200.0, -17.0 / 3600.0, -7.0 / 7200.0, 1.0 / 50.0, 29.0 / 600.0, -17.0 / 300.0, -7.0 / 600.0, -7.0 / 600.0, -203.0 / 7200.0, 119.0 / 3600.0, 49.0 / 7200.0, -1.0 / 100.0, -29.0 / 1200.0, 17.0 / 600.0, 0.0, 1.0 / 7200.0, 1.0 / 600.0, -7.0 / 7200.0, -1.0 / 1200.0, 1.0 / 600.0, 1.0 / 50.0, -7.0 / 600.0, -1.0 / 100.0, -7.0 / 7200.0, -7.0 / 600.0, 49.0 / 7200.0, 7.0 / 1200.0, -1.0 / 1200.0, -1.0 / 100.0, 7.0 / 1200.0, 1.0 / 200.0 };
double __byy_coeff[] = { 1.0 / 756.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 43.0 / 5040.0, 11.0 / 140.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 / 252.0, 311.0 / 5040.0, 11.0 / 140.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 / 15120.0, 1.0 / 252.0, 43.0 / 5040.0, 1.0 / 756.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0 / 504.0, -43.0 / 3360.0, -1.0 / 168.0, -1.0 / 10080.0, 1.0 / 252.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -43.0 / 3360.0, -33.0 / 280.0, -311.0 / 3360.0, -1.0 / 168.0, 43.0 / 1680.0, 33.0 / 140.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0 / 168.0, -311.0 / 3360.0, -33.0 / 280.0, -43.0 / 3360.0, 1.0 / 84.0, 311.0 / 1680.0, 33.0 / 140.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0 / 10080.0, -1.0 / 168.0, -43.0 / 3360.0, -1.0 / 504.0, 1.0 / 5040.0, 1.0 / 84.0, 43.0 / 1680.0, 1.0 / 252.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0 / 504.0, -43.0 / 3360.0, -1.0 / 168.0, -1.0 / 10080.0, 1.0 / 252.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -43.0 / 3360.0, -33.0 / 280.0, -311.0 / 3360.0, -1.0 / 168.0, 43.0 / 1680.0, 33.0 / 140.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0 / 168.0, -311.0 / 3360.0, -33.0 / 280.0, -43.0 / 3360.0, 1.0 / 84.0, 311.0 / 1680.0, 33.0 / 140.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0 / 10080.0, -1.0 / 168.0, -43.0 / 3360.0, -1.0 / 504.0, 1.0 / 5040.0, 1.0 / 84.0, 43.0 / 1680.0, 1.0 / 252.0, 0.0, 0.0, 0.0, 0.0, 1.0 / 1512.0, 43.0 / 10080.0, 1.0 / 504.0, 1.0 / 30240.0, 0.0, 0.0, 0.0, 0.0, -1.0 / 504.0, -43.0 / 3360.0, -1.0 / 168.0, -1.0 / 10080.0, 1.0 / 756.0, 0.0, 0.0, 0.0, 43.0 / 10080.0, 11.0 / 280.0, 311.0 / 10080.0, 1.0 / 504.0, 0.0, 0.0, 0.0, 0.0, -43.0 / 3360.0, -33.0 / 280.0, -311.0 / 3360.0, -1.0 / 168.0, 43.0 / 5040.0, 11.0 / 140.0, 0.0, 0.0, 1.0 / 504.0, 311.0 / 10080.0, 11.0 / 280.0, 43.0 / 10080.0, 0.0, 0.0, 0.0, 0.0, -1.0 / 168.0, -311.0 / 3360.0, -33.0 / 280.0, -43.0 / 3360.0, 1.0 / 252.0, 311.0 / 5040.0, 11.0 / 140.0, 0.0, 1.0 / 30240.0, 1.0 / 504.0, 43.0 / 10080.0, 1.0 / 1512.0, 0.0, 0.0, 0.0, 0.0, -1.0 / 10080.0, -1.0 / 168.0, -43.0 / 3360.0, -1.0 / 504.0, 1.0 / 15120.0, 1.0 / 252.0, 43.0 / 5040.0, 1.0 / 756.0 };

/* This function computes the upper right part of the "bending matrix".
 * - bbs: the bidimensional b-spline structure.
 * - lambdas: array storing an (nptsv-3)x(nptsu-3) matrix by column
 *   lambdas[iu*(nptsv-3)+iv] is the regularization parameter over the knot interval #(iu,iv).
 * - pr, ir, jc: internal storage for the resulting sparse bending matrix.
 */

void bending_ur_new(bbs_t* bbs, double* lambdas, double* pr, MatlabSparse::StorageIndex* ir, MatlabSparse::StorageIndex* jc) {
    int ny = bbs->nptsu; // Order of the dimensions are switched due to a previous version of the code
    int nx = bbs->nptsv;
    int i, j, a, b, c, d, e1, f1, e2, f2, curnb, total, ind, sb, nb;
    MatlabSparse::StorageIndex* pi = ir, * pp = jc;
    double* px = pr;
    double coeff_b[256];
    double lbd;
    double sy = (bbs->umax - bbs->umin) / (bbs->nptsu - 3);
    double sx = (bbs->vmax - bbs->vmin) / (bbs->nptsv - 3);

    *pp = 0;
    ++pp;
    total = 0;

    /* Location of the non-zeros coefficients of the matrix. */
    for (j = 0; j < ny; ++j) {
        for (i = 0; i < nx; ++i) {
            curnb = 0;
            for (b = maxx(j - 3, 0); b <= j - 1; ++b)
                for (a = maxx(0, i - 3); a <= min(nx - 1, i + 3); ++a) {
                    *pi = b * nx + a;
                    *px = 0.0;
                    ++pi;
                    ++px;
                    ++curnb;
                }

            for (a = maxx(0, i - 3); a <= i; ++a) {
                *pi = j * nx + a;
                *px = 0.0;
                ++pi;
                ++px;
                ++curnb;
            }
            total += curnb;
            *pp = total;
            ++pp;
        }
    }

    // Build the coefficients of the B matrix with the scales taken into account
    for (j = 0; j < 16; ++j)
        for (i = 0; i <= j; ++i) {
            ind = 16 * j + i;
            coeff_b[ind] = sy * __bxx_coeff[ind] / pow(sx, 3)
                + __bxy_coeff[ind] / (sx * sy)
                + sx * __byy_coeff[ind] / pow(sy, 3);
        }


    // Put the right coefficients at the right locations (one knot domain at a time
    // and with the scaling given by lambda)
    pp = jc;
    px = pr;
    for (b = 0; b < ny - 3; ++b) {
        for (a = 0; a < nx - 3; ++a) {
            lbd = lambdas[b * (nx - 3) + a];
            for (c = 0; c < 16; ++c) {
                for (d = c; d < 16; ++d) {
                    e1 = c / 4;
                    f1 = c % 4;
                    e2 = d / 4;
                    f2 = d % 4;
                    i = (b + e1) * nx + a + f1;
                    j = (b + e2) * nx + a + f2;

                    // See (notebook 2, p. 61) for the tricky formulas
                    nb = i / nx - maxx(j / nx - 3, 0);
                    sb = min(min(4 + (j % nx), 3 + nx - (j % nx)), min(nx, 7));
                    px[pp[j] + nb * sb + (i % nx) - maxx((j % nx) - 3, 0)] += lbd * coeff_b[16 * d + c];
                }
            }
        }
    }
}

void my_coloc_new(int nlhs, MatlabSparse &plhs, int nrhs, bbs_t bbs, VectorXd prhs2, VectorXd prhs3) {
    //bbs_t bbs;
    double* u = NULL, * v = NULL;
    int nb_u, nb_v;
    double* pr;
    MatlabSparse::StorageIndex* ir, * jc;

    // INPUT ARGUMENTS
    if (nrhs != 3)
        cout << "Three inputs required." << endl;

    // No need to check "prhs[0]" for validity since it is done in "array_to_bbs"
    //array_to_bbs(prhs1[0], &bbs);

    // Check u and v
    /*check_real(prhs2, (char*)"u");
    check_real(prhs3, (char*)"v");
    nb_u = mxGetNumberOfElements(prhs2);
    nb_v = mxGetNumberOfElements(prhs3);
    if (nb_u != nb_v)
        mexErrMsgTxt("'u' and 'v' must have the same number of elements.");*/
    // nb_u = mxGetNumberOfElements(prhs2);
    // u = mxGetPr(prhs2); // real values in prhs2
    // v = mxGetPr(prhs3); // real values in prhs3

    nb_u = prhs2.size();
    MatlabSparse prhs2_sparse = prhs2.sparseView();
    MatlabSparse prhs3_sparse = prhs3.sparseView();
    u = prhs2_sparse.valuePtr(); // real values in prhs2
    v = prhs3_sparse.valuePtr(); // real values in prhs3

    //for (int i = 0; i < 889; i++) {
    //	cout << u[i] << endl;
    //}

    // OUTPUT ARGUMENTS
    if (nlhs != 1)
        cout << "One output required." << endl;
    
    // plhs[0] = mxCreateSparse(nb_u, bbs.nptsu * bbs.nptsv, 16 * nb_u, mxREAL);
    // pr = mxGetPr(plhs[0]);
    // ir = mxGetIr(plhs[0]);
    // jc = mxGetJc(plhs[0]);

    plhs.resize(nb_u,bbs.nptsu * bbs.nptsv);
    plhs.reserve(16 * nb_u);
    pr = plhs.valuePtr();
    ir = plhs.innerIndexPtr();
    jc = plhs.outerIndexPtr();

    // ACTUAL EVALUATION
    switch (coloc_new(&bbs, u, v, nb_u, pr, ir, jc)) {
    case 1: cout << "A colocation site was outside of the spline definition domain." << endl;
    }
}

MatlabSparse my_bbs_coloc_new(bbs_t bbs, DoubleVector1D p_0, DoubleVector1D p_1)
{
    int n = p_0.size();
    int nlhs(1);
    int nrhs(3);

    // mxArray* prhs2 = mxCreateDoubleMatrix(1, n, mxREAL);
    // memcpy(mxGetData(prhs2), &p_0[0], n * sizeof(double));
    // mxArray* prhs3 = mxCreateDoubleMatrix(1, n, mxREAL);
    // memcpy(mxGetData(prhs3), &p_1[0], n * sizeof(double));

    VectorXd prhs2_new = vector_to_eigen_1D(p_0);
    VectorXd prhs3_new = vector_to_eigen_1D(p_1);

    // mxArray* plhs;
    MatlabSparse plhs_new;
    // my_coloc(nlhs, &plhs, nrhs, bbs, prhs2, prhs3);
    my_coloc_new(nlhs, plhs_new, nrhs, bbs, prhs2_new, prhs3_new);

    return plhs_new;
}

void my_bending_new(int nlhs, MatlabSparse &plhs, bbs_t bbs, MatrixXd prhs) {
    double* pr;
    MatlabSparse::StorageIndex* ir, * jc;
    double* lambdas;

    // Check u and v
    int nu = bbs.nptsu;
    int nv = bbs.nptsv;

    // lambdas = mxGetPr(prhs);

    MatlabSparse prhs_sparse = prhs.sparseView();
    lambdas = prhs_sparse.valuePtr(); // real values in prhs
    

    // OUTPUT ARGUMENTS
    if (nlhs != 1)
        cout << "One output required." << endl;
    // plhs[0] = mxCreateSparse(nu * nv, nu * nv, 25 * nu * nv - 42 * (nu + nv) + 72, mxREAL);
    // pr = mxGetPr(plhs[0]);
    // ir = mxGetIr(plhs[0]);
    // jc = mxGetJc(plhs[0]);

    plhs.resize(nu * nv,nu * nv);
    plhs.reserve(25 * nu * nv - 42 * (nu + nv) + 72);
    pr = plhs.valuePtr();
    ir = plhs.innerIndexPtr();
    jc = plhs.outerIndexPtr();

    //cout << 25 * nu * nv - 42 * (nu + nv) + 72 << endl;

    // ACTUAL COMPUTATION OF THE BENDING MATRIX
    bending_ur_new(&bbs, lambdas, pr, ir, jc);
}

MatlabSparse my_bbs_bending_new(bbs_t bbs, DoubleVector2D lambdas)
{
    int n = lambdas.size();
    int nlhs(1);
    int nrhs(2);

    // DoubleVector1D lambdas_linear(n * n);
    // int s(0);
    // for (int i = 0; i < n; i++)
    // {
    //     for (int j = 0; j < n; j++)
    //     {
    //         lambdas_linear[s] = lambdas[i][j];
    //         s++;
    //     }
    // }

    // mxArray* prhs2 = mxCreateDoubleMatrix(n, n, mxREAL);
    // memcpy(mxGetData(prhs2), &lambdas_linear[0], n * n * sizeof(double));

    MatrixXd prhs2 = vector_to_eigen(lambdas);

    /*int row(mxGetM(prhs2));
    int col(mxGetN(prhs2));
    for (int i = 0; i < row+1; i++)
    {
        for (int j = 0; j < col; j++)
        {
            double* in1pr = mxGetPr(prhs2);
            mexPrintf("row = %d col = %d in1=%f\n", i, j, in1pr[i*col + j]);
        }
    }*/

    /*for (int i = 0; i < n * n; i++) {
        double* in1pr = mxGetPr(prhs2);
        mexPrintf("in1=%f\n", in1pr[i]);
    }*/

    // mxArray* plhs;
    MatlabSparse plhs_new;
    my_bending_new(nlhs, plhs_new, bbs, prhs2);

    return plhs_new;
}

void my_eval_new(int nlhs, MatlabSparse &plhs, bbs_t bbs, MatrixXd prhs1, VectorXd prhs2, VectorXd prhs3, int du, int dv) {
    //bbs_t bbs;
    double* u = NULL, * v = NULL;
    int nb_u, nb_v;
    double* val = NULL;
    double* ctrlpts = NULL;

    // No need to check "prhs[0]" for validity since it is done in "array_to_bbs"
    //array_to_bbs(prhs1[0], &bbs);
    // check_real(prhs1, (char *)"ctrlpts");
    // check_bbs_ctrlpts_size(prhs1, &bbs);
    // ctrlpts = mxGetPr(prhs1);
    MatlabSparse prhs1_sparse = prhs1.sparseView();
    ctrlpts = prhs1_sparse.valuePtr();

    // Check u and v
    /*check_real(prhs2, (char*)"u");
    check_real(prhs3, (char*)"v");
    nb_u = mxGetNumberOfElements(prhs2);
    nb_v = mxGetNumberOfElements(prhs3);
    if (nb_u != nb_v)
        mexErrMsgTxt("'u' and 'v' must have the same number of elements.");*/
    // nb_u = mxGetNumberOfElements(prhs2);
    // u = mxGetPr(prhs2);
    // v = mxGetPr(prhs3);

    nb_u = prhs2.size();
    MatlabSparse prhs2_sparse = prhs2.sparseView();
    MatlabSparse prhs3_sparse = prhs3.sparseView();
    u = prhs2_sparse.valuePtr();
    v = prhs3_sparse.valuePtr();

    //for (int i = 0; i < 889; i++) {
    //	cout << u[i] << endl;
    //}


    // OUTPUT ARGUMENTS
    if (nlhs != 1)
        cout << "One output required." << endl;
    // plhs[0] = mxCreateDoubleMatrix(bbs.valdim, nb_u, mxREAL);
    // val = mxGetPr(plhs[0]);

    plhs.resize(bbs.valdim, nb_u);
    plhs.reserve(bbs.valdim * nb_u); 
    val = plhs.valuePtr();

    // cout << val << endl;

    // for (int i = 0; i < nb_u; i++)
    // {
    //     cout << u[i] << " " << v[i] << endl;
    // }

    // ACTUAL EVALUATION
    // eval(&bbs, ctrlpts, u, v, nb_u, val, du, dv);
    eval_new(&bbs, ctrlpts, u, v, nb_u, val, du, dv);
}

MatlabSparse my_bbs_eval_new(bbs_t bbs, DoubleVector2D ctrlpts, DoubleVector1D p_0, DoubleVector1D p_1, int du, int dv)
{
    int n = p_0.size();
    int n2 = ctrlpts[0].size();
    int nlhs(1);
    int nrhs(6);

    // mxArray* prhs2 = mxCreateDoubleMatrix(1, n, mxREAL);
    // memcpy(mxGetData(prhs2), &p_0[0], n * sizeof(double));
    // mxArray* prhs3 = mxCreateDoubleMatrix(1, n, mxREAL);
    // memcpy(mxGetData(prhs3), &p_1[0], n * sizeof(double));

    VectorXd prhs2 = vector_to_eigen_1D(p_0);
    VectorXd prhs3 = vector_to_eigen_1D(p_1);

    // DoubleVector1D ctrlpts_linear(2 * n2);
    // int s(0);
    // for (int j = 0; j < n2; j++)
    // {
    //     for (int i = 0; i < 2; i++)
    //     {
    //         ctrlpts_linear[s] = ctrlpts[i][j];
    //         s++;
    //     }
    // }

    // mxArray* prhs1 = mxCreateDoubleMatrix(2, n2, mxREAL);
    // memcpy(mxGetData(prhs1), &ctrlpts_linear[0], 2 * n2 * sizeof(double));

    MatrixXd prhs1 = vector_to_eigen(ctrlpts);

    /*int row(mxGetM(prhs4));
    int col(mxGetN(prhs4));
    cout << row << " " << col << endl;
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            double* in1pr = mxGetPr(prhs4);
            mexPrintf("i=%d j=%d in1=%f\n", i,j,in1pr[i* col + j]);
        }
    }*/

    // mxArray* plhs;
    MatlabSparse plhs;
    // my_eval(nlhs, plhs, bbs, prhs1, prhs2, prhs3, du, dv);
    my_eval_new(nlhs, plhs, bbs, prhs1, prhs2, prhs3, du, dv);

    return plhs;
}

DoubleVector2D BBS_Function(bbs_t bbs, DoubleVector2D pixles_template, DoubleVector2D pixels_image, DoubleVector2D K, int nC, double er)
{
    
    //------------------------------------------------
    DoubleVector1D p_0(pixles_template[0]);
    DoubleVector1D p_1(pixles_template[1]);
    //------------------------------------------------
    int n_visible_particles{ static_cast<int>(pixels_image[0].size()) };
    double fx{ K[0][0] };
    double fy{ K[1][1] };
    double u0{ K[0][2] };
    double v0{ K[1][2] };
    DoubleVector2D xy(2, vector<double>(n_visible_particles));
    for (int i = 0; i < n_visible_particles; i++)
    {
        xy[0][i] = (pixels_image[0][i] - u0) / fx;
        xy[1][i] = (pixels_image[1][i] - v0) / fy;
    }

    DoubleVector2D q(xy);
    DoubleVector2D q_temp(transpose(q));
    //------------------------------------------------
    // mxArray* coloc(my_bbs_coloc(bbs, p_0, p_1));
    MatlabSparse coloc(my_bbs_coloc_new(bbs, p_0, p_1));
    //------------------------------------------------

    DoubleVector2D lambdas(nC - 3, vector<double>(nC - 3));
    for (int i = 0; i < lambdas.size(); i++)
    {
        for (int j = 0; j < lambdas[0].size(); j++)
        {
            lambdas[i][j] = er;
        }
    }
    //------------------------------------------------
    // mxArray* bending(my_bbs_bending(bbs, lambdas));
    MatlabSparse bending(my_bbs_bending_new(bbs, lambdas));
    //------------------------------------------------
    //auto start_warp1 = chrono::high_resolution_clock::now();

    // MatlabSparse eigen_coloc(coloc);
    MatlabSparse eigen_coloc_transpose = coloc.transpose();
    // MatlabSparse eigen_bending_temp(bending);
    MatlabSparse eigen_bending = bending.selfadjointView<Upper>();

    //eigen_bending.triangularView<Lower>() = eigen_bending_transpose.triangularView<Upper>();

    MatrixXd mat(q_temp.size(), q_temp[0].size());
    for (int i = 0; i < q_temp.size(); i++)
        mat.row(i) = VectorXd::Map(&q_temp[i][0], q_temp[i].size());

    SparseMatrix<double> mat_sparse = mat.sparseView();

    MatlabSparse A(eigen_coloc_transpose * coloc + eigen_bending);
    MatlabSparse B(eigen_coloc_transpose * mat_sparse);

    // BiCGSTAB<Eigen::SparseMatrix<double>> chol(A);
    //MatrixXd x = solver.solve(B);
    SparseLU<Eigen::SparseMatrix<double>> chol(A);
    /*if (chol.info() != Eigen::Success)
        cout << "bad result!" << endl;*/

    MatlabSparse eigen_cpts = chol.solve(B);
    MatlabSparse eigen_ctrlpts(eigen_cpts.transpose());   

    //mxArray* ctrlpts(eigen_to_matlab_sparse(eigen_ctrlpts));
    //cout << "Sum error = " << (A * eigen_cpts - B).sum() << endl;

    //------------------------------------------------
    DoubleVector2D ctrlpts(2, vector<double>(nC * nC));
    for (int i = 0; i < ctrlpts.size(); i++)
    {
        for (int j = 0; j < ctrlpts[0].size(); j++)
        {
            ctrlpts[i][j] = eigen_ctrlpts.coeff(i, j);
        }
    }
    //------------------------------------------------

    return ctrlpts;
}

DoubleVector2D BBS_Evaluation(bbs_t bbs, DoubleVector2D ctrlpts, DoubleVector2D pixles, DoubleVector2D K, int deriv_1, int deriv_2)
{
    //------------------------------------------------
    DoubleVector1D pixles_0(pixles[0]);
    DoubleVector1D pixles_1(pixles[1]);

    for (int i = 0; i < pixles_0.size(); i++)
    {
        if(pixles_0[i] == 0) pixles_0[i] += 0.0001;
        if(pixles_1[i] == 0) pixles_1[i] += 0.0001;
    }

    //------------------------------------------------
    // mxArray* qw1(my_bbs_eval(bbs, ctrlpts, pixles_0, pixles_1, 0, 0));
    // int row1(mxGetM(qw1));
    // int col1(mxGetN(qw1));
    // DoubleVector2D qw1_myFormat(3, vector<double>(col1));
    // double* in1pr1 = mxGetPr(qw1);

    
    MatlabSparse qw(my_bbs_eval_new(bbs, ctrlpts, pixles_0, pixles_1, 0, 0));
    

    int row(qw.rows());
    int col(qw.cols());
    DoubleVector2D qw_myFormat(3, vector<double>(col));
    double* in1pr = qw.valuePtr();
    for (int j = 0; j < col; j++)
    {
        for (int i = 0; i < row; i++)
        {

            //mexPrintf("row = %d col = %d in1=%f\n", i, j, in1pr[i*col + j]);
            // cout << in1pr[j * row + i] << endl;
            qw_myFormat[i][j] = in1pr[j * row + i];
        }
        qw_myFormat[2][j] = 1;
    }

    //------------------------------------------------
    DoubleVector2D pixels_grid(projection(K, qw_myFormat));
    //------------------------------------------------
    return pixels_grid;
}

