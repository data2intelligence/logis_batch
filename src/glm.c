#include "glm.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multiroots.h>

// numerical tolerance for really small value
#define EPS 1e-10

// EPS correction for sqrt function on nearly 0 value
#define SQRT_EPS(x) (fabs(x)<EPS?0:sqrt(x))

// char values of two category response (for logistic regression, etc)
char two_category_response[2] = {0,1};


// workspace for logistic regression
typedef struct logistic_regression_workspace
{
	// global input parameters, no need to be allocated
	const gsl_matrix *X;
	const gsl_vector *Y;

	// dimension of X : n*p. length of Y: n
	size_t n, p;

	// internal variables for iterative optimization
	// For Newton Raphson
	gsl_vector *P, *W, *error, *H, *U, *delta;
	gsl_matrix *wX, *H_t, *I, *L;

} LGR_space;


// allocate the internal variable space
LGR_space * alloc_LGR_space(const gsl_matrix *X, const gsl_vector *Y)
{
	size_t n = X->size1, p = X->size2;
	LGR_space *lr = (LGR_space*)malloc(sizeof(LGR_space));

	lr->X = X;
	lr->Y = Y;

	lr->n = n;
	lr->p = p;

	lr->P = gsl_vector_alloc(n);		// P(Y=1|X)
	lr->W = gsl_vector_alloc(n);		// diagonal elements : P(1-P)
	lr->error = gsl_vector_alloc(n);	// Y-P
	lr->H = gsl_vector_alloc(n);		// diagonal elements of Hat matrix: (W^0.5)X I^-1 X'(W^0.5)
	lr->U = gsl_vector_alloc(p);		// gradient: t(Y-P)X
	lr->delta = gsl_vector_alloc(p);	// change step of beta

	lr->wX = gsl_matrix_alloc(n, p);	// (W^0.5)X
	lr->H_t = gsl_matrix_alloc(n, p);// auxiliary space for H diagonal calculation
	lr->I = gsl_matrix_alloc(p, p);	// Information matrix: X'WX
	lr->L = gsl_matrix_alloc(p, p);	// Cholesky decomposition L

	// handle the singular matrix in Cholesky decomposition by non-exit procedure
	gsl_set_error_handler_off ();

	return lr;
}


void free_LGR_space(LGR_space *ls)
{
	gsl_vector_free(ls->P);
	gsl_vector_free(ls->W);
	gsl_vector_free(ls->error);
	gsl_vector_free(ls->U);
	gsl_vector_free(ls->delta);
	gsl_vector_free(ls->H);

	gsl_matrix_free(ls->wX);
	gsl_matrix_free(ls->H_t);
	gsl_matrix_free(ls->I);
	gsl_matrix_free(ls->L);

	free(ls);
}


// swap the content of pointer for two vectors
inline void swap_gsl_vector(gsl_vector **a, gsl_vector **b)
{
	gsl_vector *c = *a;
	*a = *b;
	*b = c;
}


// logistic regression core procedure: return number of iterations if succeeded.
int logistic_regression_core(
		const size_t maxIter, const double tol, const double max_delta,
		const int correction,
		gsl_vector *beta, gsl_vector *sderr, gsl_vector *z, gsl_vector *pvalue,
		LGR_space *ls)
{
	double w, norm_delta, norm_beta, step_ratio;

	int status=GSL_CONTINUE, returnflag = 0;

	size_t i,j, n=ls->n, p=ls->p;

	// hook up to workspace
	const gsl_matrix *X = ls->X;
	const gsl_vector *Y = ls->Y;

	gsl_vector *P = ls->P, *W = ls->W, *error = ls->error, *H = ls->H,
			*U = ls->U, *delta = ls->delta;

	gsl_matrix *wX = ls->wX, *H_t = ls->H_t, *I = ls->I, *L = ls->L;

	for(i=0 ; i<maxIter ; i++)
	{
		// probability vector P = sigmod(X*beta)
		gsl_blas_dgemv(CblasNoTrans, 1, X,beta, 0, P);

		for(j=0;j<P->size;j++)
		{
			w = gsl_vector_get(P,j);
			gsl_vector_set(P, j, gsl_cdf_logistic_P(w,1));
		}

		// diagonal(W) = P*(1-P)
		gsl_vector_memcpy(W, P);
		gsl_vector_mul(W,P);
		gsl_vector_sub(W,P);
		gsl_vector_scale(W,-1);

		// error = Y-P
		gsl_vector_memcpy(error, Y);
		gsl_vector_sub(error,P);


		// Information matrix I = X'WX , Hessian matrix = -I

		// wX = W^0.5 * X
		gsl_matrix_memcpy(wX,X);

		for(j=0;j<n;j++)
		{
			gsl_vector_view X_row = gsl_matrix_row(wX, j);
			w = gsl_vector_get(W, j);
			gsl_vector_scale(&X_row.vector, SQRT_EPS(w));
		}

		// I = (W^0.5 *X)' (W^0.5 *X) = X'WX
		gsl_blas_dsyrk(CblasLower, CblasTrans, 1, wX, 0, I);

		// diagonal and lower-triangular part of the matrix are used
		if(gsl_linalg_cholesky_decomp(I) == GSL_EDOM)
		{
			//fprintf(stderr, "Cholesky decomposition failed on information matrix X'WX\n");

			// clear up all results
			gsl_vector_set_zero(beta);
			gsl_vector_set_zero(sderr);
			gsl_vector_set_zero(z);
			gsl_vector_set_all(pvalue, 1);

			return REG_FAIL;
		}

		// we will change I to I^-1, so copy the Cholesky L
		gsl_matrix_memcpy(L, I);

		// compute the inverse of information matrix, now I is I^-1
		gsl_linalg_cholesky_invert(I);

		// Correct on delta, but still use uncorrected I^-1 as first order approximation
		if(correction)
		{
			gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, wX, I, 0, H_t);

			// get diagonal elements Hii
			for(j=0;j<n;j++)
			{
				gsl_vector_view H_t_r = gsl_matrix_row(H_t, j);
				gsl_vector_view wX_r = gsl_matrix_row(wX, j);

				gsl_blas_ddot(&H_t_r.vector, &wX_r.vector, &w);
				gsl_vector_set(H, j, w);
			}

			// P - 0.5. P is changed from now
			gsl_vector_add_constant (P, -0.5);

			// H*(P - 0.5)
			gsl_vector_mul(H, P);

			// correct the error vector
			gsl_vector_sub(error, H);
		}

		//gradient U = X'error
		gsl_blas_dgemv(CblasTrans, 1, X, error, 0, U);

		// stop criterion by |U| almost 0
		status = gsl_multiroot_test_residual(U, tol);
		if(status == GSL_SUCCESS) break;


		// next round
		// Newton Raphson delta: I^-1 * U
		gsl_linalg_cholesky_solve(L, U, delta);

		norm_delta = gsl_blas_dnrm2(delta);
		norm_beta = gsl_blas_dnrm2(beta);
		step_ratio = (norm_delta + 1)/(norm_beta + 1);

		if (step_ratio > max_delta)
		{
			step_ratio = max_delta/step_ratio;
			gsl_vector_scale(delta, step_ratio);
		}

		// beta = beta + delta
		gsl_vector_add(beta, delta);
	}

	if(i==maxIter && status != GSL_SUCCESS)
	{
		//fprintf(stderr, "Exceed maximum number of iterations %lu norm U = %lf\n", maxIter, gsl_blas_dnrm2(U));
		returnflag = CONVERGE_FAIL;
	}else{
		returnflag = (int)i;	// return number of iterations
	}

	// calculate statistic measures for inference

	// take diagonal elements of information matrix as variance
	gsl_vector_view Irev_D = gsl_matrix_diagonal(I);
	gsl_vector_memcpy(sderr, &Irev_D.vector);

	for(i=0 ; i<p ; i++)
	{
		w = SQRT_EPS(gsl_vector_get(sderr, i));
		gsl_vector_set(sderr,i,w);

		if(w>0){
			w = gsl_vector_get(beta, i)/w;
			gsl_vector_set(z,i,w);

			// Wald test: Pr(X>|z|)
			w = 2*(1-gsl_cdf_gaussian_P(fabs(w),1));
			gsl_vector_set(pvalue,i,w);

		}else{
			fprintf(stderr, "standard error of beta is zero\n");

			gsl_vector_set(z,i,0);
			gsl_vector_set(pvalue,i,1);
		}
	}

	return returnflag;
}


size_t logistic_regression_batch(const gsl_matrix *X, const gsl_matrix_char *Y,
		const size_t maxIter, const double tol, const double max_delta, const size_t cntthres,
		const int correction, const int verbose,
		gsl_matrix *beta, gsl_matrix *sderr, gsl_matrix *z, gsl_matrix *pvalue)
{
	int status, guess;

	size_t i, j, total, n=X->size1, p=X->size2, k = Y->size2, success = 0, failed = 0, notconverge = 0;
	gsl_vector_view betai, sderri, zi, pvaluei;

	// double Y vector for each sample of Y matrix
	gsl_vector *Yi = gsl_vector_alloc(n);

	// make sure matrix Y only contains 0 and 1
	check_matrix_category(Y, two_category_response, 2, "Y");

	// matrix should stay in row major order for result vectors
	check_matrix_dimension((gsl_matrix*)Y, n, k, "Y", 1);

	check_matrix_dimension(beta, k, p, "beta", 1);
	check_matrix_dimension(sderr, k, p, "sderr", 1);
	check_matrix_dimension(z, k, p, "z", 1);
	check_matrix_dimension(pvalue, k, p, "pvalue", 1);

	LGR_space *ls = alloc_LGR_space(X,Yi);

	for(i=0;i<k;i++)
	{
		gsl_vector_char_const_view Y_row = gsl_matrix_char_const_column(Y,i);
		total = vector_char_to_double(Yi, &Y_row.vector);

		betai = gsl_matrix_row(beta,i);
		sderri = gsl_matrix_row(sderr,i);
		zi = gsl_matrix_row(z,i);
		pvaluei = gsl_matrix_row(pvalue,i);

		// jump over rare samples
		if(total < cntthres)
		{
			gsl_vector_set_zero(&betai.vector);
			gsl_vector_set_zero(&sderri.vector);
			gsl_vector_set_zero(&zi.vector);
			gsl_vector_set_all(&pvaluei.vector, 1);

			continue;
		}

		// if beta start from non-zero, we assume the user give it a guess
		for(guess=0, j=0;j<betai.vector.size;j++)
		{
			if(gsl_vector_get(&betai.vector, j) != 0)
			{
				guess = 1;
				break;
			}
		}

		status = logistic_regression_core(maxIter, tol, max_delta, correction,
				&betai.vector, &sderri.vector, &zi.vector, &pvaluei.vector, ls);

		// failed by and start has a guess, try from zero again
		if(status < 0 && guess != 0)
		{
			gsl_vector_set_zero(&betai.vector);
			status = logistic_regression_core(maxIter, tol, max_delta, correction,
						&betai.vector, &sderri.vector, &zi.vector, &pvaluei.vector, ls);
		}

		if(status == REG_FAIL){
			if(verbose) fprintf(stderr, "Regression failed on logistic regression %lu\n", i);
			failed++;
		}else if(status == CONVERGE_FAIL){
			if(verbose) fprintf(stderr, "Convergence failed on logistic regression %lu\n", i);
			notconverge++;
		}else{
			success++;
		}
	}

	free_LGR_space(ls);
	gsl_vector_free(Yi);

	if(verbose) fprintf(stdout, "success: %lu, failed: %lu, not converge: %lu\n", success, failed, notconverge);

	return success;
}


int logistic_regression(const gsl_matrix *X, const gsl_vector_char *Y,
		const size_t maxIter, const double tol, const double max_delta,
		const int correction, const int verbose,
		gsl_vector *beta, gsl_vector *sderr, gsl_vector *z, gsl_vector *pvalue)
{
	int status;

	size_t n=X->size1, p=X->size2;
	gsl_vector *Yi = gsl_vector_alloc(n);

	// make sure vector Y only contains 0 and 1
	check_vector_category(Y, two_category_response, 2, "Y");

	// check vector length
	check_vector_dimension((gsl_vector*)Y, n, "Y", 1);
	check_vector_dimension(beta, p, "beta", 1);
	check_vector_dimension(sderr, p, "sderr", 1);
	check_vector_dimension(z, p, "z", 1);
	check_vector_dimension(pvalue, p, "pvalue", 1);

	vector_char_to_double(Yi, Y);
	LGR_space *ls = alloc_LGR_space(X,Yi);

	status = logistic_regression_core(maxIter, tol, max_delta, correction, beta, sderr, z, pvalue, ls);

	free_LGR_space(ls);
	gsl_vector_free(Yi);

	return status;
}


// X:n*p, Y: n*k, beta:k*p
int OLS(const gsl_matrix *X, const gsl_matrix *Y, gsl_matrix *beta)
{
	size_t n=X->size1, p=X->size2, k=Y->size2;

	gsl_matrix *I = gsl_matrix_alloc(p, p), *T = gsl_matrix_alloc(p, k);

	check_matrix_dimension(Y, n, k, "Y", 1);
	check_matrix_dimension(beta, k, p, "beta", 1);

	// I = X'X
	gsl_blas_dsyrk(CblasLower, CblasTrans, 1, X, 0, I);

	if(gsl_linalg_cholesky_decomp(I) == GSL_EDOM)
	{
		fprintf(stderr, "Cholesky decomposition failed on X'X.\n");
		exit(1);
	}

	// I = (X'X)^-1
	gsl_linalg_cholesky_invert(I);

	// T = X'Y
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, X, Y, 0, T);

	// beta = (X'X)^-1 X'Y
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, T, I, 0, beta);

	gsl_matrix_free(I);
	gsl_matrix_free(T);

	return 0;
}


void print_vector(const gsl_vector *v, const char *title, FILE *fp)
{
	size_t i;

	fprintf(fp, "Vector %s:", title);

	for(i=0;i<v->size;i++)
	{
		fprintf(fp,"\t%e",gsl_vector_get(v,i));
	}
	fprintf(fp,"\n");
}

void print_matrix(const gsl_matrix *X, const char *title, FILE *fp)
{
	size_t i,j;

	fprintf(fp,"Matrix %s:\n", title);

	for(i=0;i<X->size1;i++)
	{
		for(j=0;j<X->size2;j++)
		{
			fprintf(fp,"%e\t",gsl_matrix_get(X,i,j));
		}

		fprintf(fp,"\n");
	}
}

size_t vector_char_to_double(gsl_vector *dest, const gsl_vector_char *src)
{
	char c;
	size_t i, n = dest->size, total = 0;

	if(n != src->size)
	{
		fprintf(stderr, "dest and src vectors length mismatch.\n");
		exit(1);
	}

	for(i=0; i<n; i++)
	{
		c = gsl_vector_char_get(src, i);
		gsl_vector_set(dest, i, (double)c);

		total += (size_t)c;
	}

	return total;
}

size_t matrix_char_to_double(gsl_matrix *dest, const gsl_matrix_char *src)
{
	char c;

	size_t i, j, n = dest->size1, p = dest->size2, total = 0;

	if(n != src->size1 || p != src->size2)
	{
		fprintf(stderr, "dest and src matrix length mismatch.\n");
		exit(1);
	}

	for(i=0; i<n; i++)
	{
		for(j=0; j<p; j++)
		{
			c = gsl_matrix_char_get(src, i, j);
			gsl_matrix_set(dest, i, j, (double)c);

			total += (size_t)c;
		}
	}

	return total;
}


void check_vector_dimension( const gsl_vector *X, const size_t n, const char *title, const int equal_stride)
{
	if(X->size != n)
	{
		fprintf(stderr, "Vector %s length is %lu but not expected %lu.\n", title, X->size, n);
		exit(1);
	}

	if(equal_stride != 0 && X->stride != 1)
	{
		fprintf(stderr, "Vector %s stride %lu is not equal to 1.\n", title, X->stride);
		exit(1);
	}
}

void check_matrix_dimension(const gsl_matrix *X, const size_t n, const size_t p,
		const char *title, const int equal_stride)
{
	if(X->size1 != n || X->size2 != p)
	{
		fprintf(stderr, "Matrix %s dimension is %lu * %lu but not expected %lu * %lu.\n", title, X->size1, X->size2, n, p);
		exit(1);
	}

	if(equal_stride != 0 && X->size2 != X->tda)
	{
		fprintf(stderr, "Matrix %s tda %lu is not equal to column %lu.\n", title, X->tda, X->size2);
		exit(1);
	}
}


void check_vector_category(const gsl_vector_char *X, char catarr[], size_t n, const char *title)
{
	char c;
	size_t i,k;

	for(i=0;i<X->size;i++)
	{
		c = gsl_vector_char_get(X,i);

		for(k=0;k<n;k++)
		{
			if(c == catarr[k]) break;
		}

		if(k==n)
		{
			fprintf(stderr, "Categorical vector %s contains char %c\n", title, c);
			exit(1);
		}
	}
}

void check_matrix_category(const gsl_matrix_char *X, char catarr[], size_t n, const char *title)
{
	char c;
	size_t i,j,k;

	for(i=0;i<X->size1;i++)
	{
		for(j=0;j<X->size2;j++)
		{
			c = gsl_matrix_char_get(X,i,j);

			for(k=0;k<n;k++)
			{
				if(c == catarr[k]) break;
			}

			if(k==n)
			{
				fprintf(stderr, "Categorical matrix %s contains char %c\n", title, c);
				exit(1);
			}
		}
	}
}


gsl_matrix *convert_gsl_matrix(double *data, const size_t n, const size_t p)
{
	gsl_block *block = (gsl_block*)malloc(sizeof(gsl_block));
	gsl_matrix *result = (gsl_matrix*)malloc(sizeof(gsl_matrix));

	result->block = block;
	block->data = result->data = data;
	result->size1 = n;
	result->size2 = result->tda = p;
	result->owner = 1;
	block->size = n*p;

	return result;
}

gsl_matrix_char *convert_gsl_matrix_char(char *data, const size_t n, const size_t p)
{
	gsl_block_char *block = (gsl_block_char*)malloc(sizeof(gsl_block_char));
	gsl_matrix_char *result = (gsl_matrix_char*)malloc(sizeof(gsl_matrix_char));

	result->block = block;
	block->data = result->data = data;
	result->size1 = n;
	result->size2 = result->tda = p;
	result->owner = 1;
	block->size = n*p;

	return result;
}
