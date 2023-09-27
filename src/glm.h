/*
 * glm.h
 *
 *  Created on: Apr 11, 2014
 *      Author: Peng Jiang
 *  Linear models such as OLS, logistic regression
 */

#ifndef GLM_H_
#define GLM_H_

#include <gsl/gsl_blas.h>

// all fail flags must have value < 0
#define REG_FAIL -1
#define CONVERGE_FAIL -2


// Ordinary least square on Y as matrix: n*k
int OLS(const gsl_matrix *X, const gsl_matrix *Y, gsl_matrix *beta);


// logistic regression with Y as vector. Return regression status
int logistic_regression(const gsl_matrix *X, const gsl_vector_char *Y,
		const size_t maxIter, const double tol, const double max_delta,
		const int correction, const int verbose,
		gsl_vector *beta, gsl_vector *sderr, gsl_vector *z, gsl_vector *pvalue);


// logistic regression with Y as matrix: n*k. result beta, sderr, z, pvalue as matrix k * p
// Return number of successful regressions
size_t logistic_regression_batch(const gsl_matrix *X, const gsl_matrix_char *Y,
		const size_t maxIter, const double tol, const double max_delta, const size_t cntthres,
		const int correction, const int verbose,
		gsl_matrix *beta, gsl_matrix *sderr, gsl_matrix *z, gsl_matrix *pvalue);


// some tools for debugging
void print_vector(const gsl_vector *v, const char *title, FILE *fp);
void print_matrix(const gsl_matrix *X, const char *title, FILE *fp);

// convert categorical chars to double numerical
size_t vector_char_to_double(gsl_vector *dest, const gsl_vector_char *src);
size_t matrix_char_to_double(gsl_matrix *dest, const gsl_matrix_char *src);


// check if vector X is n, if equal_stride is not 0, also check stride == 1.
void check_vector_dimension( const gsl_vector *X, const size_t n,
		const char *title, const int equal_stride);

// check if matrix X is n*p, if equal_stride is not 0, also check stride == size2.
void check_matrix_dimension( const gsl_matrix *X, const size_t n, const size_t p,
		const char *title, const int equal_stride);

// check categorical vector or matrix only contains allowed characters in catarr[] of length n
void check_vector_category(const gsl_vector_char *X, char catarr[], size_t n, const char *title);
void check_matrix_category(const gsl_matrix_char *X, char catarr[], size_t n, const char *title);


// convert some data array to gsl matrix
gsl_matrix *convert_gsl_matrix(double *data, const size_t n, const size_t p);
gsl_matrix_char *convert_gsl_matrix_char(char *data, const size_t n, const size_t p);

#endif /* GLM_H_ */
