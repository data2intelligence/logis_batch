extern "C" {
#include "glm.h"
}

#include "util.h"

#include <ctime>
#include <cmath>
#include <cstring>
#include <sstream>
#include <limits>
using namespace std;

#define EPS 1e-10

int main(int argc, char *argv[])
{
	string
		Xfile,	// X matrix file in regression
		Yfile,	// Y response matrix file
		Bfile,	// Background confounding factor matrix
		Sfile,	// Coefficient start matrix in Newton Raphson process
		output, value, type;

	size_t i,j, n,p,k, nfactor, nonzero_cnt,
		parseCnt = (argc - 1)/2,	// number of input parameters
		maxIter = 1000,	// maximum number of Newton Raphson iteration
		cntthres = 5;	// filter over rare Y out comes

	// Newton Raphson Convergence parameter
	double tol = 1e-5, max_delta = 1;

	// 0: regular, 1: Firth correction
	int correction = 1;

	bool filter_X = true,		// apply count threshold on X columns
		intercept_flag = true,	// use intercept in regression
		background_flag = false,// only profile background factors for incremental start
		verbose_flag = false;	// whether print details

	////////////////////////////////////////////////////////////////////////////////////////
	// Part 0: parameter input and check
	if ( argc < 7 )
	{
		if(argc == 2 && (value = argv[1]) == "-help")
		{
			cout << "\nBias reduced incremental logistic regression. (Peng Jiang 2014)\n" << endl;
			cout << "Usage: logis_batch -B background_matrix -X X_matrix -Y Y_matrix -out output [OPTIONS]\n"<<endl;

			cout << "Options:" <<endl;

			cout << "\n\tNetwon Raphson process:" <<endl;
			cout << "\t-S\t\tStart coefficients in iteration. Default: All zeros "<< endl;
			cout << "\t-maxiter\tMaximum number of iterations. Default: " << maxIter <<endl;
			cout << "\t-tol\t\tConvergence tolerance. Default: " << tol <<endl;
			cout << "\t-delta\t\tRelative radius of trust region. Default: " << max_delta <<endl;
			cout << "\t-correction\tFirth correction. Default: " << correction <<endl;

			cout << "\n\tOthers:" <<endl;
			cout << "\t-cntthres\tSkip over rare Y response. Default: " << cntthres <<endl;
			cout << "\t-filterX\tApply cntthres to filter X columns with too many zeros. Default: " << DISPLAY_BOOL(filter_X) <<endl;
			cout << "\t-intercept\tUse intercept in regression. Default: " << DISPLAY_BOOL(intercept_flag) <<endl;
			cout << "\t-background\tBuild background start. Default: " << DISPLAY_BOOL(background_flag) <<endl;
			cout << "\t-verbose\tPrint details. Default: " << DISPLAY_BOOL(verbose_flag) <<endl;

			cout << "\nReport bugs to peng.jiang.software@gmail.com\n" <<endl;
			exit(0);
		}else{
			cerr << "Insufficient number of arguments, do \"logis_batch -help\" for help."<<endl;
			exit(1);
		}
	}

	// read in all parameters
	for(i=0;i<parseCnt;i++)
	{
		type = argv[2*i+1];
		value = argv[2*i+2];

		if(type == "-X"){
			Xfile = value;
		}else if (type == "-Y"){
			Yfile = value;
		}else if (type == "-B"){
			Bfile = value;
		}else if (type == "-S"){
			Sfile = value;
		}else if (type == "-out"){
			output = value;

		}else if (type == "-maxiter"){
			maxIter = load_positive_value(value, type, (size_t)1, string::npos);

		}else if (type == "-tol"){
			tol = load_positive_value(value, type, 0.0, 1e-3);

		}else if (type == "-delta"){
			max_delta = load_positive_value(value, type, 0.0, numeric_limits<double>::max());

		}else if (type == "-correction"){
			correction = load_positive_value(value, type, (int)0, (int)1);

		}else if (type == "-cntthres"){
			cntthres = load_positive_value(value, type, (size_t)1, string::npos);

		}else if (type == "-filterX"){
			filter_X = (value[0]!='0');

		}else if (type == "-intercept"){
			intercept_flag = (value[0]!='0');

		}else if (type == "-background"){
			background_flag = (value[0]!='0');

		}else if (type == "-verbose"){
			verbose_flag = (value[0]!='0');

		}else if (type == "-help"){
			cerr << "Please don't use \"-help\" as parameter input." << endl;
			exit(1);

		}else{
			cerr << "Cannot recognize parameter \""<< type << "\"." << endl;
			exit(1);
		}
	}

	if(Bfile.empty())
	{
		cerr << "Cannot find background file." << endl;
		exit(1);
	}

	if(Xfile.empty())
	{
		cerr << "Cannot find X file." << endl;
		exit(1);
	}

	if(Yfile.empty())
	{
		cerr << "Cannot find Y file." << endl;
		exit(1);
	}

	if(output.empty())
	{
		cerr << "Cannot find output file." << endl;
		exit(1);
	}


	////////////////////////////////////////////////////////////////////////////////////////
	// Part 1: Matrix initialization

	vector<string> *left_pad = NULL, *right_pad = NULL;

	if(intercept_flag)
	{
		left_pad = new vector<string>();
		left_pad->push_back("Intercept");
	}

	// will put extra column for regulators
	if(background_flag == false)
	{
		right_pad = new vector<string>();
		right_pad->push_back("Regulator");
	}

	string_array
		*gene_names_X, *gene_names_Y, *gene_names_B,
		*variable_names, *sample_names, *factor_names;

	double *X_data, *B_data, *start_data;
	char *Y_data;
	vector<string> gene_names_common;

	read_matrix(Bfile, B_data, gene_names_B, variable_names, true, 0.0, left_pad, right_pad, 1.0);
	read_matrix(Xfile, X_data, gene_names_X, factor_names, true, 0.0);
	read_matrix(Yfile, Y_data, gene_names_Y, sample_names, true, '0');

	common_names(gene_names_X, gene_names_Y, gene_names_B, gene_names_common);

	if( gene_names_common.size() < cntthres)
	{
		cerr << "Number of overlap genes " << gene_names_common.size() << " is too small." << endl;
		exit(1);
	}

	compact_matrix(gene_names_X, X_data, factor_names->n, gene_names_common, false);
	compact_matrix(gene_names_Y, Y_data, sample_names->n, gene_names_common, true);
	compact_matrix(gene_names_B, B_data, variable_names->n, gene_names_common, false);

	nfactor = factor_names->n;
	n = gene_names_Y->n;
	p = variable_names->n;
	k = sample_names->n;

	gsl_matrix
		*X = convert_gsl_matrix(X_data, n, nfactor),
		*B = convert_gsl_matrix(B_data, n, p),
		*start = NULL,
		*beta = gsl_matrix_calloc(k, p),	// default to zeros
		*sderr = gsl_matrix_alloc(k, p),
		*z = gsl_matrix_alloc(k, p),
		*pvalue = gsl_matrix_alloc(k, p),

		// final result matrix to hold the last factor column result
		*beta_factor = gsl_matrix_calloc(nfactor, k),
		*z_factor = gsl_matrix_calloc(nfactor, k),
		*pvalue_factor = gsl_matrix_alloc(nfactor, k);

	gsl_matrix_char *Y = convert_gsl_matrix_char(Y_data, n, k);
	gsl_vector_view B_c = gsl_matrix_column(B, p-1);


	////////////////////////////////////////////////////////////////////////////////////////
	// Part 2: Logistic regression: start making or whole data set run

	clock_t startClock;

	if(background_flag)
	{	// only profile background coef
		if(verbose_flag) cout << "Profile background coefficients" << endl;

		logistic_regression_batch(B, Y, maxIter, tol, max_delta, cntthres, correction, verbose_flag, beta, sderr, z, pvalue);

		write_matrix(output + ".coef", beta->data, sample_names, variable_names, true);
		write_matrix(output + ".zscore", z->data, sample_names, variable_names, true);
		write_matrix(output + ".pvalue", pvalue->data, sample_names, variable_names, true);

	}else
	{
		if(Sfile.empty())
		{
			if(verbose_flag) cout << "Cannot find start file, use all zeros" << endl;

		}else{
			// load start for incremental regression
			if(verbose_flag) cout << "Profile incremental coefficients" << endl;

			free_string_array(sample_names);
			free_string_array(variable_names);

			read_matrix(Sfile, start_data, sample_names, variable_names, true, 0.0, NULL, right_pad, 0.0);

			if(beta->size1 != sample_names->n || beta->size2 != variable_names->n)
			{
				cerr << "Start matrix dimension "<<  sample_names->n << ',' << variable_names->n
					<< " is different from beta dimension " <<  beta->size1 << ',' << beta->size2 << endl;

				exit(1);
			}

			start = convert_gsl_matrix(start_data, k, p);
		}

		gsl_vector_view beta_last_c =	gsl_matrix_column(beta, p-1);
		gsl_vector_view z_last_c =		gsl_matrix_column(z, p-1);
		gsl_vector_view pvalue_last_c =	gsl_matrix_column(pvalue, p-1);

		for(i=0; i < nfactor; i++)
		{
			if(verbose_flag) cout << factor_names->array[i] << endl;

			gsl_vector_view X_c = gsl_matrix_column(X, i);
			gsl_vector_view beta_factor_r =		gsl_matrix_row(beta_factor, i);
			gsl_vector_view z_factor_r =		gsl_matrix_row(z_factor, i);
			gsl_vector_view pvalue_factor_r =	gsl_matrix_row(pvalue_factor, i);

			for(nonzero_cnt=0,j=0 ; j<X_c.vector.size ; j++)
			{
				if (fabs(gsl_vector_get(&X_c.vector, j)) > EPS) nonzero_cnt ++;
			}

			if(filter_X && nonzero_cnt < cntthres)
			{
				gsl_vector_set_zero(&beta_factor_r.vector);
				gsl_vector_set_zero(&z_factor_r.vector);
				gsl_vector_set_all(&pvalue_factor_r.vector, 1);
				continue;
			}

			gsl_vector_memcpy(&B_c.vector, &X_c.vector);

			if(start!=NULL){
				gsl_matrix_memcpy(beta, start);
			}else{
				gsl_matrix_set_zero(beta);
			}

			startClock = clock();
			logistic_regression_batch(B, Y, maxIter, tol, max_delta, cntthres, correction, verbose_flag, beta, sderr, z, pvalue);

			if(verbose_flag) {
				cout << (double)(clock() - startClock) / CLOCKS_PER_SEC << "s elapsed." << endl;
			}

			gsl_vector_memcpy(&beta_factor_r.vector,		&beta_last_c.vector);
			gsl_vector_memcpy(&z_factor_r.vector,		&z_last_c.vector);
			gsl_vector_memcpy(&pvalue_factor_r.vector,	&pvalue_last_c.vector);
		}

		write_matrix(output + ".coef",	beta_factor->data,	factor_names, sample_names, true);
		write_matrix(output + ".zscore", z_factor->data,		factor_names, sample_names, true);
		write_matrix(output + ".pvalue",	pvalue_factor->data,	factor_names, sample_names, true);
	}


	// free all memory
	gsl_matrix_char_free(Y);
	gsl_matrix_free(X);
	gsl_matrix_free(beta);
	gsl_matrix_free(sderr);
	gsl_matrix_free(z);
	gsl_matrix_free(pvalue);
	gsl_matrix_free(beta_factor);
	gsl_matrix_free(z_factor);
	gsl_matrix_free(pvalue_factor);

	if(start != NULL) gsl_matrix_free(start);

	free_string_array(gene_names_X);
	free_string_array(gene_names_Y);
	free_string_array(gene_names_B);

	free_string_array(sample_names);
	free_string_array(variable_names);
	free_string_array(factor_names);

	if (left_pad != NULL)  delete left_pad;
	if (right_pad != NULL) delete right_pad;

	return 0;
}
