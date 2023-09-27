# Bias reduced logistic regression.  
Logistic regression may fail when quasi-complete separation happens, when outcome variable separates a predictor variable or a combination of predictor variables almost completely. Also, a user may want to run multiple logistic regressions with a common set of covariates and many outcome labels. Thus, we release this package for running multiple logistic regressions with the Firth correction.  

Usage: logis_batch -B background_matrix -X X_matrix -Y Y_matrix -out output [OPTIONS]  

Options:  
  
	Netwon Raphson process:  
	-S		Start coefficients in iteration. Default: All zeros   
	-maxiter	Maximum number of iterations. Default: 1000  
	-tol		Convergence tolerance. Default: 1e-05  
	-delta		Relative radius of trust region. Default: 1  
	-correction	Firth correction. Default: 1  
  
	Others:  
	-cntthres	Skip over rare Y response. Default: 5  
	-filterX	Apply cntthres to filter X columns with too many zeros. Default: 1 (yes)  
	-intercept	Use intercept in regression. Default: 1 (yes)  
	-background	Build background start. Default: 0 (no)  
	-verbose	Print details. Default: 0 (no)  
  
Report bugs to peng.jiang.software@gmail.com  
