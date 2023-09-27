# Bias reduced logistic regression.  
The logistic regression may fail when outcome variable separates a predictor variable or a combination of predictor variables almost completely (quasi-complete separation). Also, a user may run multiple logistic regressions with a common set of covariates and many outcome labels. Thus, we release this package for running multiple logistic regressions with the Firth correction.  

**Installation**: Please read the INSTALL file  

**Test**: After compiling with configure & make, please run:  
make check  

**Usage**: logis_batch -B background_matrix -X X_matrix -Y Y_matrix -out output [OPTIONS]  
For example, in the project root folder:  
cd data  
logis_batch -B MSigDB.Hallmark.background -X CytoSig.signature -Y MSigDB.Hallmark.mat -out test.output  

**Input**:  
*X_matrix*: matrix of covariates Xi, to be tested individually.  
*Y_matrix*: matrix of binary outcome variables Yj, to be tested individually.  
*background_matrix*: background covariates to be controlled in the regression. Please do not include the intercept column, as the program will automatically append a column of 1.  

**Algorithm**:  
For each column Xi in the X_matrix, this program will first concatenate Xi with background_matrix to create X = (1, background_matrix, Xi), and run one regression in combination with every column Yj from the Y_matrix.  

**Output**:  
*output.coef, output.zscore, output.pvalue*: coefficients, Wald test z-scores, and two-sided p-values from logistic regressions. In each matrix, rows represent variables in the X_matrix, and columns represent variables in the Y_matrix.  

**Options**:  
  
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
