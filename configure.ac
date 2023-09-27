AC_INIT(logis_batch, 1.0, pengj@alumni.princeton.edu)

AM_INIT_AUTOMAKE()

AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB

#GSL check
AC_CHECK_HEADERS([gsl/gsl_blas.h gsl/gsl_cdf.h], ,AC_MSG_ERROR([Cannot find gsl headers]))
AC_CHECK_LIB(gslcblas, cblas_dgemm, ,AC_MSG_ERROR([Cannot find cblas_dgemm]))
AC_CHECK_LIB(gslcblas, cblas_dgemv, ,AC_MSG_ERROR([Cannot find cblas_dgemv]))
AC_CHECK_LIB(gslcblas, cblas_dsyrk, ,AC_MSG_ERROR([Cannot find cblas_dsyrk]))
AC_CHECK_LIB(gsl, gsl_cdf_gaussian_P, ,AC_MSG_ERROR([Cannot find gsl_cdf_gaussian_P]), [-lm])
AC_CHECK_LIB(gsl, gsl_linalg_cholesky_decomp, ,AC_MSG_ERROR([Cannot find gsl_linalg_cholesky_decomp]), [-lm])
AC_CHECK_LIB(gsl, gsl_linalg_cholesky_invert, ,AC_MSG_ERROR([Cannot find gsl_linalg_cholesky_invert]), [-lm])
AC_CHECK_LIB(gsl, gsl_linalg_cholesky_solve, ,AC_MSG_ERROR([Cannot find gsl_linalg_cholesky_solve]), [-lm])

# disable default LIBS
AC_SUBST(LIBS, "")

AC_CONFIG_FILES(Makefile src/Makefile)
AC_OUTPUT
