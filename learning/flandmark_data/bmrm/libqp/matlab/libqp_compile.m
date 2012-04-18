% This script compiles Matlab's MEX interfaces for LIBQP solvers.
%


mex libqp_splx_mex.c -largeArrayDims ../lib/libqp_splx.c
mex libqp_gsmo_mex.c -largeArrayDims ../lib/libqp_gsmo.c

