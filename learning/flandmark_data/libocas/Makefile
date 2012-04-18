MEXDETECTED := $(shell mex -h >/dev/null 2>&1 && echo yes)
MEX := mex
MEXFLAGS := $(shell if uname -m | grep -q x86_64 ; then echo -largeArrayDims ; fi)
MEXSUFFIX := $(shell if uname -m | grep -q x86_64 ; then echo mexa64 ; else echo mexglx ; fi)
CC := gcc
#CFLAGS := -g -lm -Wall -pthread
CFLAGS := -lm -msse -O3 -fPIC -fstrict-aliasing -fomit-frame-pointer -Wall -pthread
#CFLAGS := -lm -msse -O3 -fPIC -fstrict-aliasing -fomit-frame-pointer -Wall
#CFLAGS := -lm -msse -O3 -fPIC -fopenmp -fstrict-aliasing -fomit-frame-pointer -Wall


ifeq (yes,$(MEXDETECTED))
all: svmocas.$(MEXSUFFIX) svmlight_linclass.$(MEXSUFFIX) libocas.so svmocas msvmocas linclass msvmocas.$(MEXSUFFIX) compute_errors.$(MEXSUFFIX) lbppyr.$(MEXSUFFIX) svmocas_lbp.$(MEXSUFFIX) lbppyr_features.$(MEXSUFFIX) lbpfilter.$(MEXSUFFIX)

compute_errors.$(MEXSUFFIX): compute_errors_mex.c libocas.h ocas_helper.h ocas_helper.c
		$(MEX) $(MEXFLAGS) -DLIBOCAS_MATLAB -O -output compute_errors.$(MEXSUFFIX) compute_errors_mex.c ocas_helper.c lib_svmlight_format.c

svmocas.$(MEXSUFFIX): libocas.c libocas.h libqp_splx.c libqp.h svmocas_mex.c lib_svmlight_format.c lib_svmlight_format.h ocas_helper.c ocas_helper.h
		$(MEX) $(MEXFLAGS) -DLIBOCAS_MATLAB -O -output svmocas.$(MEXSUFFIX) svmocas_mex.c lib_svmlight_format.c ocas_helper.c  libocas.c libqp_splx.c
		
		
linclass_light.$(MEXSUFFIX): linclass_light_mex.c lib_svmlight_format.c lib_svmlight_format.h ocas_helper.c ocas_helper.h
		$(MEX) $(MEXFLAGS) -DLIBOCAS_MATLAB -O -output linclass_light.$(MEXSUFFIX) linclass_light_mex.c lib_svmlight_format.c ocas_helper.c

		
svmlight_linclass.$(MEXSUFFIX): svmlight_linclass_mex.c lib_svmlight_format.c lib_svmlight_format.h ocas_helper.c ocas_helper.h
		$(MEX) $(MEXFLAGS) -DLIBOCAS_MATLAB -O -output svmlight_linclass.$(MEXSUFFIX) svmlight_linclass_mex.c lib_svmlight_format.c ocas_helper.c
		
msvmocas.$(MEXSUFFIX): libocas.c libocas.h libqp_splx.c libqp.h msvmocas_mex.c lib_svmlight_format.c lib_svmlight_format.h ocas_helper.c ocas_helper.h
		$(MEX) $(MEXFLAGS) -DLIBOCAS_MATLAB -O -output msvmocas.$(MEXSUFFIX) msvmocas_mex.c lib_svmlight_format.c ocas_helper.c libocas.c libqp_splx.c
		
lbppyr.$(MEXSUFFIX): lbppyr_mex.c
		$(MEX) $(MEXFLAGS) -O -output lbppyr.$(MEXSUFFIX) lbppyr_mex.c 	

svmocas_lbp.$(MEXSUFFIX): libocas.c libocas.h libqp_splx.c libqp.h svmocas_lbp_mex.c ocas_lbp_helper.c ocas_lbp_helper.h liblbp.h liblbp.c
		$(MEX) $(MEXFLAGS) -DLIBOCAS_MATLAB -O -output svmocas_lbp.$(MEXSUFFIX) svmocas_lbp_mex.c ocas_lbp_helper.c  libocas.c libqp_splx.c liblbp.c

lbppyr_features.$(MEXSUFFIX): lbppyr_features_mex.c liblbp.c liblbp.h
		$(MEX) -g $(MEXFLAGS) -DLIBOCAS_MATLAB -O -output lbppyr_features.$(MEXSUFFIX) lbppyr_features_mex.c liblbp.c

lbpfilter.$(MEXSUFFIX): lbpfilter_mex.c
		$(MEX) $(MEXFLAGS) -DLIBOCAS_MATLAB -O -output lbpfilter.$(MEXSUFFIX) lbpfilter_mex.c


else 
all: libocas.so svmocas msvmocas
endif

svmocas:	svmocas.c lib_svmlight_format.c sparse_mat.c ocas_helper.c ocas_helper.h libocas.h sparse_mat.h libocas.c
		$(CC) $(CFLAGS) -o $@ svmocas.c lib_svmlight_format.c sparse_mat.c ocas_helper.c libocas.c libqp_splx.c

msvmocas:	msvmocas.c lib_svmlight_format.c sparse_mat.c ocas_helper.c ocas_helper.h libocas.h sparse_mat.h libocas.c
		$(CC) $(CFLAGS) -o $@ msvmocas.c lib_svmlight_format.c sparse_mat.c ocas_helper.c libocas.c libqp_splx.c

linclass:	linclass.c lib_svmlight_format.c libocas.h 
		$(CC) $(CFLAGS) -o $@ linclass.c lib_svmlight_format.c  


libocas.so:	libocas.c libocas.h libqp_splx.c libqp.h 
		$(CC) $(CFLAGS) -shared -o $@ libocas.c libqp_splx.c

clean: 
		rm -f *~ svmocas.$(MEXSUFFIX) svmlight_linclass.$(MEXSUFFIX) svmocas msvmocas linclass libocas.so msvmocas.$(MEXSUFFIX) lbpfilter.$(MEXSUFFIX) lbppyr_features.$(MEXSUFFIX) svmocas_lbp.$(MEXSUFFIX) lbppyr.$(MEXSUFFIX) compute_errors.$(MEXSUFFIX)
