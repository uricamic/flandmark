
all:
	gcc -lm -msse -O3 -fPIC -o qpsplx qpsplx.c mtf.c ../lib/libqp_splx.c
	gcc -o qpsplx_compsol mtf.c qpsplx_compsol.c

