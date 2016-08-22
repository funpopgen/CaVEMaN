#Replace with location of locally installed gsl-2.1 versions if compiling static versions
GSL = /home/andrew/gsl/
DSOURCES = src/main.d src/arg_parse.d src/read_data.d src/calculation.d src/run_analysis.d src/correct.d src/best.d src/interval.d
LDC = ldc
DMD = dmd

ldc : ${DSOURCES} src/interpolate.o
	${LDC} -Jviews -release -enable-inlining -O -w -oq ${DSOURCES} src/interpolate.o -of="bin/CaVEMaN"
	rm -f *.o
	rm -f src/*.o

test : ${DSOURCES} src/interpolate.o
	${LDC} -Jviews -d-debug -g -unittest -w ${DSOURCES} src/interpolate.o -of="unittest"
	./unittest
	rm -f unittest src/*.o *.o

static : ${DSOURCES} src/interpolate.c
	cc -c src/interpolate.c -o src/interpolate.o -I${GSL}/include
	${LDC} -Jviews -release -enable-inlining -O -w -oq -d-version=STATICLINKED -I${GSL}/include ${DSOURCES} src/interpolate.o ${GSL}/lib/libgsl.a ${GSL}/lib/libgslcblas.a -of="bin/CaVEMaN"
	rm -f *.o
	rm -f src/*.o

static_test : ${DSOURCES} src/interpolate.c
	cc -c src/interpolate.c -o src/interpolate.o -I${GSL}/include
	${LDC} -Jviews -d-debug -g -unittest -w -d-version=STATICLINKED -I${GSL}/include ${DSOURCES} src/interpolate.o ${GSL}/lib/libgsl.a ${GSL}/lib/libgslcblas.a -of="unittest"
	./unittest
	rm -f unittest src/*.o *.o

dmd : ${DSOURCES} src/interpolate.o
	${DMD} -Jviews -O -release -noboundscheck -inline ${DSOURCES} src/interpolate.o -ofbin/CaVEMaN
	rm src/*.o

dmd_test : ${DSOURCES} src/interpolate.o
	${DMD} -Jviews -debug -g -unittest -w ${DSOURCES} src/interpolate.o -ofunittest
	./unittest
	rm -f unittest src/*.o *.o

dmd_static : ${DSOURCES} src/interpolate.c
	cc -c src/interpolate.c -o src/interpolate.o -I${GSL}/include
	${DMD} -Jviews -O -release -noboundscheck -inline -version=STATICLINKED -I${GSL}/include ${GSL}/lib/libgsl.a ${GSL}/lib/libgslcblas.a ${DSOURCES} src/interpolate.o -ofbin/CaVEMaN
	rm src/*.o

dmd_static_test : ${DSOURCES} src/interpolate.c
	cc -c src/interpolate.c -o src/interpolate.o -I${GSL}/include
	${DMD} -Jviews -debug -g -unittest -w -version=STATICLINKED -I${GSL}/include ${GSL}/lib/libgsl.a ${GSL}/lib/libgslcblas.a ${DSOURCES} src/interpolate.o -ofunittest
	./unittest
	rm -f unittest src/*.o *.o

.PHONY : test static ldc dmd ldc_test dmd_test static static_test dmd_static dmd_static_test clean install uninistall

clean :
	rm -f src/*.o bin/CaVEMaN *.o

install : ${DSOURCES} CaVEMaN.1
	cp -v ${shell pwd}/bin/CaVEMaN /usr/local/bin/
	cp -v ${shell pwd}/CaVEMaN.1 /usr/local/man/man1/

uninstall : ${DSOURCES} CaVEMaN.1
	rm -v /usr/local/bin/CaVEMaN
	rm -v /usr/local/man/man1/CaVEMaN.1
