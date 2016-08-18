#Replace with location of locally installed gsl-2.1 versions if compiling cluster versions
GSL = /Home/abrown/software/gsl-2.1-lib/
DSOURCES = src/main.d src/arg_parse.d src/read_data.d src/calculation.d src/run_analysis.d src/correct.d src/best.d src/interval.d
LDC = ldc
DMD = dmd

ldc : ${DSOURCES} src/interpolate.o
	${LDC} -Jviews -release -enable-inlining -O -w -oq ${DSOURCES} src/interpolate.o -L-lgsl -L-lgslcblas -of="bin/CaVEMaN"
	rm -f *.o
	rm -f src/*.o

test : ${DSOURCES} src/interpolate.o
	${LDC} -Jviews -d-debug -g -unittest -w -L-lgsl -L-lgslcblas ${DSOURCES} src/interpolate.o -of="unittest"
	./unittest
	rm -f unittest src/*.o *.o

cluster : ${DSOURCES} src/interpolate.c
	cc -c src/interpolate.c -o src/interpolate.o -I${GSL}/include
	${LDC} -Jviews -release -enable-inlining -O -w -oq -I${GSL}/include ${DSOURCES} src/interpolate.o ${GSL}/lib/libgsl.a ${GSL}/lib/libgslcblas.a -of="bin/CaVEMaN"
	rm -f *.o
	rm -f src/*.o

cluster_test : ${DSOURCES} src/interpolate.c
	cc -c src/interpolate.c -o src/interpolate.o -I${GSL}/include
	${LDC} -Jviews -d-debug -g -unittest -w -I${GSL}/include ${DSOURCES} src/interpolate.o ${GSL}/lib/libgsl.a ${GSL}/lib/libgslcblas.a -of="unittest"
	./unittest
	rm -f unittest src/*.o *.o

dmd : ${DSOURCES} src/interpolate.o
	${DMD} -Jviews -O -release -noboundscheck -inline -L-lgsl -L-lgslcblas -I${GSL}/include ${DSOURCES} src/interpolate.o -ofbin/CaVEMaN
	rm src/*.o

dmd_test : ${DSOURCES} src/interpolate.o
	${DMD} -Jviews -debug -g -unittest -w -L-lgsl -L-lgslcblas ${DSOURCES} src/interpolate.o -ofunittest
	./unittest
	rm -f unittest src/*.o *.o

dmd_cluster : ${DSOURCES} src/interpolate.c
	cc -c src/interpolate.c -o src/interpolate.o -I${GSL}/include
	${DMD} -Jviews -O -release -noboundscheck -inline -I${GSL}/include ${GSL}/lib/libgsl.a ${GSL}/lib/libgslcblas.a ${DSOURCES} src/interpolate.o -ofbin/CaVEMaN
	rm src/*.o

dmd_test_cluster : ${DSOURCES} src/interpolate.c
	cc -c src/interpolate.c -o src/interpolate.o -I${GSL}/include
	${DMD} -Jviews -debug -g -unittest -w -I${GSL}/include ${GSL}/lib/libgsl.a ${GSL}/lib/libgslcblas.a ${DSOURCES} src/interpolate.o -ofunittest
	./unittest
	rm -f unittest src/*.o *.o

.PHONY : test static ldc dmd ldc_test dmd_test cluster cluster_test dmd_cluster dmd_test_cluster clean install uninistall

probability : utilities/probability.d src/interpolate.o
	${LDC} -release -enable-inlining -O -w -oq utilities/probability.d src/interpolate.o -L-lgsl -L-lgslcblas -of="bin/probability"
	rm -f *.o
	rm -f src/*.o

clean :
	rm -f src/*.o bin/CaVEMaN *.o

install : ${DSOURCES} CaVEMaN.1
	cp -v ${shell pwd}/bin/CaVEMaN /usr/local/bin/
	cp -v ${shell pwd}/CaVEMaN.1 /usr/local/man/man1/

uninstall : ${DSOURCES} CaVEMaN.1
	rm -v /usr/local/bin/
	rm -v /usr/local/man/man1/
