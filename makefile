LDC := ldc
DMD := dmd
GSL := ${shell pwd}/gsl
DSOURCES := ${wildcard src/*.d}

CHECK_GSL := ${shell gsl-config --version | awk '$$0 >= 2.1' 2> /dev/null}
CHECK_LDC := ${shell command -v ${LDC} 2> /dev/null}
CHECK_DMD := ${shell command -v ${DMD} 2> /dev/null}

ifneq (${CHECK_GSL},)
	C_SOURCES := src/interpolate.o
else
	GSL_FILES := ${GSL}/lib/libgsl.a ${GSL}/lib/libgslcblas.a
	C_SOURCES := src/static_interpolate.o ${GSL_FILES}
	CHECK_LOCAL_GSL := ${shell ${foreach f,${GSL_FILES},test -e "$f" || echo "$f";}}
endif

ifneq (${CHECK_LDC},)
	COMPILER := ${LDC}
	RELEASE_FLAGS := -Jviews -release -enable-inlining -O -w -oq
	DEBUG_FLAGS := -Jviews -d-debug -g -unittest -w
ifeq (${CHECK_GSL},)
	STATIC_FLAGS := -d-version=STATICLINKED -I${GSL}/include
endif
else
	COMPILER := ${DMD}
	RELEASE_FLAGS := -Jviews -release -inline -O -noboundscheck
	DEBUG_FLAGS := -Jviews -debug -g -unittest -w
ifeq (${CHECK_GSL},)
	STATIC_FLAGS := -version=STATICLINKED -I${GSL}/include
endif
endif

ifeq (${CHECK_LDC},)
ifeq (${CHECK_DMD},)
${error No D compiler found at ${LDC} or ${DMD}}
endif
endif

ifeq (${CHECK_GSL},)
ifneq (${CHECK_LOCAL_GSL},)
${error GSL not installed at ${GSL}}
endif
endif

CLEAN_OBJECTS := rm -f src/*.o bin/*.o *.o

bin/CaVEMaN	: ${DSOURCES} ${C_SOURCES} views/commit
	${COMPILER} ${RELEASE_FLAGS} ${DSOURCES} ${C_SOURCES} ${STATIC_FLAGS} -ofbin/CaVEMaN
	${CLEAN_OBJECTS}

test	: ${DSOURCES} ${C_SOURCES} views/commit
	${COMPILER} ${DEBUG_FLAGS} ${DSOURCES} ${C_SOURCES} ${STATIC_FLAGS} -ofunittest
	./unittest
	${CLEAN_OBJECTS} unittest

src/static_interpolate.o : src/interpolate.c ${GSL}/lib/libgsl.a
	cc -c src/interpolate.c -o src/static_interpolate.o -I${GSL}/include

views/commit : ${DSOURCES} src/interpolate.c
	git rev-parse --short HEAD > views/commit

.PHONY : test clean

clean :
	${CLEAN_OBJECTS} bin/CaVEMaN
