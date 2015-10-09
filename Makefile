
CC=mpicc
DEPCMD = gcc -MG -MT .objs/$(<:%.c=%.o) -MM $(CPPFLAGS)

OPTIMIZE ?= -O0 -g -fopenmp

CPPFLAGS += -I lua
LDFLAGS += -L lua

GSL_DIR   ?= #e.g. /Users/jkoda/Research/opt/gcc/gsl

DIR_PATH = $(GSL_DIR) depends/install

CPPFLAGS += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LDFLAGS += $(foreach dir, $(DIR_PATH), -L$(dir)/lib) 

SOURCES = main.c fastpm.c vpm.c pmpfft.c pmghosts.c pmpaint.c pmstore.c pm2lpt.c pmic.c \
		readparams.c msg.c power.c pmsteps.c pmtimer.c pmio-runpb.c cosmology.c
LIBSOURCES = libfastpm.c pmpfft.c pmghosts.c pmpaint.c pmstore.c pm2lpt.c pmic.c \
		msg.c pmsteps.c pmtimer.c cosmology.c

PFFTLIB = depends/install/lib/libpfft_omp.a
PFFTFLIB = depends/install/lib/libpfftf_omp.a
PFFT_VERSION=1.0.8-alpha1-fftw3
PFFT_CONFIGURE_FLAGS = --enable-sse2 --enable-avx

LUALIB = lua/liblua.a

all: fastpm libfastpm.a testlib

fastpm: $(PFFTLIB) $(PFFTFLIB) $(LUALIB) $(SOURCES:%.c=.objs/%.o) 
	$(CC) $(OPTIMIZE) -o fastpm $(SOURCES:%.c=.objs/%.o) \
			$(LDFLAGS) -llua -lgsl -lgslcblas \
			-lpfft_omp -lfftw3_mpi -lfftw3_omp -lfftw3 \
			-lpfftf_omp -lfftw3f_mpi -lfftw3f_omp -lfftw3f \
			-lm 
testlib : testlib.c libfastpm.a
	$(CC) $(OPTIMIZE) $(CPPFLAGS) -o $@ testlib.c libfastpm.a \
			$(LDFLAGS) \
			-lgsl -lgslcblas \
			-lpfft_omp -lfftw3_mpi -lfftw3_omp -lfftw3 \
			-lpfftf_omp -lfftw3f_mpi -lfftw3f_omp -lfftw3f \
			-lm 

libfastpm.a : $(PFFTLIB) $(PFFTFLIB) $(LIBSOURCES:%.c=.objs/%.o)
	$(AR) rcs $@ $(LIBSOURCES:%.c=.objs/%.o) 

$(LUALIB): lua/Makefile
	(cd lua; CC="$(CC)" make generic)

-include $(SOURCES:%.c=.deps/%.d)

.objs/%.o : %.c
	@if ! [ -d .objs ]; then mkdir .objs; fi
	$(CC) $(OPTIMIZE) -c $(CPPFLAGS) -o $@ $<

.deps/%.d : %.c
	@if ! [ -d .deps ]; then mkdir .deps; fi
	@if ! $(DEPCMD) -o $@ $< ; then \
		rm $@; \
	fi;

clean:
	rm -rf .objs
	rm -rf .deps

PFFT_CONFIGURE = $(abspath depends/src/pfft-$(PFFT_VERSION)/configure)
PFFT_CONFIGURE_FLAGS_SINGLE = $(subst --enable-sse2, --enable-sse, $(PFFT_CONFIGURE_FLAGS))
PFFTSRC = pfft-$(PFFT_VERSION).tar.gz 

$(PFFT_CONFIGURE): $(PFFTSRC)
	mkdir -p depends/src ;
	gzip -dc $(PFFTSRC) | tar xf - -C depends/src ;
	touch $@

$(PFFTLIB): $(PFFT_CONFIGURE)
	mkdir -p depends/double;
	(cd depends/double; \
	$(PFFT_CONFIGURE) --prefix=$(abspath depends/install) --disable-shared --enable-static  \
	--disable-fortran --disable-doc --enable-mpi $(PFFT_CONFIGURE_FLAGS) --enable-openmp MPICC=$(CC) \
	2>&1 ; \
	make -j 4 2>&1 ; \
	make install 2>&1; \
	) | tee pfft-double.log | tail

$(PFFTFLIB): $(PFFT_CONFIGURE)
	mkdir -p depends/single;
	(cd depends/single; \
	$(PFFT_CONFIGURE) --prefix=$(abspath depends/install) --enable-single --disable-shared --enable-static  \
	--disable-fortran --disable-doc --enable-mpi $(PFFT_CONFIGURE_FLAGS_SINGLE) --enable-openmp MPICC=$(CC) \
	2>&1 ; \
       	make -j 4 2>&1 ; \
	make install 2>&1; \
	) | tee pfft-single.log | tail
	

$(PFFTSRC): 
	wget https://github.com/rainwoodman/pfft/releases/download/$(PFFT_VERSION)/pfft-$(PFFT_VERSION).tar.gz -O $@ ; \

