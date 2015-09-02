#
# qrpm
#   

# Define CC=cc
# Define OPENMP to enable MPI+OpenMP hybrid parallelization
# OPENMP  = -fopenmp # -openmp for Intel, -fopenmp for gcc

CC      = mpicc 
WOPT    ?= -Wall
CFLAGS  := -std=c99 -O3 -g $(WOPT)
LIBS    := -lm


# Compile options


# Define paths of FFTW3 & GSL libraries if necessary.

FFTW3_DIR ?= #e.g. /Users/jkoda/Research/opt/gcc/fftw3
GSL_DIR   ?= #e.g. /Users/jkoda/Research/opt/gcc/gsl

DIR_PATH = $(FFTW3_DIR) $(GSL_DIR) depends/install

CFLAGS += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LIBS   += $(foreach dir, $(DIR_PATH), -L$(dir)/lib) 
LIBS += -Llua
CFLAGS += -Ilua

EXEC = qrpm # halo
all: $(EXEC)

OBJS := main.o
OBJS += read_param_lua.o lpt.o msg.o power.o
OBJS += pm.o stepping.o 
OBJS += readrunpb.o
OBJS += domain.o
OBJS += timer.o 
OBJS += heap.o
OBJS += mond.o

LIBS += -ldl 
LIBS += -lgsl -lgslcblas
LIBS += -lfftw3f_mpi -lfftw3f
LIBS += -lm

lua/liblua.a: lua/Makefile
	(cd lua; CC="$(CC)" make generic)

ifdef OPENMP
  LIBS += -lfftw3f_omp
  #LIBS += -lfftw3f_threads       # for thread parallelization instead of omp
endif

qrpm: $(OBJS) lua/liblua.a
	$(CC) $(CFLAGS) $(OBJS) -llua $(LIBS) -o $@

.c.o:
	$(CC) $(CFLAGS) -c -o $@ $<
	
main.o: main.c parameters.h lpt.h particle.h msg.h power.h pm.h \
  stepping.h write.h timer.h mond.h version.h
mond.o : mond.h parameters.h msg.h pm.h
domain.o: domain.c domain.h msort.c permute.c
stepping.o: stepping.c particle.h msg.h stepping.h timer.h
comm.o: comm.c msg.h
heap.o: heap.c heap.h msg.h
kd_original.o: kd_original.c kd.h
lpt.o: lpt.c msg.h power.h particle.h
lpt_original.o: lpt_original.c
msg.o: msg.c msg.h
pm.o: pm.c pm.h particle.h msg.h timer.h
pm_debug.o: pm_debug.c pm.h particle.h msg.h timer.h
pm_original.o: pm_original.c stuff.h
power.o: power.c msg.h
read_param_lua.o: read_param_lua.c parameters.h msg.h
temp.o: temp.c msg.h particle.h
timer.o: timer.c msg.h timer.h
version.h:
	(echo "#define FPM_VERSION \\"; \
        git log -1 | sed -e 's;$$;\\n"\\;' -e 's;^;";'; \
        echo \"\") > version.h

.PHONY: clean run dependence version.h
clean :
	rm -f $(EXEC) $(OBJS) $(OBJS2) move_min.?

run:
	mpirun -n 2 ./qrpm param.lua

dependence:
	gcc -MM -MG *.c
