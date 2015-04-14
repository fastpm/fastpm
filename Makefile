#
# cola_halo
#   

# Define OPENMP to enable MPI+OpenMP hybrid parallelization
# OPENMP  = -fopenmp # -openmp for Intel, -fopenmp for gcc

CC      = mpicc -std=c99 
WOPT    ?= -Wall
CFLAGS  := -O3 -g $(WOPT) $(OPENMP) -Wall
LIBS    := -lm


# Compile options


# Define paths of FFTW3 & GSL libraries if necessary.

LUA_DIR   ?= /opt/local
FFTW3_DIR ?= #e.g. /Users/jkoda/Research/opt/gcc/fftw3
GSL_DIR   ?= #e.g. /Users/jkoda/Research/opt/gcc/gsl

DIR_PATH = $(FFTW3_DIR) $(GSL_DIR) $(LUA_DIR)

CFLAGS += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LIBS   += $(foreach dir, $(DIR_PATH), -L$(dir)/lib)

EXEC = cola_halo # halo
all: $(EXEC)

OBJS := main.o
OBJS += read_param_lua.o lpt.o msg.o power.o
OBJS += pm.o cola.o fof.o comm.o #move.o move_min.o
OBJS += readrunpb.o
OBJS += domain.o
OBJS += write.o timer.o mem.o
OBJS += subsample.o coarse_grid.o

LIBS += -llua -ldl 
LIBS += -lgsl -lgslcblas
LIBS += -lfftw3f_mpi -lfftw3f


ifdef OPENMP
  LIBS += -lfftw3f_omp
  #LIBS += -lfftw3f_threads       # for thread parallelization instead of omp
endif

cola_halo: $(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $@

main.o: main.c parameters.h lpt.h particle.h msg.h power.h comm.h pm.h \
  cola.h fof.h write.h timer.h mem.h subsample.h coarse_grid.h
cola.o: cola.c particle.h msg.h cola.h timer.h
cola_original.o: cola_original.c
comm.o: comm.c msg.h comm.h
fof.o: fof.c particle.h msg.h comm.h timer.h
kd_original.o: kd_original.c kd.h
lpt.o: lpt.c msg.h power.h particle.h
lpt_original.o: lpt_original.c
mem.o: mem.c fof.h particle.h mem.h msg.h comm.h
msg.o: msg.c msg.h
pm.o: pm.c pm.h particle.h msg.h comm.h timer.h
pm_debug.o: pm_debug.c pm.h particle.h msg.h timer.h
pm_original.o: pm_original.c stuff.h
power.o: power.c msg.h
read_param_lua.o: read_param_lua.c parameters.h msg.h
temp.o: temp.c msg.h particle.h comm.h
timer.o: timer.c msg.h timer.h
write.o: write.c msg.h comm.h write.h particle.h


#
# "halo" -- cola_halo without cola, only does FoF etc.
#
OBJS2 := halo_main.o read.o
OBJS2 += read_param_lua.o msg.o fof.o comm.o
OBJS2 += write.o timer.o mem.o 
OBJS2 += subsample.o coarse_grid.o

halo: $(OBJS2)
	$(CC) $(OBJS2) $(LIBS) -o $@



.PHONY: clean run dependence
clean :
	rm -f $(EXEC) $(OBJS) $(OBJS2) move_min.?

run:
	mpirun -n 2 ./cola_halo param.lua

dependence:
	gcc -MM -MG *.c
