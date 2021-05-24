CC=pgcc 
EXE=exe
HEADERFILE=constants.h shortcut_wrappers.h cal_kvec.h global_variables.h io.h pfr_cufft1.h time_evolution.h debug_functions.h initialization.h misc.h curand_cuda.h inhom_v3.h
OBJ=initialization.o cal_kvec.o global_variables.o initial_conditions.o io.o pfr_cufft1.o time_evolution.o debug_functions.o misc.o  curand_cuda.o inhom_v3.o

#LDFLAGS = -lfftw3_omp -lfftw3 -lm 
#CFLAGS  = -fopenmp 
CFLAGS  = -fast -acc -ta=tesla:cuda10.2 -Minfo=accel -Mcudalib=cufft,curand -Minform=warn -mcmodel=medium 

PROJECT_NAME=Void_model
filename=main.o

%.o: %.c  $(HEADERFILE) 
	$(CC) -c $(CFLAGS) -o $@ $< $(LDFLAGS) 

all: ${PROJECT_NAME}

${PROJECT_NAME}: ${filename} $(OBJ) $(HEADERFILE)
	$(CC)  $(CFLAGS) -o $@.exe ${filename} $(OBJ) $(LDFLAGS)

curand_test: curand_test.c curand_cuda.o
	$(CC) $(CFLAGS) -o $@.exe $< curand_cuda.o $(LDFLAGS)

cufft_test: cufft_test.o pfr_cufft1.o
	$(CC) $(CFLAGS) -o $@.exe $< pfr_cufft1.o

test: test.c
	$(CC) $(CFLAGS) -o $@.exe $< 
clean:
	$(RM) -rf *.o HEAmodel.exe a.out 

cleanall:
	$(RM) -rf *.exe *.dat *.o HEAmodel.exe a.out *~ output_files/*.bin output_files/*.vtk output_files/*.dat

cleandata:
	$(RM) -rf output_files/*.bin output_files/*.vtk output_files/*.dat
