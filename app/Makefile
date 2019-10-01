
FFTW_ROOT = /opt/apps/intel18/impi18_0/fftw3/3.3.8
FFTW_LIB = $(FFTW_ROOT)/lib
FFTW_INCLUDE = $(FFTW_ROOT)/include

CFLAGS = -Wall -g -c -isystem $(EIGEN_DIR) -I $(MKLROOT)/include -I ../include/ -I$(FFTW_INCLUDE) -fopenmp
LDPATH = -L../include -L $(MKLROOT)/include -L $(MKLROOT)/lib/intel64
LDFLAGS =  -L$(FFTW_LIB) -lfftw3 -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group \
-lpthread -lm -ldl -fopenmp

SRC =  gmres.cpp ../src/kernel_Types.cpp ../src/H2_3D_Tree.cpp ../src/read_metadata.cpp ../src/read_sources.cpp ../src/write_Into_Binary_File.cpp

OBJ=$(SRC:.cpp=.o)


all: $(OBJ)
	g++ $(OBJ) $(LDPATH) $(LDFLAGS) -o gmres

.cpp.o:
	g++ $(CFLAGS) $< -o $@

clean:
	rm -rf $(OBJ)

