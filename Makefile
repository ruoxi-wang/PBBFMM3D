CC = g++
LD = g++

FFTW_ROOT = /opt/apps/intel18/impi18_0/fftw3/3.3.8
FFTW_LIB = $(FFTW_ROOT)/lib
FFTW_INCLUDE = $(FFTW_ROOT)/include

# CFLAGS  = -c -Wall -O3 -pg -I ./include/ -I/usr/include -I$(FFTW_INCLUDE) -fopenmp
# LDPATH = -L/usr/lib -I/usr/include -I ./include/ 
# LDFLAGS = -pg -O3 -llapack -lblas -L$(FFTW_LIB) -lfftw3 -lm -fopenmp
CFLAGS  = -c -Wall -O3 -I $(MKLROOT)/include -L $(MKLROOT)/lib/intel64 -I ./include/ -I/usr/include -I$(FFTW_INCLUDE) -fopenmp
LDPATH = -I $(MKLROOT)/include -L $(MKLROOT)/lib/intel64 -L/ -L/usr/lib -I/usr/include -I ./include/
LDFLAGS =  -L$(FFTW_LIB) -lfftw3 -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group \
-lpthread -lm -ldl -fopenmp
MFLAGS = -lm
PFLAG  =
SOURCES =  ./src/kernel_Types.cpp ./src/H2_3D_Tree.cpp ./src/read_metadata.cpp ./src/read_sources.cpp ./src/write_Into_Binary_File.cpp

SOURCEA = ./examples/get_input_through_routine_standard_kernel.cpp
SOURCEB = ./examples/binary_file_standard_kernel.cpp
SOURCEC = ./examples/get_input_through_routine_mykernel.cpp
SOURCED = ./examples/binary_file_mykernel.cpp
SOURCEHT = ./examples/3d_exp_cov.cpp

OBJECTA=$(SOURCES:.cpp=.o) $(SOURCEA:.cpp=.o)
OBJECTB=$(SOURCES:.cpp=.o) $(SOURCEB:.cpp=.o)
OBJECTC=$(SOURCES:.cpp=.o) $(SOURCEC:.cpp=.o)
OBJECTD=$(SOURCES:.cpp=.o) $(SOURCED:.cpp=.o)
OBJECTHT=$(SOURCES:.cpp=.o) $(SOURCEHT:.cpp=.o)


EXECUTABLEA= ./exec/get_input_through_routine_standard_kernel
EXECUTABLEB=  ./exec/binary_file_standard_kernel
EXECUTABLEC=  ./exec/get_input_through_routine_mykernel
EXECUTABLED=  ./exec/binary_file_mykernel
EXECUTABLEHT=  ./exec/3d_exp_cov

binary_file_standard_kernel: $(SOURCES) $(SOURCEB) $(EXECUTABLEB)
$(EXECUTABLEB): $(OBJECTB)
	$(CC)  $(OBJECTB) $(LDPATH) $(LDFLAGS)  -o $@

tomography: $(SOURCES) $(SOURCEHT) $(EXECUTABLEHT) $(EXECUTABLEHT)
$(EXECUTABLEHT): $(OBJECTHT)
	$(CC)  $(OBJECTHT) $(LDPATH) $(LDFLAGS)  -o $@

get_input_through_routine_standard_kernel: $(SOURCES) $(SOURCEA) $(EXECUTABLEA)
$(EXECUTABLEA): $(OBJECTA)
	$(CC) $(OBJECTA) $(LDPATH) $(LDFLAGS)  -o $@

get_input_through_routine_mykernel: $(SOURCES) $(SOURCEC) $(EXECUTABLEC)
$(EXECUTABLEC): $(OBJECTC)
	$(CC)  $(OBJECTC) $(LDPATH) $(LDFLAGS)  -o $@

binary_file_mykernel: $(SOURCES) $(SOURCED) $(EXECUTABLED)
$(EXECUTABLED): $(OBJECTD)
	$(CC)  $(OBJECTD) $(LDPATH) $(LDFLAGS)   -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@


clean:
	rm -rf *.o *~ ./src/*.o ./examples/*.o ./exec/*

tar:
	tar -zcvf BBFMM3D.tar.gz ./exec ./src ./include ./examples ./Makefile ./input ./output ./README.md
