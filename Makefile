#
# File: Makefile
# Description: Makefile for BBFMM3D
#

CC = g++
LD = g++
CFLAGS  = -c -Wall -O3 -I ./include/ 
LDPATH = -L/usr/lib
LDFLAGS = -llapack -lblas -lrfftw -lfftw -lm 
PFLAG  =
SOURCES =  ./src/kernel_Types.cpp ./src/H2_3D_Tree.cpp ./src/read_metadata.cpp ./src/read_sources.cpp ./src/write_Into_Binary_File.cpp

SOURCEA = ./examples/get_input_through_routine_standard_kernel.cpp
SOURCEB = ./examples/binary_file_standard_kernel.cpp
SOURCEC = ./examples/get_input_through_routine_mykernel.cpp
SOURCED = ./examples/binary_file_mykernel.cpp


OBJECTA=$(SOURCES:.cpp=.o) $(SOURCEA:.cpp=.o)
OBJECTB=$(SOURCES:.cpp=.o) $(SOURCEB:.cpp=.o)
OBJECTC=$(SOURCES:.cpp=.o) $(SOURCEC:.cpp=.o)
OBJECTD=$(SOURCES:.cpp=.o) $(SOURCED:.cpp=.o)


EXECUTABLEA= ./exec/get_input_through_routine_standard_kernel
EXECUTABLEB=  ./exec/binary_file_standard_kernel
EXECUTABLEC=  ./exec/get_input_through_routine_mykernel
EXECUTABLED=  ./exec/binary_file_mykernel


get_input_through_routine_standard_kernel: $(SOURCES) $(SOURCEA) $(EXECUTABLEA)
$(EXECUTABLEA): $(OBJECTA)
	$(CC) $(OBJECTA) $(LDPATH) $(LDFLAGS)  -o $@

binary_file_standard_kernel: $(SOURCES) $(SOURCEB) $(EXECUTABLEB)
$(EXECUTABLEB): $(OBJECTB)
	$(CC)  $(OBJECTB) $(LDPATH) $(LDFLAGS)  -o $@

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
