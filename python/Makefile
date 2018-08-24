.DEFAULT_GOAL := all

CC = g++
LD = g++

FFTW_ROOT = /usr/local
FFTW_LIB = $(FFTW_ROOT)/lib
FFTW_INCLUDE = $(FFTW_ROOT)/include

PYTHON_VERSION = 2.7
PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)

# location of the Boost Python include files and library
BOOST_INC = /usr/include
BOOST_LIB = /usr/lib

TARGET_TREE = FMMTree
TARGET_COMPUTE = FMMCompute
SRC_H2_3D_Tree = ../src/H2_3D_Tree
SRC_Tree = ../src/FMMTree
SRC_COMPUTE = ../src/FMMCompute

CFLAGS  = -Wall -O3 -I $(MKLROOT)/include -L $(MKLROOT)/lib/intel64 -I ../include -I/usr/include -I$(FFTW_INCLUDE)  \
-I/usr/include/python2.7 -I/usr/include/x86_64-linux-gnu/python2.7 -fopenmp -fPIC
LDPATH = -I $(MKLROOT)/include -L $(MKLROOT)/lib/intel64  -L/usr/lib -I/usr/include -I ../include \
 -L/usr/lib/python2.7/config-x86_64-linux-gnu 
LDFLAGS =  -L$(FFTW_LIB) -lfftw3 -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group \
-lpthread -lm -ldl -fopenmp -I/usr/include -lpython2.7 -lboost_python -lboost_system


$(TARGET_TREE).so: $(SRC_Tree).o $(SRC_H2_3D_Tree).o
	$(CC) -shared -Wl,--export-dynamic $(SRC_Tree).o  $(SRC_H2_3D_Tree).o -L$(BOOST_LIB)  \
	-l:libboost_python-py$(subst .,,$(PYTHON_VERSION)).so -L/usr/lib/python$(PYTHON_VERSION)/config \
	-lpython$(PYTHON_VERSION) $(LDPATH) $(LDFLAGS) -o $(TARGET_TREE).so $(CFLAGS)
$(SRC_Tree).o: $(SRC_Tree).cpp 
	$(CC) $(LDPATH) $(LDFLAGS) -I$(PYTHON_INCLUDE) -I$(BOOST_INC) -fPIC -c $(SRC_Tree).cpp $(CFLAGS) -o $@
$(SRC_H2_3D_Tree).o: $(SRC_H2_3D_Tree).cpp
	$(CC) $(LDPATH) $(LDFLAGS) -c $(SRC_H2_3D_Tree).cpp $(CFLAGS) -o $@
$(TARGET_COMPUTE).so: $(SRC_COMPUTE).o $(SRC_H2_3D_Tree).o
	$(CC) -shared -Wl,--export-dynamic $(SRC_COMPUTE).o $(SRC_H2_3D_Tree).o -L$(BOOST_LIB)  \
	-l:libboost_python-py$(subst .,,$(PYTHON_VERSION)).so -L/usr/lib/python$(PYTHON_VERSION)/config \
	-lpython$(PYTHON_VERSION) $(LDPATH) $(LDFLAGS) -I$(PYTHON_INCLUDE) -I$(BOOST_INC) -fPIC -o $(TARGET_COMPUTE).so $(CFLAGS)

$(SRC_COMPUTE).o: $(SRC_COMPUTE).cpp 
	$(CC) $(LDPATH) $(LDFLAGS) -I$(PYTHON_INCLUDE) -I$(BOOST_INC) -fPIC -c $(SRC_COMPUTE).cpp $(CFLAGS) -o $@
.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

all: $(TARGET_TREE).so $(TARGET_COMPUTE).so

clean:
	rm -rf *.o *~ ../src/*.o ../examples/*.o