# Makefile /dynamics/
CXX=g++
CXXFLAGS=-O2 -fopenmp -g
OBJ=nrutil.o calc_kernel.o mct_slit.o svd.o time_integration.o

all:  MCT_slit

MCT_slit: $(OBJ)
	$(CXX) $(CXXFLAGS) -o MCT_slit $(OBJ)


nrutil.o: nrutil.cpp
	$(CXX) $(CXXFLAGS) -c nrutil.cpp
svd.o: svd.cpp
	$(CXX) $(CXXFLAGS) -c svd.cpp
calc_kernel.o: calc_kernel.cpp
	$(CXX) $(CXXFLAGS) -c calc_kernel.cpp
time_integration.o: time_integration.cpp
	$(CXX) $(CXXFLAGS) -c time_integration.cpp
mct_slit.o: mct_slit.cpp
	$(CXX) $(CXXFLAGS) -c mct_slit.cpp

clean:
	rm *.o
