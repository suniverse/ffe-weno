-include Makefile.in

COW = Cow/src/libcow.a
SRC = $(filter-out src/main.cpp, $(wildcard src/*.cpp))
OBJ = $(SRC:%.cpp=%.o)

default : ffe-weno

$(COW) :
	$(MAKE) -C Cow

ffe-weno: src/main.cpp $(SRC) $(COW)
	$(CXX) $(CFLAGS) -o $@ $^ $(COW) $(H5I) $(H5L)

#ffe-weno : main.cpp SolutionData.cpp Solver.cpp UserParameters.cpp SimulationControl.cpp Roe.cpp Reconstruction.cpp Cart.cpp $(COW)
#	$(CXX) $(CFLAGS) -o $@ $^ $(COW) $(H5L)

clean :
	$(MAKE) -C Cow clean
	$(RM) $(OBJ) src/*.o ffe-weno

