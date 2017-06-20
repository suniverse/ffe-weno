-include Makefile.in

COW = Cow/src/libcow.a

default : ffe-weno

$(COW) :
	$(MAKE) -C Cow

ffe-weno : main.cpp SolutionData.cpp Solver.cpp UserParameters.cpp SimulationControl.cpp Roe.cpp Reconstruction.cpp $(COW)
	$(CXX) $(CFLAGS) -o $@ $^ $(COW) $(H5L)

clean :
	$(MAKE) -C Cow clean
	$(RM) *.o ffe-weno

