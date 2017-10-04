#include "Grid.hpp"
#include "../Cow/src/HDF5.hpp"
#include "../Cow/src/MPI.hpp"

#include <iostream>
#include <array>
#include <cmath>
#include <sstream>
#include <iomanip>

using namespace Cow;

int main()
{
	// initialize MPI
	MpiSession Mpi;
	auto world = MpiCommunicator::world();
	auto rank = world.rank();

	// Read Parameters
	UserParameters userParameters;

	// initialize MPI
	Cart mycart(userParameters);

    userParameters.SetPatchDimension(mycart);

	// Initialize Simulation Control
	SimulationControl simulationControl(userParameters);

	// Initialize SolutionData class: field of primitive variables
	SolutionData solution(userParameters, mycart);

	// Initialize Solver Class
	Solver solver(userParameters);

	// set the initial condition
	solution.InitialData();

	// get the array for prim field
	Array& P = solution.GetPField();
	// get the array for resistivity
	Array& resistivity = solution.GetResistivity();

	// set boundary for each patches
	solver.SetBoundaryValue(mycart, P);
	solver.SetBoundaryValue(mycart, resistivity);


	// compute and print parameters from initial data
	solution.ComputeEnergy();
	solution.ComputeMaxEdotB();
	solution.ComputeMaxdivB();

	double MaxdivB = world.maximum(solution.MaxdivB);
	double MaxEdotB = world.maximum(solution.MaxEdotB);
	double TotalEnergy = world.dsum(solution.TotalEnergy);
	double OhmHeat = world.dsum(solution.TotalOhmHeat);


	if(0==rank)
	{
		printf("INITIAL DATA: Total Energy %.12e \n", TotalEnergy);
		printf("INITIAL DATA: Max EdotB %.12e \n", MaxEdotB);
		if(MaxEdotB > 1e-14)
		{
			printf("WARNING !!! INITIAL DATA NOT SATISFY FORCE FREE CONDITION !!!\n");
		}
		printf("INITIAL DATA: Max divB %.12e \n", MaxdivB);
		if(MaxdivB > 1e-2)
		{
			printf("WARNING !!! INITIAL DATA HAS NON-ZERO DIV B !!!\n");
		}
	}

	// Output the initial conditions
	simulationControl.OutputData(world, mycart, solution);
	simulationControl.UpdateTimeSeries(world, solution);
	simulationControl.UpdateCheckPoint();

	// INITIALIZATION FINISHED

	// main loop for evolution

	while(simulationControl.ShouldContinue())
	{
		double dt = simulationControl.Getdt();

		if(rank==0){
			printf("Cycle = %d, Time = %f, dt = %f \n", simulationControl.Cycle, simulationControl.SimTime, dt);
		}
		solver.Setdt(dt);
		solution.ComputeDampedEnergy(dt);

		solver.Advance(mycart, P, resistivity);

		simulationControl.UpdateTimeandCycle(dt);
		simulationControl.UpdateTimeSeries(world, solution);

		// Do some output at check point
		if(simulationControl.ShouldOutput())
		{
			if(rank==0){
				printf("Output Data # %d \n", simulationControl.OutputCheckPoint);
			}

			// output some info to cmd line	
			solution.ComputeEnergy();
			solution.ComputeMaxEdotB();
			solution.ComputeMaxdivB();

			MaxdivB = world.maximum(solution.MaxdivB);
			MaxEdotB = world.maximum(solution.MaxEdotB);
			TotalEnergy = world.dsum(solution.TotalEnergy);
		    OhmHeat = world.dsum(solution.TotalOhmHeat);


			if(0==rank){
				printf("Total Energy %.12e \n", TotalEnergy);
				printf("Max EdotB %.12e \n", MaxEdotB);
				printf("Max divB %.12e \n", MaxdivB);
				printf("Ohm Heat %.12e \n", OhmHeat);

			}

			// output to file
			simulationControl.OutputData(world, mycart, solution);
			// update control info
			simulationControl.UpdateCheckPoint();
			simulationControl.UpdateNextDumpTime();
		}// finish checkpoint
	}// finish main evolution loop

	// Output Final Data

	simulationControl.OutputData(world, mycart, solution);
	solution.ComputeEnergy();
	solution.ComputeMaxEdotB();
	solution.ComputeMaxdivB();

	MaxdivB = world.maximum(solution.MaxdivB);
	MaxEdotB = world.maximum(solution.MaxEdotB);
	/*TotalEnergy = world.dsum(solution.TotalEnergy);
			world.inSequence([&] (int rank)
			{
				printf("Final Patch Total Energy %f \n", solution.TotalEnergy);
			});*/


	if(0==rank)
	{
		printf("Final Total Energy %.12e \n", TotalEnergy);
		printf("Final Max EdotB %.12e \n", MaxEdotB);
		printf("Final Max divB %.12e \n", MaxdivB);	
	}



	return 0;
}