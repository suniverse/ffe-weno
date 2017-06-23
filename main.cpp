#include "Grid.hpp"
#include "Cow/src/HDF5.hpp"

#include <iostream>
#include <array>
#include <cmath>
#include <sstream>
#include <iomanip>

using namespace Cow;

int main()
{

	// Read Parameters
	UserParameters userParameters;
	//int* N = userParameters.GridLength;

	// Initialize Simulation Control
	SimulationControl simulationControl(userParameters);

	// Initialize SolutionData class: field of primitive variables
	SolutionData solution(userParameters);

	// Initialize Solver Class
	Solver solver(userParameters);

	// set the initial condition
	solution.InitialData();

	Array& P = solution.GetPField();

	solver.SetBoundaryValue(P);

	solution.ComputeTotalEnergy();
	solution.ComputeTotalEdotB();
	solution.ComputeMaxdivB();


	printf("Initial Total Energy %f \n",solution.TotalEnergy);
	printf("Initial Total EdotB %f \n",solution.TotalEdotB);
	printf("Initial Max divB %.12f \n",solution.MaxdivB);

	Array& B = solution.GetBField();
	Array& E = solution.GetEField();

	// Output the initial conditions
	simulationControl.OutputData(B,E);
	simulationControl.UpdateCheckPoint();

	while(simulationControl.ShouldContinue())
		{
			double dt = simulationControl.Getdt();
			printf("Cycle = %d, Time = %f, dt = %f \n", simulationControl.Cycle, simulationControl.SimTime, dt);

			solver.Setdt(dt);

			solver.Advance(P);
			
			simulationControl.UpdateTimeandCycle(dt);

			if(simulationControl.ShouldOutput())
			{
				printf("Output Data # %d \n", simulationControl.OutputCheckPoint);
				B = solution.GetBField();
				E = solution.GetEField();

				simulationControl.OutputData(B,E);
				simulationControl.UpdateCheckPoint();

				simulationControl.UpdateNextDumpTime();
				
				solution.ComputeTotalEnergy();
				solution.ComputeTotalEdotB();
				solution.ComputeMaxdivB();


				printf("Total Energy %f \n",solution.TotalEnergy);
				printf("Total EdotB %f \n",solution.TotalEdotB);
				printf("Max divB %.12f \n",solution.MaxdivB);

			}

		}

	B = solution.GetBField();
	E = solution.GetEField();

	simulationControl.OutputData(B,E);
	solution.ComputeTotalEnergy();
	solution.ComputeTotalEdotB();
	solution.ComputeMaxdivB();


	printf("Final Total Energy %f \n",solution.TotalEnergy);
	printf("Final Total EdotB %f \n",solution.TotalEdotB);
	printf("Final Max divB %.12f \n",solution.MaxdivB);


	return 0;
}