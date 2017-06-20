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

	solution.ComputeTotalEnergy();
	solution.ComputeTotalEdotB();

	printf("total energy %f \n",solution.TotalEnergy);
	printf("total EdotB %f \n",solution.TotalEdotB);


	solver.SetBoundaryValue(P);

	Array& B = solution.GetBField();
	Array& E = solution.GetEField();

	// Output the initial conditions
	simulationControl.OutputData(B,E);
	simulationControl.UpdateCheckPoint();


	while(simulationControl.ShouldContinue())
		{
			double dt = simulationControl.Getdt();
			//printf("time %f %f %d\n", simulationControl.SimTime, dt, simulationControl.Cycle);

			solver.Setdt(dt);
			solver.Advance(P);
			//printf("output %d \n", simulationControl.ShouldOutput());
			
			simulationControl.UpdateTimeandCycle(dt);

			if(simulationControl.ShouldOutput())
			{
				//printf("output \n");
				B = solution.GetBField();
				E = solution.GetEField();

				simulationControl.OutputData(B,E);
				simulationControl.UpdateCheckPoint();

				simulationControl.UpdateNextDumpTime();
			}
		}

	B = solution.GetBField();
	E = solution.GetEField();

	simulationControl.OutputData(B,E);
	solution.ComputeTotalEnergy();
	solution.ComputeTotalEdotB();

	printf("total energy %f \n",solution.TotalEnergy);
	printf("total EdotB %f \n",solution.TotalEdotB);

	return 0;
}