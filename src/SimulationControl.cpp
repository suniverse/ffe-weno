#include "../Cow/src/Array.hpp"
#include "../Cow/src/HDF5.hpp"
#include "../Cow/src/MPI.hpp"

#include "Grid.hpp"
#include <cmath>
#include <sstream>
#include <iomanip>

using namespace Cow;

SimulationControl::SimulationControl(UserParameters userParameters)
{
	double NumCFL = userParameters.CFLSafetyNumber;

	// loop over three dimensions
	for (int d=0; d<3; ++d)
	{
		// set dx in that direction
		dx[d] = userParameters.dx[d];
		dtFixed = fmin(dx[d]*NumCFL, dtFixed);
	}

	NextDataDumpTime = userParameters.DataDumpInterval;
	DataDumpInterval = userParameters.DataDumpInterval;
	StopTime = userParameters.StopTime;
	StopCycle = userParameters.StopCycle;
}

double SimulationControl::Getdt() const
{
	double dt = dtFixed;
	if ((SimTime+dtFixed) > NextDataDumpTime)
	{
		dt = NextDataDumpTime-SimTime;
	}

	return dt;
}

bool SimulationControl::ShouldOutput() const
{
	return (SimTime >= NextDataDumpTime);
}

bool SimulationControl::ShouldContinue() const
{
	bool isAfterStopTime = (SimTime >= StopTime) && (StopTime > 0.0);
	bool isAfterStopCycle = (Cycle >= StopCycle) && (StopCycle > 0);

	return !(isAfterStopCycle || isAfterStopTime);
	
}

void SimulationControl::UpdateTimeandCycle(double dt)
{
	SimTime += dt;
	Cycle += 1;
}

void SimulationControl::UpdateNextDumpTime()
{
	NextDataDumpTime += DataDumpInterval;
	NextDataDumpTime = fmin(NextDataDumpTime, StopTime);
}

void SimulationControl::UpdateCheckPoint()
{
	OutputCheckPoint += 1;
}

void SimulationControl::OutputData(MpiCommunicator world, Cart cart, SolutionData solution) const
{
//	auto world = MpiCommunicator::world();
	int rank = world.rank();

	std::ostringstream stream;
	stream << "./chkpt" << std::setfill('0') << std::setw(4) << OutputCheckPoint <<".cpu"
		   << std::setfill('0') << std::setw(4)<< rank << ".h5";

	Cow::Array B = solution.GetBField();
	Cow::Array E = solution.GetEField();
	//Cow::Array P = solution.GetPField();

	auto F = H5::File(stream.str(),"w");
	F.writeArray("magnetic",B);
	F.writeArray("electric",E);

	
	if ( solution.SpongeLayer > 0 ) {
		Cow::Array DE = solution.GetDampedEnergy();
		F.writeArray("DampedEnergy", DE);
	}

	//F.writeArray("prim",P);
	auto G = F.createGroup("TimeSeriesData");

	// perdiodically:
	G.writeVectorDouble ("time", TimeSeries);
	G.writeVectorDouble ("TotalEnergy", TotalEnergy);
	G.writeVectorDouble ("ElectricEnergy", ElectricEnergy);
	G.writeVectorDouble ("MagneticEnergy", MagneticEnergy);


	auto G1 = F.createGroup("StartIndex");
	G1.writeInt("startI", cart.StartPos[0]);
	G1.writeInt("startJ", cart.StartPos[1]);
	G1.writeInt("startK", cart.StartPos[2]);

	auto G2 = F.createGroup("EndIndex");
	G2.writeInt("endI", cart.EndPos[0]);
	G2.writeInt("endJ", cart.EndPos[1]);
	G2.writeInt("endK", cart.EndPos[2]);

	auto G3 = F.createGroup("GlobalDim");
	G3.writeInt("GridDimI", cart.GlobalDim[0]);
	G3.writeInt("GridDimJ", cart.GlobalDim[1]);
	G3.writeInt("GridDimK", cart.GlobalDim[2]);
}

void SimulationControl::UpdateTimeSeries(MpiCommunicator world, SolutionData solution)
{
	// Compute Energy of solutions
	solution.ComputeEnergy();

	TimeSeries.push_back(SimTime);
	double Beng = world.dsum(solution.MagneticEnergy);
	double Eeng = world.dsum(solution.ElectricEnergy);
	double Teng = world.dsum(solution.TotalEnergy);

	MagneticEnergy.push_back(Beng);
	ElectricEnergy.push_back(Eeng);
	TotalEnergy.push_back(Teng);

}
