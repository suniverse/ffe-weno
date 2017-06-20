#include "Cow/src/Array.hpp"
#include "Cow/src/HDF5.hpp"
#include "Grid.hpp"
#include <cmath>
#include <sstream>
#include <iomanip>

using namespace Cow;

SimulationControl::SimulationControl(UserParameters userParameters)
{
	int NumGhost = userParameters.NumberOfGhostZones;
	double NumCFL = userParameters.CFLSafetyNumber;

	// loop over three dimensions
	for (int d=0; d<3; ++d)
	{
		if (userParameters.GridLength[d]==1) continue;
		// set dx in that direction
		dx[d] = 1./(userParameters.GridLength[d]-2*NumGhost);
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

void SimulationControl::OutputData(Array &B, Array &E) const
{
	std::ostringstream stream;
	stream << "chkpt" << std::setfill('0') << std::setw(4) << OutputCheckPoint << ".h5";
	auto F = H5::File(stream.str(),"w");
	F.writeArray("magnetic",B);
	F.writeArray("electric",E);
}
