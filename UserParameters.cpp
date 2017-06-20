#include "Grid.hpp"
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <sstream>

using namespace Cow;

UserParameters::UserParameters()
{
	ReadParameters();

	ErrorCheck();

	for (int n = 0; n < 3; ++n)
	{
		if (GridDimensions[n] > 1)
		{
			GridLength[n] = GridDimensions[n] + 2 * NumberOfGhostZones;
		}
	}

}

void UserParameters::ReadParameters()
{
	FILE * parameterfile;
	int MAX_LINE_LENGTH = 128;
	char line[MAX_LINE_LENGTH];

	parameterfile = fopen("parameters.dat", "r");
	//rewind(parameterfile);

	while (fgets(line, MAX_LINE_LENGTH, parameterfile) != NULL) 
	{
		int ret = 0;
		ret += sscanf(line, "NumberOfGhostZones = %d", &NumberOfGhostZones);
		ret += sscanf(line, "GridDimensions = %d %d %d", GridDimensions, GridDimensions+1, GridDimensions+2);
		ret += sscanf(line, "CFLSafetyNumber = %lf", &CFLSafetyNumber);
		ret += sscanf(line, "DataDumpInterval = %lf", &DataDumpInterval);
		ret += sscanf(line, "StopTime = %lf", &StopTime);
		ret += sscanf(line, "StopCycle = %d", &StopCycle);
		ret += sscanf(line, "ReconstructionMethod = %d", &ReconstructionMethod);

	}
}

void UserParameters::ErrorCheck()
{
	// Error check: GridDimensions must be larger than Ghost Zones
	for (int d=0; d<3; ++d)
	{
		if (GridDimensions[d] == 1) continue;
		if (GridDimensions[d] <= NumberOfGhostZones)
		{
			std::ostringstream stream;
			stream << "PARAMETER ERROR: Too Few Grids in Dimension [";
			stream << d << "]!\n";
			throw std::runtime_error(stream.str());
		}
	}

	// Error check: StopTime and StopCycle can't both be negative
	if ((StopCycle < 0) && (StopTime < 0.))
	{
		std::ostringstream stream;
		stream << "PARAMETER ERROR: Either StopTime or StopCycle should be assigned! \n";
		throw std::runtime_error(stream.str());
	}
}