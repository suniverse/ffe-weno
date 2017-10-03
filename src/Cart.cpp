#include "Grid.hpp"
#include "Roe.hpp"
#include "Reconstruction.hpp"
#include <cmath>
#include <array>
#include "../Cow/src/MPI.hpp"

using namespace Cow;

Cart::Cart(UserParameters userParameters)
{
	SetGlobalDim(userParameters);
	// Create Processor Lattice
	auto world = MpiCommunicator::world();
	// Create cart, split in directions where GlobalDim>1
	cart = world.createCartesian (3,  {GlobalDim[0] > 1, GlobalDim[1] > 1, GlobalDim[2] > 1});

	// Get lattice Dimension
	LatticeDim = cart.getNumberOfDimensions();
	// Get 3D size of the lattice
    auto pdim = cart.getDimensions();
    // Get the coordiates of the block patch
    auto coord = cart.getCoordinates();

	for (int d=0; d<3; ++d)
	{
		// Compute the Dimension of Small Patches
		// Distribute residue among first several processors
		int res = GlobalDim[d]%pdim[d];
		int ind = ( res > coord[d]);
		PatchSize[d] = GlobalDim[d]/pdim[d] + ind ;

		Coordinates[d] = coord[d];
		StartPos[d] = GlobalDim[d]/pdim[d]*coord[d] + fmin(coord[d], res);
		EndPos[d] = StartPos[d] + PatchSize[d];

		if ((GlobalDim[d]!= 1) && (PatchSize[d] < NumGhost))
		{
			throw std::runtime_error("MPI PatchSize ERROR: Interior Region is smaller than Ghost Zones! \n");
		}

	}

}

// Read global size of box from UserParameters
void Cart::SetGlobalDim(UserParameters userParameters)
{
	NumGhost = userParameters.NumberOfGhostZones;
	for (int d=0; d<3; ++d)
	{
		GlobalDim[d] = userParameters.GridDimensions[d];
	}
}

// Communicate to apply periodic boundary condition
void Cart::apply (Array &A, int NumGhost, int N[3])
{
    for (int d=0; d<3; ++d){

    	if (N[d] == 1) continue;

		Region guardL, guardR, interiorL, interiorR;

		guardL.lower[d] = 0;
		guardL.upper[d] = NumGhost;

		interiorL.lower[d] = NumGhost;
		interiorL.upper[d] = NumGhost * 2;

		interiorR.lower[d] = -2*NumGhost;
		interiorR.upper[d] = -NumGhost;

		guardR.lower[d] = -NumGhost;
		guardR.upper[d] = 0;

		cart.shiftExchange (A, d, 'L', interiorL, guardR);
    	cart.shiftExchange (A, d, 'R', interiorR, guardL);

    }
}