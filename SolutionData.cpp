#include "Grid.hpp"
#include <cmath>
using namespace Cow;


SolutionData::SolutionData (UserParameters userParameters)
{
	//BField = Array(Nx, Ny, Nz, 3);
	//EField = Array(Nx, Ny, Nz, 3);
	int NumGhost = userParameters.NumberOfGhostZones;

	// loop over three dimensions
	for (int d=0; d<3; ++d)
	{
		// set number of grids in that dimension
		N[d] = userParameters.GridLength[d];

		if (N[d]==1) continue;

		// dimension of interior (no ghost zone)
		Nint[d] = N[d]-2*NumGhost;

		// set start and end index in that dimension (exclude ghost zones)
		StartIndex[d] = NumGhost;
		EndIndex[d] = N[d] - NumGhost;
		dx[d] = 1./Nint[d];
	}

	// coordinates of each grid point
	x = Array(Nint[0], Nint[1], Nint[2], 3);

	for (int i = 0; i < Nint[0]; ++i)
	{
		for (int j = 0; j < Nint[1]; ++j)
		{
			for (int k = 0; k < Nint[2]; ++k)
			{
				x(i,j,k,0) = (Nint[0]==1)? 0 : i*dx[0]+dx[0]/2.;
				x(i,j,k,1) = (Nint[1]==1)? 0 : j*dx[1]+dx[1]/2.;
				x(i,j,k,2) = (Nint[2]==1)? 0 : k*dx[2]+dx[2]/2.;
			} // end loop over k
		}// end loop over j
	}// end loop over i

	PField = Array(N[0], N[1], N[2], 6);
	Region BRegion, ERegion;

	for (int d=0; d<3; ++d)
	{
		ERegion.lower[d] = StartIndex[d];
		ERegion.upper[d] = EndIndex[d];
		BRegion.lower[d] = StartIndex[d];
		BRegion.upper[d] = EndIndex[d];
	}

	BRegion.lower[3]  = 0;
	BRegion.upper[3]  = 3;
	ERegion.lower[3]  = 3;
	ERegion.upper[3]  = 6;

	BField = PField[BRegion];
	EField = PField[ERegion];
}

Array& SolutionData::GetBField()
{
	Region BRegion;

	for (int d=0; d<3; ++d)
	{
		BRegion.lower[d] = StartIndex[d];
		BRegion.upper[d] = EndIndex[d];
	}

	BRegion.lower[3]  = 0;
	BRegion.upper[3]  = 3;

	BField = PField[BRegion];

	return BField;
}

Array& SolutionData::GetEField()
{
	Region ERegion;

	for (int d=0; d<3; ++d)
	{
		ERegion.lower[d] = StartIndex[d];
		ERegion.upper[d] = EndIndex[d];
	}

	ERegion.lower[3]  = 3;
	ERegion.upper[3]  = 6;

	EField = PField[ERegion];

	return EField;
}

Array& SolutionData::GetPField()
{
	return PField;
}

void SolutionData::CombineBEtoP()
{
	Region BRegion, ERegion;

	for (int d=0; d<3; d++)
	{
		BRegion.lower[d] = StartIndex[d];
		BRegion.upper[d] = EndIndex[d];
		ERegion.lower[d] = StartIndex[d];
		ERegion.upper[d] = EndIndex[d];
	}

	BRegion.lower[3]  = 0;
	BRegion.upper[3]  = 3;
	ERegion.lower[3]  = 3;
	ERegion.upper[3]  = 6;

	PField[BRegion] = BField;
	PField[ERegion] = EField;
}

void SolutionData::ComputeTotalEnergy()
{
	TotalEnergy = 0;
	for (int i = StartIndex[0]; i < EndIndex[0]; ++i)
	{
		for (int j = StartIndex[1]; j < EndIndex[1]; ++j)
		{
			for (int k = StartIndex[2]; k < EndIndex[2]; ++k)
			{
				for (int d=0; d<6; ++d) 
				{
					TotalEnergy  += PField(i,j,k,d)*PField(i,j,k,d);
				}

			} // end loop over k
		}// end loop over j
	}// end loop over i
	TotalEnergy = TotalEnergy/2.;
}

void SolutionData::ComputeTotalEdotB()
{
	TotalEdotB = 0;
	for (int i = StartIndex[0]; i < EndIndex[0]; ++i)
	{
		for (int j = StartIndex[1]; j < EndIndex[1]; ++j)
		{
			for (int k = StartIndex[2]; k < EndIndex[2]; ++k)
			{
				for (int d=0; d<3; ++d) 
				{
					TotalEdotB  += PField(i,j,k,d)*PField(i,j,k,d+3);
				}

			} // end loop over k
		}// end loop over j
	}// end loop over i
}

void SolutionData::ComputeMaxdivB()
{
	MaxdivB = 0;
	for (int i = StartIndex[0]; i < EndIndex[0]; ++i)
	{
		for (int j = StartIndex[1]; j < EndIndex[1]; ++j)
		{
			for (int k = StartIndex[2]; k < EndIndex[2]; ++k)
			{
				double divB = 0;
				if ( N[0] > 1 )
				{
					divB += ((PField(i-2,j,k,0)-PField(i+2,j,k,0))/12. - (PField(i-1,j,k,0)-PField(i+1,j,k,0))*2./3.)/dx[0];
				}	
				if ( N[1] > 1 )
				{
					divB += ((PField(i,j-2,k,1)-PField(i,j+2,k,1))/12. - (PField(i,j-1,k,1)-PField(i,j+1,k,1))*2./3.)/dx[1];
				}
				if ( N[2] > 1 )
				{
					divB += ((PField(i,j,k-2,2)-PField(i,j,k+2,2))/12. - (PField(i,j,k-1,2)-PField(i,j,k+1,2))*2./3.)/dx[2];
				}	
				MaxdivB = fmax(divB, MaxdivB);
			} // end loop over k
		}// end loop over j
	}// end loop over i
}


void SolutionData::InitialData()
{
	double kvec[3] = {0, 2*3.1415927, 2*3.1415927};

	for (int i = 0; i < Nint[0]; ++i)
	{
		for (int j = 0; j < Nint[1]; ++j)
		{
			for (int k = 0; k < Nint[2]; ++k)
			{	
				double kx = 0;

				for (int d=0; d<3; ++d)
				{
					kx += kvec[d]*x(i,j,k,d);
				}

				//Stationary Alfve Wave
				/*BField(i,j,k,0) = 1.;
				BField(i,j,k,1) = 1.;

				if ((x(i,j,k,0)-0.5)<=-0.1) BField(i,j,k,2) = 1;
				if ((x(i,j,k,0)-0.5)>-0.1 && (x(i,j,k,0)-0.5)<0.1) BField(i,j,k,2) = 1+0.15*(1+sin(5*3.1415927*(x(i,j,k,0)-0.5)));
				if ((x(i,j,k,0)-0.5)>0.1) BField(i,j,k,2) = 1.3;

				EField(i,j,k,2) = 1;
				EField(i,j,k,0) = -BField(i,j,k,2);*/

				// Fast Wave
				/*BField(i,j,k,0) = 0.1*cos(kx);
				BField(i,j,k,2) = 1;
				EField(i,j,k,1) = -BField(i,j,k,0);*/

				// Alfven Wave
				BField(i,j,k,0) = 0.1*cos(kx);
				BField(i,j,k,2) = 1;
				EField(i,j,k,1) = -BField(i,j,k,0);

				//Current Sheet
				/*BField(i,j,k,0) = 1;
				if ((x(i,j,k,0)-0.5)<0)
				{
					BField(i,j,k,1) = 2.;
				} 
				if ((x(i,j,k,0)-0.5)>0)
				{
					BField(i,j,k,1) = -2.;
				} */

				//Three Wave
				/*if ((x(i,j,k,0)-0.5)<0)
				{
					BField(i,j,k,0) = 1.0;
					BField(i,j,k,1) = 1.5;
					BField(i,j,k,2) = 3.5;
					EField(i,j,k,0) = -1.0;
					EField(i,j,k,1) = -0.5;
					EField(i,j,k,2) = 0.5;

				} 
				if ((x(i,j,k,0)-0.5)>0)
				{
					BField(i,j,k,0) = 1.0;
					BField(i,j,k,1) = 2.0;
					BField(i,j,k,2) = 2.2;
					EField(i,j,k,0) = -1.5;
					EField(i,j,k,1) = 1.3;
					EField(i,j,k,2) = -0.5;
				} */

				// Non-degenerate Alfven wave
				/*double beta = -0.5;
				double gamma = 1./sqrt(1-beta*beta);
				//field in the wave frame
				double Bx, By, Bz, Ex, Ey, Ez;

				Bx = 1.;
				By = 1.;

				if ((x(i,j,k,0)-0.5)<=-0.1) Bz = 1;
				if ((x(i,j,k,0)-0.5)>-0.1 && (x(i,j,k,0)-0.5)<0.1) Bz = 1+1.5*(x(i,j,k,0)-0.5+0.1);
				if ((x(i,j,k,0)-0.5)>0.1) Bz = 1.3;

				Ex = 0;
				Ey = 0;
				Ez = 1.;

				BField(i,j,k,0) = Bx;
				BField(i,j,k,1) = gamma*(By - beta*Ez);
				BField(i,j,k,2) = gamma*Bz;

				EField(i,j,k,0) = -Bz;
				EField(i,j,k,1) = beta*BField(i,j,k,2);
				EField(i,j,k,2) = BField(i,j,k,1);*/

				// Degenerate Alfven wave
				/*double beta = -0.5;
				double gamma = 1./sqrt(1-beta*beta);
				//field in the wave frame
				double Bx, By, Bz, Ex, Ey, Ez;
				double phi;

				Bx = 0.;

				if ((x(i,j,k,0)-0.5)<=-0.1) phi = 0;
				if ((x(i,j,k,0)-0.5)>-0.1 && (x(i,j,k,0)-0.5)<0.1) phi = 2.5*3.1415927*(x(i,j,k,0)-0.5+0.1);
				if ((x(i,j,k,0)-0.5)>0.1) Bz = 3.1415927/2.;

				By = 2*cos(phi);
				Bz = 2*sin(phi);

				Ex = 0;
				Ey = 0;
				Ez = 0.;

				BField(i,j,k,0) = Bx;
				BField(i,j,k,1) = gamma*By;
				BField(i,j,k,2) = gamma*Bz;

				EField(i,j,k,0) = 0;
				EField(i,j,k,1) = -beta*BField(i,j,k,2);
				EField(i,j,k,2) = beta*BField(i,j,k,1);*/

			} // end loop over k
		}// end loop over j
	}// end loop over i
	CombineBEtoP();
}