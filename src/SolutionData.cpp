#include "Grid.hpp"
#include <cmath>
using namespace Cow;
#include <cstdlib>


SolutionData::SolutionData (UserParameters userParameters, Cart cart)
{
	//BField = Array(Nx, Ny, Nz, 3);
	//EField = Array(Nx, Ny, Nz, 3);
	int NumGhost = userParameters.NumberOfGhostZones;

	auto StartPos = cart.StartPos;
	auto GlobalDim = cart.GlobalDim;

	// loop over three dimensions
	for (int d=0; d<3; ++d)
	{
		// set number of grids in that dimension
		N[d] = userParameters.GridLength[d];

		if (N[d]==1) continue;
		dx[d] = userParameters.dx[d];

		// Calculate Grid Dimension
		GridDim += 1;

		// dimension of interior (no ghost zone)
		Nint[d] = N[d]-2*NumGhost;

		// set start and end index in that dimension (exclude ghost zones)
		StartIndex[d] = NumGhost;
		EndIndex[d] = N[d] - NumGhost;
	}

	// On/off sponge layer
	SpongeLayer = userParameters.SpongeLayer;

	AmplitudeOfAlfvenPacket = userParameters.AmplitudeOfAlfvenPacket;
	ntheta = userParameters.ntheta;

	// coordinates of each grid point
	x = Array(Nint[0], Nint[1], Nint[2], 3);

	for (int i = 0; i < Nint[0]; ++i)
	{
		for (int j = 0; j < Nint[1]; ++j)
		{
			for (int k = 0; k < Nint[2]; ++k)
			{
				x(i,j,k,0) = (Nint[0]==1)? 0 : i*dx[0]+dx[0]/2.+double(StartPos[0])/GlobalDim[0];
				x(i,j,k,1) = (Nint[1]==1)? 0 : j*dx[1]+dx[1]/2.+double(StartPos[1])/GlobalDim[1];
				x(i,j,k,2) = (Nint[2]==1)? 0 : k*dx[2]+dx[2]/2.+double(StartPos[2])/GlobalDim[2];
			} // end loop over k
		}// end loop over j
	}// end loop over i

	// Initialize EM field values
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

	// Initialize Resistivity
	//if ( SpongeLayer > 0 ) {
	Resistivity = Array(N[0], N[1], N[2]);
	DampedEnergy = Array(N[0],N[1],N[2]);
	//}
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

Array& SolutionData::GetResistivity()
{
	return Resistivity;
}

Array SolutionData::GetDampedEnergy()
{
	Region region;

	for (int d=0; d<3; ++d)
	{
		region.lower[d] = StartIndex[d];
		region.upper[d] = EndIndex[d];
	}
	return DampedEnergy[region];
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

void SolutionData::ComputeEnergy()
{
	TotalEnergy = 0;
	MagneticEnergy = 0;
	ElectricEnergy = 0;
	TotalOhmHeat = 0;
	for (int i = StartIndex[0]; i < EndIndex[0]; ++i)
	{
		for (int j = StartIndex[1]; j < EndIndex[1]; ++j)
		{
			for (int k = StartIndex[2]; k < EndIndex[2]; ++k)
			{
				for (int d=0; d<3; ++d) 
				{
					MagneticEnergy  += PField(i,j,k,d)*PField(i,j,k,d)*dx[0]*dx[1]*dx[2];
					ElectricEnergy  += PField(i,j,k,d+3)*PField(i,j,k,d+3)*dx[0]*dx[1]*dx[2];
				}
				if (SpongeLayer > 0 ){
					TotalOhmHeat += DampedEnergy(i,j,k);
				}

			} // end loop over k
		}// end loop over j
	}// end loop over i
	MagneticEnergy = MagneticEnergy/2.; 
	ElectricEnergy = ElectricEnergy/2.;
	TotalEnergy = MagneticEnergy + ElectricEnergy;
}

void SolutionData::ComputeMaxEdotB()
{
	MaxEdotB = 0;
	for (int i = StartIndex[0]; i < EndIndex[0]; ++i)
	{
		for (int j = StartIndex[1]; j < EndIndex[1]; ++j)
		{
			for (int k = StartIndex[2]; k < EndIndex[2]; ++k)
			{
				double EdotB = 0;
				for (int d=0; d<3; ++d) 
				{
					EdotB  += PField(i,j,k,d)*PField(i,j,k,d+3)*dx[0]*dx[1]*dx[2];
				}

				MaxEdotB = fmax(fabs(EdotB), MaxEdotB);
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
					divB += ((PField(i-2,j,k,0)-PField(i+2,j,k,0))/12. - (PField(i-1,j,k,0)-PField(i+1,j,k,0))*2./3.);
				}	
				if ( N[1] > 1 )
				{
					divB += ((PField(i,j-2,k,1)-PField(i,j+2,k,1))/12. - (PField(i,j-1,k,1)-PField(i,j+1,k,1))*2./3.);
				}
				if ( N[2] > 1 )
				{
					divB += ((PField(i,j,k-2,2)-PField(i,j,k+2,2))/12. - (PField(i,j,k-1,2)-PField(i,j,k+1,2))*2./3.);
				}	

				MaxdivB = fmax(fabs(divB), MaxdivB);
			} // end loop over k
		}// end loop over j
	}// end loop over i
}

void SolutionData::InitialData()
{
	// Initial Condition for test cases
	//TestCases();

	// 2D dipole field in x-y plane
	//Dipole2D();

	// Colliding Alfven packates
	//Equilibrium();
	AlfvenPacket3D();

	// Combine E and B field to Primitive field P
	CombineBEtoP();
	if ( SpongeLayer > 0){
		InitializeResistivity();
	}
}


void SolutionData::TestCases()
{
	double pi = 3.1415927;
	double kvec[3] = {0, 2*pi, 0};


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

				//Stationary Alfven Wave
				/*BField(i,j,k,0) = 1.;
				BField(i,j,k,1) = 1.;

				if ((x(i,j,k,0)-0.5)<=-0.1) BField(i,j,k,2) = 1;
				if ((x(i,j,k,0)-0.5)>-0.1 && (x(i,j,k,0)-0.5)<0.1) BField(i,j,k,2) = 1+0.15*(1+sin(5*3.1415927*(x(i,j,k,0)-0.5)));
				if ((x(i,j,k,0)-0.5)>0.1) BField(i,j,k,2) = 1.3;

				EField(i,j,k,2) = 1;
				EField(i,j,k,0) = -BField(i,j,k,2);*/

				// Fast Wave
				/*BField(i,j,k,0) = 0.1*cos(kx);
				BField(i,j,k,1) = 1;
				EField(i,j,k,2) = -BField(i,j,k,0);*/

				// Alfven Wave
				/*BField(i,j,k,0) = 0.1*cos(kx);
				BField(i,j,k,2) = 1;
				EField(i,j,k,1) = -BField(i,j,k,0);*/

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

				// Tearing Instability
				double phi = pi*(1+tanh((x(i,j,k,2)-0.5)/0.1));
				BField(i,j,k,0) = sin(phi);
				BField(i,j,k,1) = -cos(phi);

				BField(i,j,k,2) += 1e-1*sin(4*3.14159265*(x(i,j,k,0)-0.5));



			} // end loop over k
		}// end loop over j
	}// end loop over i
}

void SolutionData::AlfvenPacket1D()
{
	double pi = 3.1415927;
	//double kvec[3] = {2*pi, 0, 2*pi};
	double amp = AmplitudeOfAlfvenPacket;
	double spd = 0.01;
	double c1[3] = {0.5, 0.5, 0.25};
	double c2[3] = {0.5, 0.5, 0.75};
	double theta = pi/20.*ntheta;

	for (int i = 0; i < Nint[0]; ++i)
	{
		for (int j = 0; j < Nint[1]; ++j)
		{
			for (int k = 0; k < Nint[2]; ++k)
			{	
				double E1 = exp(-(x(i,j,k,2) - c1[2]) * (x(i,j,k,2) - c1[2])/spd);
				double E2 = exp(-(x(i,j,k,2) - c2[2]) * (x(i,j,k,2) - c2[2])/spd);

				E1 *= amp;
				E2 *= amp;
				// Create A Gaussian packate for two colliding Alfven waves

				// 3D case two cylinder
				
				
				EField(i,j,k,0) = E1;
                EField(i,j,k,0) += E2*cos(theta);
                EField(i,j,k,1) += E2*sin(theta);

				BField(i,j,k,1) = E1;
                BField(i,j,k,1) += -E2*cos(theta);
                BField(i,j,k,0) += E2*sin(theta);

				BField(i,j,k,2) = 1;

			} // end loop over k
		}// end loop over j
	}// end loop over i
}

void SolutionData::Random2D(){
	// generate radom field in 2D x-z plane
	//double pi = 3.1415927;
	//double kvec[3] = {2*pi, 0, 2*pi};
	int num = 100;
	double amp[num];
	double spd[num];
	double rx[num];
	double rz[num];
	int sign[num];

	for (int i=0; i<num; ++i)
	{
		amp[i] = (rand()/((double)RAND_MAX)-0.5)*0.8;
		rx[i] = rand()/((double)RAND_MAX);
		rz[i] = rand()/((double)RAND_MAX);
		spd[i] = 0.01+rand()/((double)RAND_MAX)*0.04;
		double p = rand()/((double)RAND_MAX);
		if (p>= 0.5)
			{ 
				sign[i] = 1;
			}else{
				sign[i] = -1;
			} 
	}


	for (int i = 0; i < Nint[0]; ++i)
	{
		for (int j = 0; j < Nint[1]; ++j)
		{
			for (int k = 0; k < Nint[2]; ++k)
			{	
				double E = 0;
				double B = 0;

				for (int p = 0; p<num; ++p){
					double xmin = fmin(fabs(x(i,j,k,0) - rx[p]),fmin(fabs(x(i,j,k,0) - rx[p] + 1),fabs(x(i,j,k,0) - rx[p] - 1)));
					double zmin = fmin(fabs(x(i,j,k,2) - rz[p]),fmin(fabs(x(i,j,k,2) - rz[p] + 1),fabs(x(i,j,k,2) - rz[p] - 1)));

					E += amp[p] * exp( - xmin * xmin/spd[p] )
								* exp( - zmin * zmin/spd[p] );
					B += sign[p] * amp[p] * exp( - xmin * xmin/spd[p] )
										  * exp( - zmin * zmin/spd[p] );
				}
				//printf("%f %f \n", E,B);

				EField(i,j,k,0) = E;
                BField(i,j,k,1) += B;
				BField(i,j,k,2) = 1;

			} // end loop over k
		}// end loop over j
	}// end loop over i
}

void SolutionData::AlfvenPacket2D()
{
	//double pi = 3.1415927;
	//double kvec[3] = {2*pi, 0, 2*pi};
	double amp = AmplitudeOfAlfvenPacket;
	double spd = 0.01;
	double c1[3] = {0.5, 0.5, 0.25};
	double c2[3] = {0.5, 0.5, 0.75};


	for (int i = 0; i < Nint[0]; ++i)
	{
		for (int j = 0; j < Nint[1]; ++j)
		{
			for (int k = 0; k < Nint[2]; ++k)
			{	
				double E1 = exp(-(x(i,j,k,2) - c1[2]) * (x(i,j,k,2) - c1[2])/spd)+(rand()/((double)RAND_MAX)-0.5)*0.1;
							//20*(x(i,j,k,0) - c1[0])
							//*exp(-(x(i,j,k,0) - c1[0]) * (x(i,j,k,0) - c1[0])/spd)
							//*exp(-(x(i,j,k,2) - c1[2]) * (x(i,j,k,2) - c1[2])/spd);
				double E2 = exp(-(x(i,j,k,2) - c2[2]) * (x(i,j,k,2) - c2[2])/spd)+(rand()/((double)RAND_MAX)-0.5)*0.1;
							//20*(x(i,j,k,0) - c2[0])
							//*exp(-(x(i,j,k,0) - c2[0]) * (x(i,j,k,0) - c2[0])/spd)
							//*exp(-(x(i,j,k,2) - c2[2]) * (x(i,j,k,2) - c2[2])/spd);

				E1 *= amp;
				E2 *= amp;
				// Create A Gaussian packate for two colliding Alfven waves

				// 3D case two cylinder
				
				
				EField(i,j,k,0) = -E1;
                EField(i,j,k,0) +=-E2;
				BField(i,j,k,1) = E2;
                BField(i,j,k,1) += -E1;
				BField(i,j,k,2) = 1 ;

			} // end loop over k
		}// end loop over j
	}// end loop over i
}

void SolutionData::AlfvenPacket3D()
{
	//double pi = 3.1415927;
	//double kvec[3] = {2*pi, 0, 2*pi};
	double amp = AmplitudeOfAlfvenPacket;
	double width = 0.2;
	double spd = width * width;
	double c1[3] = {0.5, 0.5, 0.25};
	double c2[3] = {0.5, 0.5, 0.75};

	for (int i = 0; i < Nint[0]; ++i)
	{
		for (int j = 0; j < Nint[1]; ++j)
		{
			for (int k = 0; k < Nint[2]; ++k)
			{	
				double r1 = 0;
				double r2 = 0;

				for (int d=0; d<3; ++d)
				{
					r1 += (x(i,j,k,d) - c1[d]) * (x(i,j,k,d) - c1[d]);
					r2 += (x(i,j,k,d) - c2[d]) * (x(i,j,k,d) - c2[d]);
				}

				// Create A Gaussian packate for two colliding Alfven waves

				// 3D case two cylinder
				double g1 = amp*exp(-r1/spd)/width;
				double g2 = amp*exp(-r2/spd)/width;
				
				BField(i,j,k,0) = 2*(x(i,j,k,1)-c1[1])*g1 + 2*(x(i,j,k,1)-c2[1])*g2;
                BField(i,j,k,1) = -2*(x(i,j,k,0)-c1[0])*g1 - 2*(x(i,j,k,0)-c2[0])*g2;
				BField(i,j,k,2) = 1;

				EField(i,j,k,0) = -2*(x(i,j,k,0)-c1[0])*g1 + 2*(x(i,j,k,0)-c2[0])*g2;
                EField(i,j,k,1) = -2*(x(i,j,k,1)-c1[1])*g1 + 2*(x(i,j,k,1)-c2[1])*g2;
			} // end loop over k
		}// end loop over j
	}// end loop over i
}

void SolutionData::Equilibrium()
{
	double pi = 3.1415926535897932384626;
	double kvec = 2*pi;

	for (int i = 0; i < Nint[0]; ++i)
	{
		for (int j = 0; j < Nint[1]; ++j)
		{
			for (int k = 0; k < Nint[2]; ++k)
			{					
				BField(i,j,k,2) = sqrt(2.) * sin(kvec*x(i,j,k,2)) * cos(kvec*x(i,j,k,0));
				BField(i,j,k,0) = -sqrt(2.) * cos(kvec*x(i,j,k,2)) * sin(kvec*x(i,j,k,0));
				BField(i,j,k,1) = -2. * sin(kvec*x(i,j,k,2)) * sin(kvec*x(i,j,k,0)) + sin(5*kvec*x(i,j,k,2))*0.5;

				//BField(i,j,k,1) +=  exp(-( x(i,j,k,2) - 0.25 ) * ( x(i,j,k,2) - 0.25 ) / 0.01)
				//				  * exp(-( x(i,j,k,0) - 0.5 ) * ( x(i,j,k,0) - 0.5 ) / 0.05);
				//BField(i,j,k,1) +=  exp(-( x(i,j,k,2) - 0.75 ) * ( x(i,j,k,2) - 0.75 ) / 0.01)
				//				  * exp(-( x(i,j,k,0) - 0.5 ) * ( x(i,j,k,0) - 0.5 ) / 0.05);

				//EField(i,j,k,0) +=  0.1*exp(-( x(i,j,k,2) - 0.25 ) * ( x(i,j,k,2) - 0.25 ) / 0.01)
				//				  * exp(-( x(i,j,k,0) - 0.5 ) * ( x(i,j,k,0) - 0.5 ) / 0.05);
				//EField(i,j,k,0) +=  exp(-( x(i,j,k,2) - 0.75 ) * ( x(i,j,k,2) - 0.75 ) / 0.01)
				//				  * exp(-( x(i,j,k,0) - 0.5 ) * ( x(i,j,k,0) - 0.5 ) / 0.05);


			} // end loop over k
		}// end loop over j
	}// end loop over i
}


void SolutionData::InitializeResistivity()
{
	for (int i = StartIndex[0]; i < EndIndex[0]; ++i)
	{
		for (int j = StartIndex[1]; j < EndIndex[1]; ++j)
		{
			for (int k = StartIndex[2]; k < EndIndex[2]; ++k)
			{	
				// 2D case: absorbing in X direction
				double xsig = (fabs(x(i-StartIndex[0],j-StartIndex[1],k-StartIndex[2],0)-0.5)-0.3)/0.2;
				if ( GridDim == 3){
					// 3D case: absorbing in X-Y direction
					xsig = pow(x(i-StartIndex[0],j-StartIndex[1],k-StartIndex[2],0)-0.5,2)
				            	+ pow(x(i-StartIndex[0],j-StartIndex[1],k-StartIndex[2],1)-0.5,2);
					xsig = (xsig-0.3)/0.2;
					
				} 
				if ( xsig >= 0 ){
					Resistivity(i,j,k) = 1. - exp(-8. * pow( xsig, 4));
				}

			} // end loop over k
		}// end loop over j
	}// end loop over i
	// Set resivitity, currently used only for sponge layer
				
}

void SolutionData::ComputeDampedEnergy(double dt)
{
	for (int i = StartIndex[0]; i < EndIndex[0]; ++i)
	{
		for (int j = StartIndex[1]; j < EndIndex[1]; ++j)
		{
			for (int k = StartIndex[2]; k < EndIndex[2]; ++k)
			{
				for (int d=0; d<3; ++d) 
				{
					DampedEnergy(i,j,k)  += Resistivity(i,j,k)*PField(i,j,k,d+3)*PField(i,j,k,d+3)*dx[0]*dx[1]*dx[2]*dt;
				}

			} // end loop over k
		}// end loop over j
	}// end loop over i

}
