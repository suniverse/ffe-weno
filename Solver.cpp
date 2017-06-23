#include "Grid.hpp"
#include "Roe.hpp"
#include "Reconstruction.hpp"
#include <cmath>
#include <array>

using namespace Cow;

Solver::Solver(UserParameters userParameters)
{
	int NumGhost = userParameters.NumberOfGhostZones;
	double NumCFL = userParameters.CFLSafetyNumber;

	ReconstructionMethod = userParameters.ReconstructionMethod;
	SetReconstructionMethod();

	// loop over three dimensions
	for (int d=0; d<3; ++d)
	{
		// set number of grids in that dimension
		N[d] = userParameters.GridLength[d];

		if (N[d]==1) continue;

		// set start and end index in that dimension (exclude ghost zones)
		StartIndex[d] = NumGhost;
		EndIndex[d] = N[d] - NumGhost;

		// set dx in that direction
		dx[d] = 1./(N[d]-2*NumGhost);
		dt = fmin(dx[d]*NumCFL, dt);

	}
	/* Fluxes: they are indexed on the left face of a cell
	   FluxL(i,j,k,0) lives at the left edge of the cell (i,j,k) */
	// left flux
	//FluxL  = Array(N[0], N[1], N[2],6);
	// right flux
	//FluxR  = Array(N[0], N[1], N[2],6);
	// Godunov flux
	FluxG     = Array(N[0], N[1], N[2],6);
	// change of the primitive field
	DeltaP    = Array(N[0], N[1], N[2],6);
	// auxillary field to store primitive varible (B,E)
	P0        = Array(N[0], N[1], N[2],6);
	// auxillary field to clean divB
	Psi    	  = Array(N[0], N[1], N[2]);
	Psi0      = Array(N[0], N[1], N[2]);
	DeltaPsi  = Array(N[0], N[1], N[2]);

	//reconstruct = Reconstruction();// implied!
	//reconstruct = Reconstruction(arg1, arg2);
}

void Solver::ComputeDeltaP(const Array& P)
{	
	AddSourceTerm(P);
	DednerDamp(P);

	if ( N[0] > 1)
	{
		SweepX(P);
	}
	if ( N[1] > 1)
	{	
		SweepY(P);
	}	
	if ( N[2] > 1)
	{
		SweepZ(P);
	}
}

void Solver::AddSourceTerm(const Array& P)
{	

	double rho, B2;

	// update j_pert in the source
	for (int i = StartIndex[0]; i < EndIndex[0]; ++i)
	{
		for (int j = StartIndex[1]; j < EndIndex[1]; ++j)
		{
			for (int k = StartIndex[2]; k < EndIndex[2]; ++k)
			{
				// Calculate rho
				rho = 0;
				if ( N[0] > 1 )
				{
					rho += ((P(i-2,j,k,3)-P(i+2,j,k,3))/12. - (P(i-1,j,k,3)-P(i+1,j,k,3))*2./3.)/dx[0];
					//rho +=  - (P(i-1,j,k,3)-P(i+1,j,k,3))/2./dx[0];
				}	
				if ( N[1] > 1 )
				{
					rho += ((P(i,j-2,k,4)-P(i,j+2,k,4))/12. - (P(i,j-1,k,4)-P(i,j+1,k,4))*2./3.)/dx[1];
				}
				if ( N[2] > 1 )
				{
					rho += ((P(i,j,k-2,5)-P(i,j,k+2,5))/12. - (P(i,j,k-1,5)-P(i,j,k+1,5))*2./3.)/dx[2];
				}	

				//printf("rho %d %f\n", i,rho);
				// calculabe B*B
				B2 = 0;
				for (int d=0; d<3; d++)
				{
					B2 += P(i,j,k,d)*P(i,j,k,d);
				}
				//printf("rho B2 =%f %f\n", rho, B2);
				// Calculate E X B
				// put everything together
				DeltaP(i,j,k,0) = 0.;
				DeltaP(i,j,k,1) = 0.;
				DeltaP(i,j,k,2) = 0.;
				DeltaP(i,j,k,3) = -rho / B2 * ( P(i,j,k,4)*P(i,j,k,2) - P(i,j,k,5)*P(i,j,k,1));
				DeltaP(i,j,k,4) = -rho / B2 * ( P(i,j,k,5)*P(i,j,k,0) - P(i,j,k,3)*P(i,j,k,2));
				DeltaP(i,j,k,5) = -rho / B2 * ( P(i,j,k,3)*P(i,j,k,1) - P(i,j,k,4)*P(i,j,k,0));
			}
		}
	}
}

void Solver::DednerDamp(const Array& P)
{	
	double ch2 = 1.0;
	double tau = 20*dt;

	double divB;

	// compute divB
	for (int i = StartIndex[0]; i < EndIndex[0]; ++i)
	{
		for (int j = StartIndex[1]; j < EndIndex[1]; ++j)
		{
			for (int k = StartIndex[2]; k < EndIndex[2]; ++k)
			{
				// Calculate divB
				divB = 0;
				if ( N[0] > 1 )
				{
					divB += ((P(i-2,j,k,0)-P(i+2,j,k,0))/12. - (P(i-1,j,k,0)-P(i+1,j,k,0))*2./3.)/dx[0];
					//rho +=  - (P(i-1,j,k,3)-P(i+1,j,k,3))/2./dx[0];
				}	
				if ( N[1] > 1 )
				{
					divB += ((P(i,j-2,k,1)-P(i,j+2,k,1))/12. - (P(i,j-1,k,1)-P(i,j+1,k,1))*2./3.)/dx[1];
				}
				if ( N[2] > 1 )
				{
					divB += ((P(i,j,k-2,2)-P(i,j,k+2,2))/12. - (P(i,j,k-1,2)-P(i,j,k+1,2))*2./3.)/dx[2];
				}	

				DeltaPsi(i,j,k) = -ch2*divB - Psi(i,j,k)/tau;

				if ( N[0] > 1 )
				{
					DeltaP(i,j,k,0) = -((Psi(i-2,j,k)-Psi(i+2,j,k))/12. - (Psi(i-1,j,k)-Psi(i+1,j,k))*2./3.)/dx[0];
					//rho +=  - (P(i-1,j,k,3)-P(i+1,j,k,3))/2./dx[0];
				}	
				if ( N[1] > 1 )
				{
					DeltaP(i,j,k,1) = -((Psi(i,j-2,k)-Psi(i,j+2,k))/12. - (Psi(i,j-1,k)-Psi(i,j+1,k))*2./3.)/dx[1];
				}
				if ( N[2] > 1 )
				{
					DeltaP(i,j,k,2) = -((Psi(i,j,k-2)-Psi(i,j,k+2))/12. - (Psi(i,j,k-1)-Psi(i,j,k+1))*2./3.)/dx[2];
				}	

			}
		}
	}
}

void Solver::SweepX(const Array& P)
{
	// loop over each cell to get the Godunov flux on the x-axis

	std::array<double,6> leftstate;
	std::array<double,6> rightstate;
	std::array<double,6> fluxGod;

	for (int i = StartIndex[0]; i <= EndIndex[0]; ++i) // Fluxes have one more cell than the Primitive Vars
	{
		for (int j = StartIndex[1]; j < EndIndex[1]; ++j)
		{
			for (int k = StartIndex[2]; k < EndIndex[2]; ++k)
			{
				// convert primitive vars to std::array to pass into Roe solver
				for (int d=0; d<6; ++d) 
				{
					// prepare stencil
					double arr[6];
					for (int n=-3; n<3; ++n)
					{
						arr[n+3] = P(i+n,j,k,d);
					}
					// reconstruct left and right state
					leftstate[d]  = reconstruct.reconstruct(&arr[2], reconstructModeC2R);
					rightstate[d] = reconstruct.reconstruct(&arr[3], reconstructModeC2L);

				}

				fluxGod = xRoe(leftstate, rightstate);

				// convert godnov flux back to Cow::Array
				for (int d=0; d<6; ++d) 
				{
					FluxG(i,j,k,d)  = fluxGod[d];
				}
			} // end loop over k
		}// end loop over j
	}// end loop over i

	for (int i = StartIndex[0]; i < EndIndex[0]; ++i)
	{
		for (int j = StartIndex[1]; j < EndIndex[1]; ++j)
		{
			for (int k = StartIndex[2]; k < EndIndex[2]; ++k)
			{
				// convert primitive vars to std::array to pass into Roe solver
				for (int d=0; d<6; ++d) 
				{
					DeltaP(i,j,k,d) += (FluxG(i,j,k,d)-FluxG(i+1,j,k,d))/dx[0];
				}

			} // end loop over k
		}// end loop over j
	}// end loop over i
}

void Solver::SweepY(const Array& P)
{
	// loop over each cell to get the Godunov flux on the y-axis

	std::array<double,6> leftstate;
	std::array<double,6> rightstate;
	std::array<double,6> fluxGod;

	for (int i = StartIndex[0]; i < EndIndex[0]; ++i) // Fluxes have one more cell than the Primitive Vars
	{
		for (int j = StartIndex[1]; j <= EndIndex[1]; ++j)
		{
			for (int k = StartIndex[2]; k < EndIndex[2]; ++k)
			{
				// convert primitive vars to std::array to pass into Roe solver
				for (int d=0; d<6; ++d) 
				{
					// prepare stencil
					double arr[6];
					for (int n=-3; n<3; ++n)
					{
						arr[n+3] = P(i,j+n,k,d);
					}

					// reconstruct left and right state
					leftstate[d]  = reconstruct.reconstruct(&arr[2], reconstructModeC2R);
					rightstate[d] = reconstruct.reconstruct(&arr[3], reconstructModeC2L);
				}

				fluxGod = yRoe(leftstate, rightstate);

				// convert godnov flux back to Cow::Array
				for (int d=0; d<6; ++d) 
				{
					FluxG(i,j,k,d)  = fluxGod[d];
				}
			} // end loop over k
		}// end loop over j
	}// end loop over i

	for (int i = StartIndex[0]; i < EndIndex[0]; ++i)
	{
		for (int j = StartIndex[1]; j < EndIndex[1]; ++j)
		{
			for (int k = StartIndex[2]; k < EndIndex[2]; ++k)
			{
				// convert primitive vars to std::array to pass into Roe solver
				for (int d=0; d<6; ++d) 
				{
					DeltaP(i,j,k,d)  += (FluxG(i,j,k,d)-FluxG(i,j+1,k,d))/dx[1];
				}

			} // end loop over k
		}// end loop over j
	}// end loop over i
}

void Solver::SweepZ(const Array& P)
{
	// loop over each cell to get the Godunov flux on the z-axis

	std::array<double,6> leftstate;
	std::array<double,6> rightstate;
	std::array<double,6> fluxGod;

	for (int i = StartIndex[0]; i < EndIndex[0]; ++i) // Fluxes have one more cell than the Primitive Vars
	{
		for (int j = StartIndex[1]; j < EndIndex[1]; ++j)
		{
			for (int k = StartIndex[2]; k <= EndIndex[2]; ++k)
			{
				// convert primitive vars to std::array to pass into Roe solver
				for (int d=0; d<6; ++d) 
				{
					// prepare stencil
					double arr[6];
					for (int n=-3; n<3; ++n)
					{
						arr[n+3] = P(i,j,k+n,d);
					}

					// reconstruct left and right state
					leftstate[d]  = reconstruct.reconstruct(&arr[2], reconstructModeC2R);
					rightstate[d] = reconstruct.reconstruct(&arr[3], reconstructModeC2L);
				}

				fluxGod = zRoe(leftstate, rightstate);

				// convert godnov flux back to Cow::Array
				for (int d=0; d<6; ++d) 
				{
					FluxG(i,j,k,d)  = fluxGod[d];
				}
			} // end loop over k
		}// end loop over j
	}// end loop over i

	for (int i = StartIndex[0]; i < EndIndex[0]; ++i)
	{
		for (int j = StartIndex[1]; j < EndIndex[1]; ++j)
		{
			for (int k = StartIndex[2]; k < EndIndex[2]; ++k)
			{
				// convert primitive vars to std::array to pass into Roe solver
				for (int d=0; d<6; ++d) 
				{
					DeltaP(i,j,k,d)  += (FluxG(i,j,k,d)-FluxG(i,j,k+1,d))/dx[2];
				}

			} // end loop over k
		}// end loop over j
	}// end loop over i
}

void Solver::CleanEdotB( Array& P )
{
	double EdotB, B2;
	// update j_pert in the source
	for (int i = StartIndex[0]; i < EndIndex[0]; ++i)
	{
		for (int j = StartIndex[1]; j < EndIndex[1]; ++j)
		{
			for (int k = StartIndex[2]; k < EndIndex[2]; ++k)
			{
				// Calculate EdotB
				EdotB = 0;
				// Calculate B^2
				B2 = 0;

				for (int d=0; d<3; ++d)
				{
					EdotB += P(i,j,k,d)*P(i,j,k,d+3);
					B2 += P(i,j,k,d)*P(i,j,k,d);				
				}

				for (int d=0; d<3; ++d)
				{
				// update parallel E
					P(i,j,k,d+3) -= EdotB/B2*P(i,j,k,d);

				}
			}
		}
	}
}

void Solver::ShrinkE( Array& P)
{
	for (int i = StartIndex[0]; i < EndIndex[0]; ++i)
	{
		for (int j = StartIndex[1]; j < EndIndex[1]; ++j)
		{
			for (int k = StartIndex[2]; k < EndIndex[2]; ++k)
			{
				// Calculate B^2
				double B2 = 0;
				// Calculate E^2
				double E2 = 0;

				for (int d=0; d<3; ++d)
				{
					B2 += P(i,j,k,d)*P(i,j,k,d);
					E2 += P(i,j,k,d+3)*P(i,j,k,d+3);				
				}

				if (E2 > B2)
				{
					for (int d=0; d<3; ++d)
					{
						P(i,j,k,d+3) = sqrt(B2/E2)*P(i,j,k,d+3);				
					}
				}
			}
		}
	}
}

void Solver::Advance(Array& P)
{ 
	//TVD-rk3 evolution
	Cache(P, Psi);


	// 1st step
	ComputeDeltaP(P);

	TVDstep(P, P0, DeltaP, 0.0, 1.0);
	TVDstep(Psi, Psi0, DeltaPsi, 0.0, 1.0);

	// add J_parallel by removing EdotB
	CleanEdotB(P);


	SetBoundaryValue(P);

	// 2nd step
	ComputeDeltaP(P);

	TVDstep(P, P0, DeltaP, 3./4., 1./4.);
	TVDstep(Psi, Psi0, DeltaPsi, 3./4., 1./4.);

	// add J_parallel by removing EdotB
	CleanEdotB(P);


	//TVDstep(P, 1./2., 1./2.);


	SetBoundaryValue(P);

	// 3rd step
	ComputeDeltaP(P);

	TVDstep(P, P0, DeltaP, 1./3., 2./3.);
	TVDstep(Psi, Psi0, DeltaPsi, 1./3., 2./3.);

	// add J_parallel by removing EdotB
	CleanEdotB(P);


	ShrinkE(P);

	SetBoundaryValue(P);

}

void Solver::Cache(const Array& P, const Array& Psi)
{
	// copy B and E to auxillary variables.
	P0 = P;
	Psi0 = Psi;
}

void Solver::TVDstep(Array& P, const Array& P0, const Array& DeltaP, double alpha, double beta)
{
	for (int n = 0; n < P.size(); ++n)
	{
		P[n] = alpha * P0[n] + beta * (P[n] + dt * DeltaP[n]);
	}

}

void Solver::SetBoundaryValue(Array& P)
{
	for (int d = 0; d < 3; ++d)
	{
		if (N[d] == 1) continue;

		Region guardL, guardR, interiorL, interiorR;

		guardL.lower[d] = 0;
		guardL.upper[d] = StartIndex[d];

		interiorL.lower[d] = StartIndex[d];
		interiorL.upper[d] = StartIndex[d] * 2;

		guardR.lower[d] = EndIndex[d];
		guardR.upper[d] = N[d];

		interiorR.lower[d] = EndIndex[d] - StartIndex[d];
		interiorR.upper[d] = EndIndex[d];

		P[guardL] = P[interiorR];
		P[guardR] = P[interiorL];

	}
}

void Solver::Setdt(double dt_input)
{
	dt = dt_input;
}

void Solver::SetReconstructionMethod()
{	
	switch (ReconstructionMethod)
	{
		case 1: // PLM
			reconstructModeC2R = Reconstruction::PLM_C2R;
			reconstructModeC2L = Reconstruction::PLM_C2L;
			break;
		case 2: // WENO5
			reconstructModeC2R = Reconstruction::WENO5_FD_C2R;
			reconstructModeC2L = Reconstruction::WENO5_FD_C2L;
			break;
	}
					
}
