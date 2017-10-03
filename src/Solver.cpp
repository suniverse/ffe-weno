#include "Grid.hpp"
#include "Roe.hpp"
#include "Reconstruction.hpp"
#include <cmath>
#include <array>

using namespace Cow;

Solver::Solver(UserParameters userParameters)
{
	NumGhost = userParameters.NumberOfGhostZones;
	double NumCFL = userParameters.CFLSafetyNumber;

	ReconstructionMethod = userParameters.ReconstructionMethod;
	SetReconstructionMethod();

	// loop over three dimensions
	for (int d=0; d<3; ++d)
	{
		// set number of grids in that dimension
		GridLength[d] = userParameters.GridLength[d];

		if (GridLength[d]==1) continue;

		// set start and end index in that dimension (exclude ghost zones)
		StartIndex[d] = NumGhost;
		EndIndex[d] = GridLength[d] - NumGhost;

		// set dx in that direction
		dx[d] = userParameters.dx[d];
		dt = fmin(dx[d]*NumCFL, dt);
	}
	
	// On/off sponge layer
	SpongeLayer = userParameters.SpongeLayer;

	/* Fluxes: they are indexed on the left face of a cell
	   FluxL(i,j,k,0) lives at the left edge of the cell (i,j,k) */

	// Godunov flux
	FluxG     = Array(GridLength[0], GridLength[1], GridLength[2],6);
	// change of the primitive field
	DeltaP    = Array(GridLength[0], GridLength[1], GridLength[2],6);
	// auxillary field to store primitive varible (B,E)
	P0        = Array(GridLength[0], GridLength[1], GridLength[2],6);
	// auxillary field to clean divB
	Psi    	  = Array(GridLength[0], GridLength[1], GridLength[2]);
	Psi0      = Array(GridLength[0], GridLength[1], GridLength[2]);
	DeltaPsi  = Array(GridLength[0], GridLength[1], GridLength[2]);

}

void Solver::ComputeDeltaP(const Array& P, const Array& Resistivity)
{	/* compute detla P:
	1. Initialize to set deltaP 0
	2. Compute curl E and curl B, the linear operator in the maxwell equation.
	   The Gudnov flux is calculated using various interpolation method and Roe solver
	3. Compute J_parallel, using curl E and curl B computed above
	4. Compute J_pert, finite differencing to get div E
	5. Add Dedner Damping to do hyperbolic cleaning*/
	InitializeDeltaP(P);

	if ( GridLength[0] > 1)
	{
		SweepX(P);
	}
	if ( GridLength[1] > 1)
	{	
		SweepY(P);
	}	
	if ( GridLength[2] > 1)
	{
		SweepZ(P);
	}
	AddJparallel(P);
	AddJpert(P);
	DednerDamp(P);

	// Add resistive term for sponge layer
	if ( SpongeLayer > 0 ){
		for (int i = StartIndex[0]; i < EndIndex[0]; ++i)
		{
			for (int j = StartIndex[1]; j < EndIndex[1]; ++j)
			{
				for (int k = StartIndex[2]; k < EndIndex[2]; ++k)
				{
					DeltaP(i,j,k,3) += - Resistivity(i,j,k) * P(i,j,k,3);
					DeltaP(i,j,k,4) += - Resistivity(i,j,k) * P(i,j,k,4);
					DeltaP(i,j,k,5) += - Resistivity(i,j,k) * P(i,j,k,5);
				}
			}
		}
	}


}

void Solver::InitializeDeltaP(const Array& P)
{
	for (int n = 0; n < DeltaP.size(); ++n)
	{
		DeltaP[n] = 0;
	}
}

void Solver::AddJparallel(const Array& P)
{	
/* Add J parallel to the DeltaP, this must be done after sweep in all 3 direction, 
   since it uses curlB and curlE calculated in the sweep routine. 
   J parallel is modified fro mideal FFE to damp E dot B, tau is the damping time. */
	double b, e, B2, EB;
	double tau = 2*dt;

	// update j_pert in the source
	for (int i = StartIndex[0]; i < EndIndex[0]; ++i)
	{
		for (int j = StartIndex[1]; j < EndIndex[1]; ++j)
		{
			for (int k = StartIndex[2]; k < EndIndex[2]; ++k)
			{
				b = 0;
				e = 0;
				// E dot B
				EB = 0;
				//printf("rho %d %f\n", i,rho);
				B2 = 0;
				for (int d=0; d<3; ++d)
				{
					// B^2
					B2 += P(i,j,k,d)*P(i,j,k,d);
					// E dot B
					EB += P(i,j,k,d)*P(i,j,k,d+3);
					// B dot curl B
					b += P(i,j,k,d)*DeltaP(i,j,k,d+3);
					// -E dot curl E
					e += P(i,j,k,d+3)*DeltaP(i,j,k,d);
				}
				//printf("rho B2 =%f %f\n", rho, B2);
				// Calculate E X B
				// put everything together
				//DeltaP(i,j,k,0) += 0.;
				//DeltaP(i,j,k,1) += 0.;
				//DeltaP(i,j,k,2) += 0.;
				DeltaP(i,j,k,3) += - ((b+e) + EB / tau ) / B2 * P(i,j,k,0);
				DeltaP(i,j,k,4) += - ((b+e) + EB / tau ) / B2 * P(i,j,k,1);
				DeltaP(i,j,k,5) += - ((b+e) + EB / tau ) / B2 * P(i,j,k,2);
			}
		}
	}
}

void Solver::AddJpert(const Array& P)
{	// compute J_pert, modified for E dot B non-zero, see Parfrey et al. 2017

	double rho, E2, B2, EB, den;

	// update j_pert in the source
	for (int i = StartIndex[0]; i < EndIndex[0]; ++i)
	{
		for (int j = StartIndex[1]; j < EndIndex[1]; ++j)
		{
			for (int k = StartIndex[2]; k < EndIndex[2]; ++k)
			{
				// Calculate rho
				rho = 0;
				if ( GridLength[0] > 1 )
				{
					rho += ((P(i-2,j,k,3)-P(i+2,j,k,3))/12. - (P(i-1,j,k,3)-P(i+1,j,k,3))*2./3.)/dx[0];
					//rho +=  - (P(i-1,j,k,3)-P(i+1,j,k,3))/2./dx[0];
				}	
				if ( GridLength[1] > 1 )
				{
					rho += ((P(i,j-2,k,4)-P(i,j+2,k,4))/12. - (P(i,j-1,k,4)-P(i,j+1,k,4))*2./3.)/dx[1];
				}
				if ( GridLength[2] > 1 )
				{
					rho += ((P(i,j,k-2,5)-P(i,j,k+2,5))/12. - (P(i,j,k-1,5)-P(i,j,k+1,5))*2./3.)/dx[2];
				}	

				//printf("rho %d %f\n", i,rho);
				// calculabe B*B
				B2 = 0;
				E2 = 0;
				EB = 0;
				for (int d=0; d<3; d++)
				{
					B2 += P(i,j,k,d)*P(i,j,k,d);
					E2 += P(i,j,k,d+3)*P(i,j,k,d+3);
					EB += P(i,j,k,d)*P(i,j,k,d+3);
				}
				den = E2 + 0.5 * ( B2-E2 + sqrt( (B2-E2)*(B2-E2) + 4*EB*EB ));
				//printf("rho B2 =%f %f\n", rho, B2);
				// Calculate E X B
				// put everything together
				DeltaP(i,j,k,0) += 0.;
				DeltaP(i,j,k,1) += 0.;
				DeltaP(i,j,k,2) += 0.;
				DeltaP(i,j,k,3) += -rho * ( P(i,j,k,4)*P(i,j,k,2) - P(i,j,k,5)*P(i,j,k,1)) / den;
				DeltaP(i,j,k,4) += -rho * ( P(i,j,k,5)*P(i,j,k,0) - P(i,j,k,3)*P(i,j,k,2)) / den;
				DeltaP(i,j,k,5) += -rho * ( P(i,j,k,3)*P(i,j,k,1) - P(i,j,k,4)*P(i,j,k,0)) / den;
			}
		}
	}
}

void Solver::DednerDamp(const Array& P)
{	
	// Hyperbolic clean of div B
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
				if ( GridLength[0] > 1 )
				{
					divB += ((P(i-2,j,k,0)-P(i+2,j,k,0))/12. - (P(i-1,j,k,0)-P(i+1,j,k,0))*2./3.)/dx[0];
					//rho +=  - (P(i-1,j,k,3)-P(i+1,j,k,3))/2./dx[0];
				}	
				if ( GridLength[1] > 1 )
				{
					divB += ((P(i,j-2,k,1)-P(i,j+2,k,1))/12. - (P(i,j-1,k,1)-P(i,j+1,k,1))*2./3.)/dx[1];
				}
				if ( GridLength[2] > 1 )
				{
					divB += ((P(i,j,k-2,2)-P(i,j,k+2,2))/12. - (P(i,j,k-1,2)-P(i,j,k+1,2))*2./3.)/dx[2];
				}	

				DeltaPsi(i,j,k) = -ch2*divB - Psi(i,j,k)/tau;

				if ( GridLength[0] > 1 )
				{
					DeltaP(i,j,k,0) += -((Psi(i-2,j,k)-Psi(i+2,j,k))/12. - (Psi(i-1,j,k)-Psi(i+1,j,k))*2./3.)/dx[0];
					//rho +=  - (P(i-1,j,k,3)-P(i+1,j,k,3))/2./dx[0];
				}	
				if ( GridLength[1] > 1 )
				{
					DeltaP(i,j,k,1) += -((Psi(i,j-2,k)-Psi(i,j+2,k))/12. - (Psi(i,j-1,k)-Psi(i,j+1,k))*2./3.)/dx[1];
				}
				if ( GridLength[2] > 1 )
				{
					DeltaP(i,j,k,2) += -((Psi(i,j,k-2)-Psi(i,j,k+2))/12. - (Psi(i,j,k-1)-Psi(i,j,k+1))*2./3.)/dx[2];
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
{ // clean E dot B, not used
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
{ // damp E where B^2 < E^2
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

void Solver::Advance(Cart cart, Array& P, const Array& Resistivity)
{ 
	//TVD-rk3 evolution
	Cache(P, Psi);


	// 1st step
	ComputeDeltaP(P, Resistivity);
	TVDstep(P, P0, DeltaP, 0.0, 1.0);
	TVDstep(Psi, Psi0, DeltaPsi, 0.0, 1.0);

	SetBoundaryValue(cart, P);
	SetBoundaryValue(cart, Psi);
	
	// 2nd step
	ComputeDeltaP(P, Resistivity);

	TVDstep(P, P0, DeltaP, 3./4., 1./4.);
	TVDstep(Psi, Psi0, DeltaPsi, 3./4., 1./4.);

	SetBoundaryValue(cart, P);
	SetBoundaryValue(cart, Psi);

	// 3rd step
	ComputeDeltaP(P, Resistivity);
	
	TVDstep(P, P0, DeltaP, 1./3., 2./3.);
	TVDstep(Psi, Psi0, DeltaPsi, 1./3., 2./3.);

	ShrinkE(P);

	SetBoundaryValue(cart, P);
	SetBoundaryValue(cart, Psi);

}

void Solver::Cache(const Array& P, const Array& Psi)
{
	// copy B and E to auxillary variables.
	P0 = P;
	Psi0 = Psi;
}

void Solver::TVDstep(Array& P, const Array& P0, const Array& DeltaP, double alpha, double beta)
{// TVD runge kutta step
	for (int n = 0; n < P.size(); ++n)
	{
		P[n] = alpha * P0[n] + beta * (P[n] + dt * DeltaP[n]);
	}

}

void Solver::SetBoundaryValue(Cart cart, Array& P)
{// Set values in the ghost zones for each patch
	cart.apply(P, NumGhost, cart.GlobalDim);

	/*for (int d = 0; d < 3; ++d)
	{
		if (GridLength[d] == 1) continue;

		Region guardL, guardR, interiorL, interiorR;

		guardL.lower[d] = 0;
		guardL.upper[d] = StartIndex[d];

		interiorL.lower[d] = StartIndex[d];
		interiorL.upper[d] = StartIndex[d] * 2;

		guardR.lower[d] = EndIndex[d];
		guardR.upper[d] = GridLength[d];

		interiorR.lower[d] = EndIndex[d] - StartIndex[d];
		interiorR.upper[d] = EndIndex[d];

		P[guardL] = P[interiorL];
		P[guardR] = P[interiorR];
	}*/
}

void Solver::Setdt(double dt_input)
{ // set dt for evolution
	dt = dt_input;
}

void Solver::SetReconstructionMethod()
{	// set the interpolation method
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

