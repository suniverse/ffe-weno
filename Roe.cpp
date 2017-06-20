#include <cmath>
#include <array>

std::array<double,6> xRoe(std::array<double,6> pl, std::array<double, 6> pr)
{
	// Convert prim to left and right flux
	std::array<double,6> fl =  {{0.,-pl[5],pl[4],0.,pl[2],-pl[1]}};
	std::array<double,6> fr =  {{0.,-pr[5],pr[4],0.,pr[2],-pr[1]}};

	// six eigenvalues
	std::array<double,6> lambda =  {{0.,1.,1.,0.,-1.,-1.}};
	std::array<double,6> alpha = {{0.,0.,0.,0.,0.,0.}};

	// six eigenvectors correspondingly
	std::array<std::array<double,6> ,6> eigvec = 
	{{
		{{ 1,           0,           0, 0,           0,          0 }},
		{{ 0, -1./sqrt(2),           0, 0,           0, 1./sqrt(2) }},
		{{ 0,           0,  1./sqrt(2), 0,  1./sqrt(2),          0 }},
		{{ 0,           0,           0, 1,           0,          0 }},
		{{ 0,           0, -1./sqrt(2), 0,  1./sqrt(2),          0 }},
		{{ 0,  1./sqrt(2),           0, 0,           0, 1./sqrt(2) }}
		}};
	
	// Gudonov Flux
	std::array<double,6> fg = {{0,0,0,0,0,0}};

	// project flux difference to the eigenvectors
	for (int i=0; i<6; ++i)
	{
		for (int j=0; j<6; ++j)
		{
			alpha[i] += eigvec[i][j] * (pr[j]-pl[j]);
		}
	}

 	// Take Roe average
	for (int i=0; i<6; ++i)
	{
		for (int j=0; j<6; ++j)
		{
			fg[i] -= fabs(lambda[j])*alpha[j]*eigvec[j][i];
		}
		fg[i] = (fg[i] + fl[i] + fr[i])/2.;
		//printf("left flux %f, right flux %f, godnov flux %f\n", fl[i], fr[i], fg[i]);

	}

	return fg;
}

std::array<double,6> yRoe(std::array<double,6> pl, std::array<double, 6> pr)
{
	// Convert prim to left and right flux
	std::array<double,6> fl =  {{pl[5], 0., -pl[3], -pl[2], 0., pl[0]}};
	std::array<double,6> fr =  {{pr[5], 0., -pr[3], -pr[2], 0., pr[0]}};

	// six eigenvalues
	std::array<double,6> lambda =  {{1.,0.,1.,-1.,0.,-1.}};
	std::array<double,6> alpha = {{0.,0.,0.,0.,0.,0.}};

	// six eigenvectors correspondingly
	std::array<std::array<double,6> ,6> eigvec = 
	{{
		{{           0, 0, -1./sqrt(2), 1./sqrt(2), 0,          0 }},
		{{           0, 1,           0,          0, 0,          0 }},
		{{  1./sqrt(2), 0,           0,          0, 0, 1./sqrt(2) }},
		{{           0, 0,  1./sqrt(2), 1./sqrt(2), 0,          0 }},
		{{           0, 0,           0,          0, 1,          0 }},
		{{ -1./sqrt(2), 0,           0,          0, 0, 1./sqrt(2) }}
		}};
	
	// Gudonov Flux
	std::array<double,6> fg = {{0,0,0,0,0,0}};

	// project flux difference to the eigenvectors
	for (int i=0; i<6; ++i)
	{
		for (int j=0; j<6; ++j)
		{
			alpha[i] += eigvec[i][j] * (pr[j]-pl[j]);
		}
	}

 	// Take Roe average
	for (int i=0; i<6; ++i)
	{
		for (int j=0; j<6; ++j)
		{
			fg[i] -= fabs(lambda[j])*alpha[j]*eigvec[j][i];
		}
		fg[i] = (fg[i] + fl[i] + fr[i])/2.;
	}

	return fg;
}

std::array<double,6> zRoe(std::array<double,6> pl, std::array<double, 6> pr)
{
	// Convert prim to left and right flux
	std::array<double,6> fl =  {{-pl[4], pl[3], 0, pl[1], -pl[0], 0}};
	std::array<double,6> fr =  {{-pr[4], pr[3], 0, pr[1], -pr[0], 0}};

	// six eigenvalues
	std::array<double,6> lambda =  {{-1.,-1.,1.,1.,0.,0.}};
	std::array<double,6> alpha = {{0.,0.,0.,0.,0.,0.}};

	// six eigenvectors correspondingly
	std::array<std::array<double,6> ,6> eigvec = 
	{{
		{{  1./sqrt(2),           0, 0,          0, 1./sqrt(2), 0 }},
		{{           0, -1./sqrt(2), 0, 1./sqrt(2),          0, 0 }},
		{{ -1./sqrt(2),           0, 0,          0, 1./sqrt(2), 0 }},
		{{           0,  1./sqrt(2), 0, 1./sqrt(2),          0, 0 }},
		{{           0,           0, 0,          0,          0, 1 }},
		{{           0,           0, 1,          0,          0, 0 }}
		}};
	
	// Gudonov Flux
	std::array<double,6> fg = {{0,0,0,0,0,0}};

	// project flux difference to the eigenvectors
	for (int i=0; i<6; ++i)
	{
		for (int j=0; j<6; ++j)
		{
			alpha[i] += eigvec[i][j] * (pr[j]-pl[j]);
		}
	}

 	// Take Roe average
	for (int i=0; i<6; ++i)
	{
		for (int j=0; j<6; ++j)
		{
			fg[i] -= fabs(lambda[j])*alpha[j]*eigvec[j][i];
		}
		fg[i] = (fg[i] + fl[i] + fr[i])/2.;
	}

	return fg;
}