#include "Cow/src/Array.hpp"
#include "Reconstruction.hpp"




class UserParameters
{
public:
	int GridDimensions[3] = {1, 1, 1};
	int NumberOfGhostZones = 3;
	int GridLength[3] = {1, 1, 1};

	double CFLSafetyNumber = 0.5;

	double DataDumpInterval = 1.;
	double StopTime = -1;
	int StopCycle = -1;

	int ReconstructionMethod = 1;

	UserParameters();
	void ReadParameters();
	void ErrorCheck();
};


class SolutionData
{
private:
	Cow::Array BField;
	Cow::Array EField;
	Cow::Array PField;
	Cow::Array x;

	// Number of grids in each direction (including ghost zone)
	int N[3] = {1,1,1};
	// Number of grids in each direction (no ghost zone)
	int Nint[3] = {1,1,1};

	// start and end index in each dimension (exclude ghost zones)
	int StartIndex[3] = {0, 0, 0};
	int EndIndex[3] = {1, 1, 1};

	double dx[3] = {1, 1, 1};

public:
	SolutionData(UserParameters userParameters);

	Cow::Array& GetBField();
	Cow::Array& GetEField();
	Cow::Array& GetPField();
	void CombineBEtoP();

	double TotalEnergy;
	double TotalEdotB;
	double MaxdivB;
	
	void ComputeTotalEnergy();
	void ComputeTotalEdotB();
	void ComputeMaxdivB();


	void InitialData();

	// resource allocation is initialization
	// heap
};




class Solver
{
private:
	/* auxilary variables */
	//Cow::Array FluxL;
	//Cow::Array FluxR;
	Cow::Array FluxG;
	Cow::Array P0;
	Cow::Array DeltaP;
	Cow::Array Psi;
	Cow::Array Psi0;

	Cow::Array DeltaPsi;

	// Number of grids in each direction (including ghost zone)
	int N[3] = {1,1,1};

	// start and end index in each dimension (exclude ghost zones)
	int StartIndex[3] = {0, 0, 0};
	int EndIndex[3] = {1, 1, 1};

	double dx[3] = {1, 1, 1};

	double dt = 1.;

	int ReconstructionMethod = 1;

	Reconstruction reconstruct;
	Reconstruction::Operation reconstructModeC2R;
	Reconstruction::Operation reconstructModeC2L;


public:

	Solver(UserParameters userParameters);
	void SetReconstructionMethod();

	void ComputeDeltaP(const Cow::Array& P);
	void AddSourceTerm(const Cow::Array& P);
	void SweepX(const Cow::Array& P);
	void SweepY(const Cow::Array& P);
	void SweepZ(const Cow::Array& P);
	void CleanEdotB(Cow::Array& P);

	void Advance(Cow::Array& P);
	void Cache(const Cow::Array& P, const Cow::Array& Psi);
	void TVDstep(Cow::Array& P, const Cow::Array& P0, const Cow::Array& DeltaP, double alpha, double beta);
	void SetBoundaryValue(Cow::Array& P);
	void Setdt(double dt);
	void ShrinkE(Cow::Array& P);

	void DednerDamp(const Cow::Array& P);
};

class SimulationControl
{
private:
	double dx[3] = {1, 1, 1};
	double dtFixed = 1.;

	double NextDataDumpTime = 0;
	double DataDumpInterval = 0;

	double StopTime = 1;
	int StopCycle = 100;

public:
	SimulationControl(UserParameters userParameters);

	double SimTime = 0;
	int Cycle = 0;
	int OutputCheckPoint = 0;

	double Getdt() const;
	bool ShouldOutput() const;
	bool ShouldContinue() const;

	void UpdateTimeandCycle(double dt);
	void UpdateNextDumpTime();
	void UpdateCheckPoint();

	void OutputData(Cow::Array& B, Cow::Array& E) const;

};

