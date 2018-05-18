#include "../Cow/src/Array.hpp"
#include "../Cow/src/MPI.hpp"
#include "Reconstruction.hpp"

class UserParameters;


class Cart
{
/* initialize MPI and apply boundary condition on a parallel machine */

public:
	Cow::MpiCartComm cart;
  	int StartPos[3] = { 0, 0, 0 };           // Start position of the block
  	int EndPos[3] = { 0, 0, 0 };             // End position of the block
  	int GlobalDim[3] = { 1, 1, 1 };          // Total number of grids in the global box
  	int PatchSize[3] = { 1, 1, 1 };          // Size of the block patch
  	int Coordinates[3] = { 0, 0, 0 };        // Coordinate of the block patch in the lattice
  	int LatticeDim = 1;                      // Dimension of the computer lattice
  	int WorldRank = 0;                       // Rank of the block patch in MPICommWorld

  	int NumGhost;                            // Number of Ghost Zones

  	// Constructor
  	Cart(UserParameters userParameters);
  	// Communitate to neighboring cells to apply periodic boundary condition
  	void apply (Cow::Array &A, int ng, int N[3]);
  	// Read global dimension of box from userparameters
  	void SetGlobalDim(UserParameters userParameters);

};

class UserParameters
{
public:
	char* file;
	int GridDimensions[3] = {1, 1, 1};        // Computation box dimension
	int NumberOfGhostZones = 3;               // Number of ghost zones
	int GridLength[3] = {1, 1, 1};            // Total size, computation domain + ghost zones

	double dx[3] = {1, 1, 1};                 // Grid spacing in each direction

	double CFLSafetyNumber = 0.5;             // Courant safety number

	double DataDumpInterval = 1.;	          // Interval for data dumping
	double StopTime = -1;                     // Stop time
	int StopCycle = -1;                       // Stop cycle

	int ReconstructionMethod = 1;             // Reconstruction method for Riemann solver

	double AmplitudeOfAlfvenPacket = 0.1;     // Amplitude Of Alfven Packet, used for alfven wave collision
	int ntheta = 0;

	int SpongeLayer = 0;           // On/off of sponge layer for each direction

	// Constructor
	UserParameters(char* file);
	// Read parameters from file
	void ReadParameters();
	// Check for parameter errors
	void ErrorCheck();
	// Set computational dimension for each block
	void SetPatchDimension(Cart cart);
};


class SolutionData
{
private:
	// Field variables
	Cow::Array BField;
	Cow::Array EField;
	Cow::Array PField;
	// Coordinates
	Cow::Array x;
	// Resistivity
	Cow::Array Resistivity;
	Cow::Array DampedEnergy;

	// Number of grids in each direction (including ghost zone)
	int N[3] = {1,1,1};
	// Number of grids in each direction (no ghost zone)
	int Nint[3] = {1,1,1};

	// Number dimensions
	int GridDim = 0;

	// start and end index in each dimension (exclude ghost zones)
	int StartIndex[3] = {0, 0, 0};
	int EndIndex[3] = {1, 1, 1};

	// Grid Spacings
	double dx[3] = {1, 1, 1};

	//Amplitude Of Alfven Packet, used for alfven wave collision
	double AmplitudeOfAlfvenPacket = 0.1;
	int ntheta = 0;

public:
	// Constructor
	SolutionData(UserParameters userParameters, Cart cart);

	// On/off of sponge layer for each direction
	int SpongeLayer = 0;

	// Return Arrays of Field Variables
	Cow::Array& GetBField();
	Cow::Array& GetEField();
	Cow::Array& GetPField();
	// Return Arrays of Resistivity
	Cow::Array& GetResistivity();
	Cow::Array GetDampedEnergy();


	// Combine E and B field to the primitive variable P
	void CombineBEtoP();

	double TotalEnergy;
	double MaxEdotB;
	double MaxdivB;
	double MagneticEnergy;
	double ElectricEnergy;
	double TotalOhmHeat;
	
	void ComputeEnergy();
	void ComputeDampedEnergy(double dt);

	void ComputeMaxEdotB();
	void ComputeMaxdivB();

	// Initialize Field Configurations
	void InitialData();
	void TestCases();
	void AlfvenPacket1D();
	void AlfvenPacket2D();
	void AlfvenPacket3D();
	void Random2D();
	void Equilibrium();
	void InitializeResistivity();

};

class Solver
{
private:
    // Gudonov Flux
	Cow::Array FluxG;
	// Field before evolution step
	Cow::Array P0;
	// Change to the Field
	Cow::Array DeltaP;
	// Auxillary field for Dedner damping
	Cow::Array Psi;
	Cow::Array Psi0;
	Cow::Array DeltaPsi;

	// Number of grids in each direction (including ghost zone)
	int GridLength[3] = {1,1,1};
	// Number of guard zones
	int NumGhost = 0;

	// start and end index in each dimension (exclude ghost zones)
	int StartIndex[3] = {0, 0, 0};
	int EndIndex[3] = {1, 1, 1};

	double dx[3] = {1, 1, 1};

	double dt = 1.;

	// On/off of sponge layer for each direction
	int SpongeLayer = 0;

	int ReconstructionMethod = 1;

	Reconstruction reconstruct;
	Reconstruction::Operation reconstructModeC2R;
	Reconstruction::Operation reconstructModeC2L;

public:

	Solver(UserParameters userParameters);
	void SetReconstructionMethod();

	void ComputeDeltaP(const Cow::Array& P, const Cow::Array& Resistivity);
	void AddJpert(const Cow::Array& P);
	void AddJparallel(const Cow::Array& P);
	void InitializeDeltaP(const Cow::Array& P);

	void SweepX(const Cow::Array& P);
	void SweepY(const Cow::Array& P);
	void SweepZ(const Cow::Array& P);
	void CleanEdotB(Cow::Array& P);

	void Advance(Cart cart, Cow::Array& P, const Cow::Array& Resistivity);
	void Cache(const Cow::Array& P, const Cow::Array& Psi);
	void TVDstep(Cow::Array& P, const Cow::Array& P0, const Cow::Array& DeltaP, double alpha, double beta);
	void SetBoundaryValue(Cart cart, Cow::Array& P);
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

	void OutputData(Cow::MpiCommunicator world, Cart cart, SolutionData solution) const;

	std::vector<double> TimeSeries;
	std::vector<double> ElectricEnergy;
	std::vector<double> MagneticEnergy;
	std::vector<double> TotalEnergy;
	std::vector<double> OhmHeat;


	void UpdateTimeSeries(Cow::MpiCommunicator world, SolutionData solution);

};
