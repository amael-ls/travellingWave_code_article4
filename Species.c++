
#ifndef SPECIES_C
#define SPECIES_C

// My headers
#include "Species.h++"
#include "Params.h++"

/****************************************/
/******        Constructors        ******/
/****************************************/
Species::Species(std::string const& species_filename, std::string const& species_path, const std::string& delim)
{
	// Read species' main file from the parameters got at the execution
	std::string file_sp = species_path + species_filename;
	par::Params speciesParams(file_sp.c_str(), delim);

	// Species name
	m_speciesName = speciesParams.get_val<std::string>("species");

	// Create the files' name for a given species
	std::string file_G = species_path + m_speciesName;
	file_G.append("_G.txt");

	std::string file_M = species_path + m_speciesName;
	file_M.append("_M.txt");

	std::string file_allometries = species_path + m_speciesName;
	file_allometries.append("_allometries.txt");

	std::string file_scaling_G = species_path + m_speciesName;
	file_scaling_G.append("_scaling_G.txt");

	std::string file_scaling_M = species_path + m_speciesName;
	file_scaling_M.append("_scaling_M.txt");

	std::string file_scaling_dispersal = species_path + m_speciesName;
	file_scaling_dispersal.append("_dispersal.txt");

	// Load parameters from files
	par::Params speciesParams_G(file_G.c_str(), delim);
	par::Params speciesParams_M(file_M.c_str(), delim);
	par::Params speciesParams_allometries(file_allometries.c_str(), delim);
	par::Params speciesParams_scaling_G(file_scaling_G.c_str(), delim);
	par::Params speciesParams_scaling_M(file_scaling_M.c_str(), delim);
	par::Params speciesParams_dispersal(file_scaling_dispersal.c_str(), delim);

	// Allometry parameters (height from dbh)
	a = speciesParams_allometries.get_val<double>("a");
	b = speciesParams_allometries.get_val<double>("b");

	// Allometry parameters (crown area from dbh)
	T_param = speciesParams_allometries.get_val<double>("T_param");
	R0_C0 = speciesParams_allometries.get_val<double>("R0_C0");
	R0_C1 = speciesParams_allometries.get_val<double>("R0_C1");
	R40_C0 = speciesParams_allometries.get_val<double>("R40_C0");
	R40_C1 = speciesParams_allometries.get_val<double>("R40_C1");
	M_C0 = speciesParams_allometries.get_val<double>("M_C0");
	M_C1 = speciesParams_allometries.get_val<double>("M_C1");
	B_C0 = speciesParams_allometries.get_val<double>("B_C0");
	B_C1 = speciesParams_allometries.get_val<double>("B_C1");

	// Fecundity
	// --- Get provided parameters
	fecundity = speciesParams.get_val<double>("fecundity");
	minAgeReproduction = speciesParams.get_val<double>("minAgeReproduction");

	// --- Compute minimal size to reproduce
	// ...... Define solver parameters
	std::string timeSpan = "[0, " + std::to_string(minAgeReproduction) + "]";
	alglib::real_1d_array timeSpanArray(timeSpan.c_str());
	alglib::real_1d_array y = "[0]"; // Initial condition
	double eps = 0.00001;
	double initialTimeStep = 0;
	alglib::odesolverstate s;
	alglib::ae_int_t m;
	alglib::real_1d_array xtbl;
	alglib::real_2d_array ytbl;
	alglib::odesolverreport rep;
	alglib::odesolverrkck(y, timeSpanArray, eps, initialTimeStep, s);

	// ...... Solve ds/dt = G(s, t)
	alglib::odesolversolve(s, Species::growth_callback);
	odesolverresults(s, m, xtbl, ytbl, rep);

	// ...... Set minimal size to reproduce
	if (*(ytbl[1]) < 0 || m != 2)
		throw Except_Species(m, *(ytbl[1]));
	
	minDbhReproduction = *(ytbl[1]);
	minHeightReproduction = 6;

	// Dispersal parameters, not necessarily all provided by the user. Note: I could have used a map to optimise the coding...
	keysToRead = speciesParams_dispersal.get_val<std::string>("keysToRead");
	min_dispersalProba = false;
	max_dispersalDist = false;
	if (keysToRead.find("dispersalProbaThreshold") != std::string::npos)
	{
		dispersalProbaThreshold = speciesParams_dispersal.get_val<double>("dispersalProbaThreshold");
		min_dispersalProba = true;
		std::cout << "Using min dispersal probability" << std::endl;
	}

	if (keysToRead.find("refKernel_doi") != std::string::npos)
		refKernel_doi = speciesParams_dispersal.get_val<std::string>("refKernel_doi");

	if (keysToRead.find("propLDD") != std::string::npos)
		propLDD = speciesParams_dispersal.get_val<double>("propLDD");

	if (keysToRead.find("relLDDtoSDD") != std::string::npos)
		relLDDtoSDD = speciesParams_dispersal.get_val<double>("relLDDtoSDD");

	if (keysToRead.find("dispersalDistThreshold") != std::string::npos)
	{
		dispersalDistThreshold = speciesParams_dispersal.get_val<double>("dispersalDistThreshold");
		max_dispersalDist = true;
		std::cout << "Using max dispersal distance" << std::endl;
	}

	if (keysToRead.find("twoDt_a") != std::string::npos)
		twoDt_a = speciesParams_dispersal.get_val<double>("twoDt_a");

	if (keysToRead.find("twoDt_b") != std::string::npos)
		twoDt_b = speciesParams_dispersal.get_val<double>("twoDt_b");

	// Others
	maxDiameter = speciesParams.get_val<double>("maxDiameter");
}

/**************************************/
/******        Demography        ******/
/**************************************/
// Individual growth rate, s is the diameter, s_star is dbh_star (obtained from height_star)
double Species::v(double s, double const s_star) const
{
	/*
		The values are for Acer saccharum (mesic soil) from Purves 2008:
		Predicting and understanding forest dynamics using a simple tractable model.

		Note that beyond smax, trees stop growing. This is a Neumann boundary condition
	*/
	bool belowSmax = s <= maxDiameter;
	if (belowSmax)
	{
		bool cs = s_star <= s; // ? false : true;
		if (cs)
			return 0.307;
		else
			return 0.086;
	}
	else
		return 0;
}

// Individual death rate, s is the diameter, s_star is dbh_star (obtained from height_star)
double Species::d(double s, double const s_star) const
{
	/*
		The values are for Acer saccharum (mesic soil) from Purves 2008:
		Predicting and understanding forest dynamics using a simple tractable model.
	*/
	bool cs = s_star <= s; // ? false : true;
	if (cs)
		return 0.0034;
	else
		return 0.0195;
}

// Differentiate individual growth rate, s is the diameter, s_star is dbh_star (obtained from height_star)
double Species::dv_ds(double s, double const s_star) const
{
	return 0;
}

// Differentiate individual mortality rate, s is the diameter, s_star is dbh_star (obtained from height_star)
double Species::dd_ds(double s, double const s_star) const
{
	return 0;
}

/*************************************/
/******        Dispersal        ******/
/*************************************/
// From Moorcroft & Lewis 2006: Potential role of natural enemies during tree range expansions following climate change.
// Rk: Their kernel has small mistakes that have been corrected here
double Species::K(double const distance) const
{
	double proba = 0;
	if (refKernel_doi == "10.1016/j.jtbi.2005.12.019") // Moorcroft2006
		proba = (1.0 - propLDD)/2.0*exp(-std::abs(distance)) + propLDD*relLDDtoSDD/2.0*exp(-relLDDtoSDD*distance);

	if (refKernel_doi == "laplacian") // Laplace kernel with parameter = 100 (for testing), Cousens2008, p. 82
		proba = 1.0/(2*M_PI*100*100) * std::exp(- distance/100);

	if (refKernel_doi == "10.2307/176541") // 2Dt Clark1999, Boisvert-Marsh2020 (might be 2021)
		proba = twoDt_a/(M_PI*twoDt_b) * std::exp((-twoDt_a - 1)*std::log(1 + distance*distance/twoDt_b));

	if (refKernel_doi == "gaussian") // Gaussian, Cousens2008, p. 82, with a = 30
		proba = 1.0/(M_PI*30.0*30.0) * std::exp(-distance*distance/(30.0*30.0));

	if (refKernel_doi == "gaussian1D") // Gaussian, with a = 30
		proba = 2.0/(30*30) * distance * std::exp(-distance*distance/(30.0*30.0));

	if (refKernel_doi == "dirac")
		proba = (distance == 0) ? 1 : 0;
	
	return proba;
}

double Species::K(double delta_lon, double delta_lat) const
{
	double proba = 0;
	double distance = sqrt(delta_lon*delta_lon + delta_lat*delta_lat);
	if (refKernel_doi == "10.1016/j.jtbi.2005.12.019") // Moorcroft2006
		proba = (1.0 - propLDD)/2.0*exp(-std::abs(distance)) + propLDD*relLDDtoSDD/2.0*exp(-relLDDtoSDD*distance);
	
	if (refKernel_doi == "laplacian") // Laplacian with parameter = 100 (for testing)
		proba = 1.0/(2*M_PI*100*100) * std::exp(- distance/100);

	if (refKernel_doi == "10.2307/176541") // Clark1999, Boisvert-Marsh2020 (might be 2021)
		proba = twoDt_a/(M_PI*twoDt_b) * std::exp((-twoDt_a - 1)*std::log(1 + distance*distance/twoDt_b));

	if (refKernel_doi == "gaussian") // Gaussian, Cousens2008, p. 82, with a = 30
		proba = 1.0/(M_PI*30.0*30.0) * std::exp(-distance*distance/(30.0*30.0));

	if (refKernel_doi == "gaussian1D") // Gaussian, Cousens2008, p. 82, with a = 30
		proba = proba = 2.0/(30.0*30.0) * distance * std::exp(-distance*distance/(30.0*30.0));

	if (refKernel_doi == "dirac")
		proba = (distance == 0) ? 1 : 0;

	return proba;
}

double Species::K(double const longitude1, double const latitude1, double const longitude2, double const latitude2) const
{
	double proba = 0;
	double distance = sqrt((longitude1 - longitude2)*(longitude1 - longitude2) + (latitude1 - latitude2)*(latitude1 - latitude2));
	if (refKernel_doi == "10.1016/j.jtbi.2005.12.019") // Moorcroft2006
		proba = (1.0 - propLDD)/2.0*exp(-std::abs(distance)) + propLDD*relLDDtoSDD/2.0*exp(-relLDDtoSDD*distance);

	if (refKernel_doi == "dirac")
		proba = (distance == 0) ? 1 : 0;
	
	return proba;
}

bool Species::isDirac() const
{
	return (refKernel_doi == "dirac");
}

/************************************/
/******        Overload        ******/
/************************************/
std::ostream &operator<<(std::ostream &os, Species const& species)
{
	os << std::setprecision(5);
	os << species.m_speciesName << std::endl;
	os << std::endl;

	os << "Allometry parameters height = f(dbh):" << std::endl;
	os << species.a << "\t" << species.b << std::endl;
	os << std::endl;

	os << "Allometry parameters crownArea = f(dbh):" << std::endl;
	os << species.T_param << "\t" << species.R0_C0 << "\t" << species.R0_C1 << "\t" << species.R40_C0
	<< "\t" << species.R40_C1 << "\t" << species.M_C0 << "\t" << species.M_C1
	<< "\t" << species.B_C0 << "\t" << species.B_C1 << std::endl;
	os << std::endl;

	os << "Growth parameters (oversotrey and then understorey):" << std::endl;
	os << 0.307 << "\t" << 0.086 << std::endl;
	os << std::endl;

	os << "Mortality parameters (oversotrey and then understorey):" << std::endl;
	os << 0.0034 << "\t" << 0.0195 << std::endl;
	os << std::endl;


	os << "Fecundity parameters:" << std::endl;
	os << species.fecundity << "\t" << species.minHeightReproduction << std::endl;
	os << std::endl;

	os << "Dispersal parameters:" << std::endl;
	if (species.keysToRead.find("dispersalProbaThreshold") != std::string::npos)
		os << "dispersalProbaThreshold: " << species.dispersalProbaThreshold << std::endl;
	if (species.keysToRead.find("refKernel_doi") != std::string::npos)
		os << "refKernel_doi: " << species.refKernel_doi << std::endl;
	if (species.keysToRead.find("propLDD") != std::string::npos)
		os << "propLDD: " << species.propLDD << std::endl;
	if (species.keysToRead.find("relLDDtoSDD") != std::string::npos)
		os << "relLDDtoSDD: " << species.relLDDtoSDD << std::endl;
	if (species.keysToRead.find("dispersalDistThreshold") != std::string::npos)
		os << "dispersalDistThreshold: " << species.dispersalDistThreshold << std::endl;

	os << "Maximal diameter (correspond to a 45m height):" << std::endl;
	os << species.maxDiameter << std::endl;
	os << std::endl;

	return os;
}

// Sort by alphabetical order
bool operator<(Species const& species1, Species const& species2)
{
	return (species1.m_speciesName < species2.m_speciesName);
}

void Species::printName(std::ostream& os) const
{
	os << m_speciesName;
}

// Getter
std::string Species::getName() const
{
	return m_speciesName;
}

void Species::growth_callback(const alglib::real_1d_array &y, double x, alglib::real_1d_array &dy, void *ptr) 
{
	/*
		x = time
		y = function to find dy/dx = growth_callback(x, y(x))

		It is assumed that the tree is in the understorey, i.e., canopy status = false
	*/

	// Growth function
	dy[0] = 0.086;
}

#endif
