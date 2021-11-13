
#ifndef DISPERSAL_C
#define DISPERSAL_C

//* Maybe a good opportunity to cease to switch from ALGLIB to Boost!
//? https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/quadrature.html

// My headers
#include "Dispersal.h++"

/****************************************/
/******        Constructors        ******/
/****************************************/
Dispersal::Dispersal(Species const* const sp, std::string const climateFilename):
	m_species(sp), m_min_dispersalProba(sp->min_dispersalProba),
	m_max_dispersalDist(sp->max_dispersalDist), m_totalIntegral(0)
{
	/**** Read landscape parameters from file climateFilename ****/
	par::Params climateParams(climateFilename.c_str(), "=");

	m_nRow_land = climateParams.get_val<unsigned int>("nRow");
	m_nCol_land = climateParams.get_val<unsigned int>("nCol");
	m_dim_land = m_nRow_land*m_nCol_land;

	m_deltaLon = climateParams.get_val<double>("deltaLon");
	m_deltaLat = climateParams.get_val<double>("deltaLat");

	if (m_min_dispersalProba)
		m_dispersalProbaThreshold = sp->dispersalProbaThreshold;
	
	if (m_max_dispersalDist)
		m_dispersalDistThreshold = sp->dispersalDistThreshold;
}

/*************************************************/
/******        Dispersal integration        ******/
/*************************************************/
void Dispersal::kernel(double r, double xminusa, double bminusx, double &y, void *ptr) const
{
	y = m_species->K(r);
}

void Dispersal::wrapper_r_integral(double x, double xminusa, double bminusx, double &y, void *ptr)
{
	// explicitly cast global variable <pt2Object> to a pointer to Dispersal
	// warning: <pt2Object> MUST point to an appropriate object!
	Dispersal* mySelf = (Dispersal*) pt2Object;

	// call member
	mySelf->kernel(x, xminusa, bminusx, y, ptr);
}

void landscapeIntegrals(Dispersal& disp)
{
	// Integral bounds
	// --- Cartesian
	double y1, y2;

	// Dimension patch (interval y1 ---> y2, with y2 = y1 + deltaLat)
	double deltaLon(disp.m_deltaLon);
	double deltaLat(disp.m_deltaLat);
	double euclideanDistanceToZero;

	// Values for integral
	pt2Object = (void*) &disp;

	double integPatch(0);
	disp.m_totalIntegral = 0;

	alglib::autogkstate ss;
	alglib::autogkreport reprep;

	for (unsigned int row = 0; row < disp.m_nRow_land; ++row) // latitude direction
	{
		// Reset integPatch for current patch
		integPatch = 0;

		// Coordinates and distance, row = y = latitude
		y1 = row*deltaLat;
		Distance distanceToZero(0, 0, row, 0, deltaLat, deltaLon);
		euclideanDistanceToZero = sqrt(y1*y1);

		y2 = y1 + deltaLat;

		alglib::autogksmooth(y1, y2, ss);
		alglib::autogkintegrate(ss, Dispersal::wrapper_r_integral);
		alglib::autogkresults(ss, integPatch, reprep);
		
		disp.m_totalIntegral += integPatch;

		if ((disp.m_min_dispersalProba) &&(integPatch >= disp.m_dispersalProbaThreshold)) // If using proba threshold
		{
			if (disp.m_map_distance_integral.find(distanceToZero) == disp.m_map_distance_integral.end())
				disp.m_map_distance_integral[distanceToZero] = integPatch;
		}

		if ((disp.m_max_dispersalDist) && (euclideanDistanceToZero <= disp.m_dispersalDistThreshold)) // If using distance threshold
		{
			if (disp.m_map_distance_integral.find(distanceToZero) == disp.m_map_distance_integral.end())
				disp.m_map_distance_integral[distanceToZero] = integPatch;
		}
	}

	// Because we integrated only on half the line (from 0 to +Inf), we need to multiply the total integral by 2
	// disp.m_totalIntegral *= 2; // If the landscape is big enough, this should be close to 1

	std::map<Distance, double>::const_iterator it = disp.m_map_distance_integral.cbegin();
	for (; it != disp.m_map_distance_integral.cend(); ++it)
		std::cout << it->first << "    " << it->second << std::endl;
	std::cout << std::endl;

	if ((disp.m_totalIntegral > 1.01) || (disp.m_totalIntegral < 0))
		throw Except_Dispersal(disp.m_totalIntegral, disp.m_species->getName());
}

/************************************/
/******        Overload        ******/
/************************************/
std::ostream& operator<<(std::ostream& os, Dispersal const &dispersal)
{
	// Landscape
	os << std::endl << "Dimensions (row x col): " << dispersal.m_nRow_land << " x " << dispersal.m_nCol_land << std::endl;
	os << "Resolution (lon x lat): " << dispersal.m_deltaLon << " x " << dispersal.m_deltaLat << std::endl;
	
	os << "Integral on Γ: " << dispersal.m_totalIntegral << std::endl;

	os << "Dimension of the dispersal map: " << (dispersal.m_map_distance_integral).size();
	return os;
}

#endif
