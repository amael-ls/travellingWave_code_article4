
/* Description of the class
The class species defines a tree species by its parameters and vital rates.

It contains:
	- all the species specific parameters required to define a species (including allometries)
	- vital rates (functions)

I list the functions here, but describe them in the associated c++ file:
	- Constructor
	- Overloading for ease
	- Vital rates and their differentiate:
		* v (growth speed)
		* d (deat rate)
		* dv_ds differentiate of v with respect to s (size variable)
		* dd_ds differentiate of d with respect to s (size variable)

Allometries:
1. Calculate the height of a tree given its dbh:
	- Parameters for 106 species, source table S3, Purves2007, a and b
	- Function define at the end of the program, source appendix S2, end 1st page

2. Calculate the radius of the crown at a distance y from the top, given dbh of a tree:
	- Parameters for 106 species, source table S3, Purves2007, T
	- Function define ate the end of the program, source appendix S1, eq S1.6, S1.7
	
3. Data C0-C1:
	- Source = table S2 Purves2007

4. Bibliography:
		- Crown plasticity and competition for canopy space: a new spatially implicit model parameterized for 250 North American tree species
*/

#ifndef SPECIES_H
#define SPECIES_H

// Official headers
#include <iomanip> // std::setw, std::left, std::setprecision
#include <string>
#include <cmath> // for log

// Alglib header
#include "diffequations.h" // To compute minSizeReproduction

// My headers
#include "Error_classes.h++"
#include "Environment.h++"

/*************************************/
/******      CLASS SPECIES      ******/
/*************************************/

class Species
{
	// Friendship
	friend class Population;
	friend class Dispersal;
	friend class Forest;
	friend class Cohort;
	friend class Patch;

	public :
		// constructors
		Species(std::string const& species_filename, std::string const& species_path, const std::string& delim);

		// friend functions and overload
		friend std::ostream& operator<<(std::ostream& os, const Species& species);
		friend bool operator<(Species const& species1, Species const& species2);

		// Demographic functions
		double v(double s, double const s_star) const;
		double d(double s, double const s_star) const;
		double dv_ds(double s, double const s_star) const;
		double dd_ds(double s, double const s_star) const;

		// Dispersal
		double K(double const distance) const;
		double K(double delta_lon, double delta_lat) const;
		double K(double const longitude1, double const latitude1, double const longitude2, double const latitude2) const;
		bool isDirac() const;

		// others
		static void growth_callback(const alglib::real_1d_array &y, double x, alglib::real_1d_array &dy, void *ptr);
		void printName(std::ostream& os) const;
		std::string getName() const;

	private :
		// Species' name
		std::string m_speciesName;

		// Allometry parameters (height from dbh)
		double a, b;

		// Allometry parameters (crown area from dbh)
		double T_param, R0_C0, R0_C1, R40_C0, R40_C1, M_C0, M_C1, B_C0, B_C1;

		// Fecundity
		double fecundity;
		double minAgeReproduction;
		double minDbhReproduction;
		double minHeightReproduction;

		// Dispersal parameters, not necessarily all provided by the user
		// --- General parameters
		std::string keysToRead;
		double dispersalProbaThreshold; // Value between 0 and 1. Threshold from which probability K(x, y) is considered 0
		bool min_dispersalProba; // True if dispersalProbaThreshold is defined
		std::string refKernel_doi; // The DOI from which the kernel is from
		double dispersalDistThreshold; // Threshold distance from which probability K(x, y) is considered 0 (i.e., max dist)
		bool max_dispersalDist; // True if dispersalDistThreshold is defined

		// --- Moorcroft 2006
		double relLDDtoSDD; // Average long-distance dispersal relative to short-distance dispersal ---> 1/beta_l = 10 => beta_l = 0.1
		double propLDD; // Fraction of fecundity going to long-distance dispersal ---> p in Moorcroft, 0.1

		// Clark 1999, cf Boisvert-Marsh 2020 for some values
		double twoDt_a;
		double twoDt_b;

		// Others
		double maxDiameter;
};

#endif
