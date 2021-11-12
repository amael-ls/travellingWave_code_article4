
#ifndef DISTANCE_H
#define DISTANCE_H

// Official headers
#include <iostream>
#include <vector>
#include <limits> // for std::numeric_limits
#include <cmath>

// My headers
#include "Environment.h++"

class Distance
{
	friend class Dispersal;
	public:
	// Constructors, latitude = row x deltaLat, longitude = col x deltaLon
		Distance(int const row1, int const col1, int const row2, int const col2, double const deltaLat, double deltaLon);
		Distance(Environment const& env1, Environment const& env2, double const deltaLat, double const deltaLon);

	// Overloads
		friend bool operator<(Distance const d1, Distance const d2);
		friend bool operator==(Distance const d1, Distance const d2);
		friend std::ostream& operator<<(std::ostream& os, Distance const &dist);

	// Other
		std::string to_string() const;

	private:
		// double m_orthodromic;
		double m_euclidean;
		std::vector<int> m_manhattan;
};

#endif
