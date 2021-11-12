
#include <stdexcept>
#include <vector>
#include <string>

#include "Environment.h++"
#include "Distance.h++"

#ifndef ERROR_CLASSES_H
#define ERROR_CLASSES_H

/**********************************/
/******        Forest        ******/
/**********************************/
class Except_Forest : public std::exception
{
	public:
		Except_Forest(unsigned int const freqSave, unsigned int const nIter, unsigned int const dimLandscape, bool const overPopulated);
		Except_Forest(unsigned int const nRow, unsigned int const nCol);
		Except_Forest(std::string const& path, std::string const& key);
		const char* what() const throw();

	private:
		std::string m_error_msg = "Error from Forest: ";
		unsigned int m_freqSave;
		unsigned int m_dimLandscape;
		bool m_overPopulated;
};

/*********************************/
/******        Patch        ******/
/*********************************/
class Except_Patch : public std::exception
{
	public:
		Except_Patch(unsigned int const patch_id, std::vector<std::string> const& speciesNames);
		Except_Patch(Distance const& distance, std::string const& speciesNames, unsigned int const target_id, unsigned int const source_id);
		const char* what() const throw();

	private:
		std::string m_error_msg = "Error from Patch: ";
};

/**************************************/
/******        Population        ******/
/**************************************/
class Except_Population : public std::exception
{
	public:
		Except_Population(int const s_inf, int const tallestTree, std::string const& filename);
		Except_Population(int const maxCohorts, std::string const& filename);
		Except_Population(int const maxCohorts, int const nbCohorts, double const t);
		Except_Population(int const maxCohorts, int const nbCohorts);
		Except_Population(unsigned int const iter, std::vector<bool> const& merged_deleted, int const patch_id);
		const char* what() const throw();

	private:
		int m_maxCohorts;
		int m_s_inf;
		std::string m_error_msg = "Error from Population: ";
};

/***********************************/
/******        Species        ******/
/***********************************/
class Except_Species : public std::exception
{
	public:
		Except_Species(int const nbVals, double const integVal);
		const char* what() const throw();

	private:
		std::string m_error_msg = "Error from Species: ";
};

/*************************************/
/******        Landscape        ******/
/*************************************/
class Except_Landscape : public std::exception
{
	public:
		Except_Landscape(int const dim, int const i);
		Except_Landscape(int const dim, std::string const& filename);
		Except_Landscape(double const plotArea, double const deltaLon, double const deltaLat, std::string const climateFile);
		const char* what() const throw();

	private:
		std::string m_error_msg = "Error from Forest when creating the landscape: ";
};

/*************************************/
/******        Dispersal        ******/
/*************************************/
class Except_Dispersal : public std::exception
{
	public:
		Except_Dispersal(double const totalIntegral, std::string const species);
		const char* what() const throw();

	private:
		std::string m_error_msg = "Error from Dispersal: ";
};

#endif
