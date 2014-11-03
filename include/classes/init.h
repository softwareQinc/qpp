/*
 * init.h
 *
 *  Created on: Nov 2, 2014
 *      Author: vlad
 */

#ifndef INCLUDE_CLASSES_INIT_H_
#define INCLUDE_CLASSES_INIT_H_

namespace qpp
{

/**
 * \class qpp::Init
 * \brief const Singleton class that performs
 * additional initializations/cleanups
 */
class Init: public internal::Singleton<const Init> // const Singleton
{
	friend class internal::Singleton<const Init>;
public:
	/**
	 * \brief Additional initializations
	 */
	Init()
	{
		// on entry message
		std::cout << ">>> " << "Starting quantum++..." << std::endl;

		// gets and displays current system time
		time_t current_date;
		time(&current_date);
		std::cout << ">>> " << std::ctime(&current_date);

		// set default output format and precision
		std::cout << std::fixed; // use fixed format for nice formatting
		// std::cout << std::scientific;
		std::cout << std::setprecision(4); // only for fixed or scientific modes
	}
private:
	/**
	 * \brief Cleanups
	 */
	~Init()
	{
		// on exit message
		std::cout << std::endl << ">>> " << "Exiting quantum++..." << std::endl;

		// gets and displays current system time
		time_t current_date;
		time(&current_date);
		std::cout << ">>> " << std::ctime(&current_date);
	}
};
/* class Init */

} /* namespace qpp */

#endif /* INCLUDE_CLASSES_INIT_H_ */
