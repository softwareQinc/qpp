/*
 * init.h
 *
 *  Created on: Nov 2, 2014
 *      Author: vlad
 */

#ifndef INIT_H_
#define INIT_H_

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
		// On entry message
		std::cout << "Starting qpp..." << std::endl;

		// set default output format and precision
		std::cout << std::fixed; // use fixed format for nice formatting
		std::cout << std::setprecision(4); // only for fixed or scientific modes
	}
private:
	/**
	 * \brief Cleanups
	 */
	~Init()
	{
		// On exit message
		std::cout << std::endl << "Exiting qpp..." << std::endl;
	}
};
/* class Init */

} /* namespace qpp */

#endif /* INIT_H_ */
