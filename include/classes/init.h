/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2016 Vlad Gheorghiu (vgheorgh@gmail.com)
 *
 * This file is part of Quantum++.
 *
 * Quantum++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Quantum++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Quantum++.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
* \file classes/init.h
* \brief Initialization
*/

#ifndef CLASSES_INIT_H_
#define CLASSES_INIT_H_

namespace qpp
{

/**
* \class qpp::Init
* \brief const Singleton class that performs
* additional initializations/cleanups
*/
class Init final : public internal::Singleton<const Init> // const Singleton
{
    friend class internal::Singleton<const Init>;

private:
    /**
    * \brief Additional initializations
    */
    Init()
    {
//        // on entry message
//        std::cout << ">>> Starting Quantum++..." << std::endl;
//
//        // gets and displays current system time
//        std::time_t current_date = std::time(nullptr);
//        std::cout << ">>> " << std::ctime(&current_date) << std::endl;
//
//        // set default output format and precision
//        std::cout << std::fixed; // use fixed format for nice formatting
//        // std::cout << std::scientific;
//        std::cout << std::setprecision(4); // only for fixed/scientific modes
    }

    /**
    * \brief Cleanups
    */
    ~Init()
    {
//        // on exit message
//        std::cout << std::endl << ">>> Exiting Quantum++..." << std::endl;
//
//        // gets and displays current system time
//        std::time_t current_date = std::time(nullptr);
//        std::cout << ">>> " << std::ctime(&current_date);
    }
}; /* class Init */

} /* namespace qpp */

#endif /* CLASSES_INIT_H_ */
