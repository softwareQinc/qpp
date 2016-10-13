/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2017 Vlad Gheorghiu (vgheorgh@gmail.com)
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
* \file macros.h
* \brief Preprocessor macros
*/

#ifndef MACROS_H_
#define MACROS_H_

// Useful macros
#ifndef NDEBUG // activate only in DEBUG version
/*! Prints a message */
#define PRINT(x)    std::cout << (x)
/*! Prints a message and adds a new line */
#define PRINTLN(x)  std::cout << (x) << std::endl
/*! Prints an error message to std::cerr */
#define ERROR(x)    std::cerr << (x)
/*! Prints an error message to std::cerr and adds a new line */
#define ERRORLN(x)  std::cerr << (x) << std::endl
#else // deactivate in release version
/*! Prints a message */
#define PRINT(x)
/*! Prints a message and adds a new line */
#define PRINTLN(x)
/*! Prints an error message to std::cerr */
#define ERROR(x)
/*! Prints an error message to std::cerr and adds a new line */
#define ERRORLN(x)
#endif

#endif /* MACROS_H_ */
