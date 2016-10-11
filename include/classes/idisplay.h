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
* \file classes/idisplay.h
* \brief Display interface via the non-virtual interface (NVI)
*/

#ifndef CLASSES_IDISPLAY_H_
#define CLASSES_IDISPLAY_H_

namespace qpp
{

/**
* \class qpp::IDisplay
* \brief Abstract class (interface) that mandates the definition of
* virtual std::ostream& display(std::ostream& os) const
*
* This class defines friend inline std::ostream& operator<<
* (std::ostream& os, const qpp::IDisplay& rhs). The latter delegates
* the work to the pure private virtual function qpp::IDisplay::display()
* which has to be overridden by all derived classes.
*/
class IDisplay
{
private:
    /**
    * \brief Must be overridden by all derived classes
    *
    * The actual stream extraction processing is performed by the overriden
    * member function in the derived class. This function is automatically
    * invoked by friend inline std::ostream& operator<<(std::ostream& os,
    * const IDisplay& rhs).
    */
    virtual std::ostream& display(std::ostream& os) const = 0;

public:
    /**
    * \brief Default constructor
    */
    IDisplay() = default;

    /**
    * \brief Default copy constructor
    */
    IDisplay(const IDisplay&) = default;

    /**
    * \brief Default move constructor
    */
    IDisplay(IDisplay&&) = default;

    /**
    * \brief Default copy assignment operator
    */
    IDisplay& operator=(const IDisplay&) = default;

    /**
    * \brief Default move assignment operator
    */
    IDisplay& operator=(IDisplay&&) = default;

    /**
    * \brief Default virtual destructor
    */
    virtual ~IDisplay() = default;

    /**
    * \brief Overloads the extraction operator
    *
    * Delegates the work to the virtual function qpp::IDisplay::display()
    */
    friend inline
    std::ostream& operator<<(std::ostream& os, const IDisplay& rhs)
    {
        return rhs.display(os);
    }
}; /* class IDisplay */

} /* namespace qpp */

#endif /* CLASSES_IDISPLAY_H_ */
