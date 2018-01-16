/*
 * This file is part of Quantum++.
 *
 * MIT License
 *
 * Copyright (c) 2013 - 2018 Vlad Gheorghiu (vgheorgh@gmail.com)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/**
* \file classes/idisplay.h
* \brief Display interface via the non-virtual interface (NVI)
*/

#ifndef CLASSES_IDISPLAY_H_
#define CLASSES_IDISPLAY_H_

namespace qpp {
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
class IDisplay {
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
    friend inline std::ostream& operator<<(std::ostream& os,
                                           const IDisplay& rhs) {
        return rhs.display(os);
    }
}; /* class IDisplay */

} /* namespace qpp */

#endif /* CLASSES_IDISPLAY_H_ */
