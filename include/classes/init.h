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
* \file classes/init.h
* \brief Initialization
*/

#ifndef CLASSES_INIT_H_
#define CLASSES_INIT_H_

namespace qpp {
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
    Init() {
        // set default output format and precision
        //        std::cout << std::fixed; // use fixed format for nice
        //        formatting
        // std::cout << std::scientific;
        //        std::cout << std::setprecision(4); // only for
        //        fixed/scientific modes
    }

    /**
    * \brief Cleanups
    */
    ~Init() {}
}; /* class Init */

} /* namespace qpp */

#endif /* CLASSES_INIT_H_ */
