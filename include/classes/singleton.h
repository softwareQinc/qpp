/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2014 Vlad Gheorghiu (vgheorgh@gmail.com)
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

#ifndef INCLUDE_CLASSES_SINGLETON_H_
#define INCLUDE_CLASSES_SINGLETON_H_

namespace qpp
{
namespace internal // internal class, do not modify
{

/**
* \class qpp::internal::Singleton
* \brief Singleton policy class, used internally to implement
* the singleton pattern via CRTP (Curiously recurring template pattern)
*
* To implement a singleton, derive your class from qpp::internal::Singleton,
* make qpp::internal::Singleton a friend of your class, then declare
* the constructor of your class as private. To get an instance, use the static
* member function qpp::internal::Singleton::get_instance(), which returns a
* reference to your newly created singleton (thread-safe in C++11).
*
* Example:
* \code
* class MySingleton: public qpp::internal::Singleton<MySingleton>
* {
* 	   friend class qpp::internal::Singleton<MySingleton>;
* public:
*     // Declare all public members here
* private:
*     MySingleton()
*     {
*         // Implement the constructor here
*     }
* };
*
* MySingleton& mySingleton = MySingleton::get_instance(); // Get an instance
*
* \endcode
*
* \see Code of qpp::Codes, qpp::Gates, qpp::RandomDevices,
* qpp::States or qpp.h for real world examples of usage.
*/
template<typename T>
class Singleton
{
protected:
    Singleton() = default;

    virtual ~Singleton()
    {
    }

    // = default yields "looser throw specifer in g++ <= 4.7"
    // see http://stackoverflow.com
    // /questions/11497252/default-destructor-nothrow
    Singleton(const Singleton &) = delete;

    Singleton &operator=(const Singleton &) = delete;

public:
    static T &get_instance()
    {
        // Guaranteed to be destroyed.
        // Instantiated on first use.
        // Thread safe in C++11
        static T instance;
        return instance;
    }
}; /* class Singleton */

} /* namespace internal */
} /* namespace qpp */

#endif /* INCLUDE_CLASSES_SINGLETON_H_ */
