/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2015 Vlad Gheorghiu (vgheorgh@gmail.com)
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
* \file internal/classes/singleton.h
* \brief Singleton pattern via CRTP
*/

#ifndef INTERNAL_CLASSES_SINGLETON_H_
#define INTERNAL_CLASSES_SINGLETON_H_

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
* member function qpp::internal::Singleton::get_instance()
* (qpp::internal::Singleton::get_thread_local_instance()), which returns a
* reference (thread_local reference) to your newly created singleton
* (thread-safe in C++11).
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
* thread_local MySingleton& tls = MySingleton::get_thread_local_instance();
* // Get a thread_local instance
*
* \endcode
*
* \see Code of qpp::Codes, qpp::Gates, qpp::Init, qpp::RandomDevices,
* qpp::States or qpp.h for real world examples of usage.
*/
template<typename T>
class Singleton
{
protected:
    // prevents deleting pointer to instance
    static void operator delete(void*)
    {
    }

    Singleton() noexcept = default;

    Singleton(const Singleton&) = delete;

    Singleton& operator=(const Singleton&) = delete;

    virtual ~Singleton() = default; // to silence base class Singleton<T> has a
    // non-virtual destructor [-Weffc++]

public:
    static T& get_instance() noexcept(std::is_nothrow_constructible<T>::value)
    {
        // Guaranteed to be destroyed.
        // Instantiated on first use.
        // Thread safe in C++11
        static T instance;

        return instance;
    }

    static thread_local T& get_thread_local_instance()
    noexcept(std::is_nothrow_constructible<T>::value)
    {
        // Guaranteed to be destroyed.
        // Instantiated on first use.
        // Thread safe in C++11
        static thread_local T instance;

        return instance;
    }
}; /* class Singleton */

} /* namespace internal */
} /* namespace qpp */

#endif /* INTERNAL_CLASSES_SINGLETON_H_ */
