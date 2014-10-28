/*
 * Singleton.h
 *
 *  Created on: May 20, 2014
 *      Author: vlad
 */

#ifndef SINGLETON_H_
#define SINGLETON_H_

namespace qpp
{
namespace internal // internal class, do not modify
{

// Singleton policy class
/**
 * \class qpp::internal::Singleton
 * \brief Singleton policy class, used internally to implement
 * the singleton pattern via CRTP (Curiously recurring template pattern)
 *
 * To implement a singleton, derive your class from \a qpp::internal::Singleton,
 * make \a qpp::internal::Singleton a friend of your class, then declare
 * the constructor of your class as private. To get an instance, use the static
 * member function \a qpp::internal::Singleton::get_instance(), which returns a
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
 * \see Code of \a qpp::Gates, \a qpp::RandomDevices, \a qpp::States or
 * \a qpp.h for real world examples of usage.
 */
template<typename T>
class Singleton
{
protected:
	Singleton() = default;
	virtual ~Singleton(){};
	Singleton(const Singleton&) = delete;
	Singleton& operator=(const Singleton&) = delete;
public:
	static T& get_instance()
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

#endif /* SINGLETON_H_ */
