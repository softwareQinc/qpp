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

template<typename T> // Singleton policy class
class Singleton
{
protected:
	Singleton()=default;
	Singleton(const Singleton&) = delete;
	Singleton& operator=(const Singleton&) = delete;
	~Singleton()=default;
public:
	static T& getInstance() // singleton
	{
		// Guaranteed to be destroyed.
		// Instantiated on first use.
		// Thread safe in C++11
		static T instance;
		return instance;
	}
};

} /* namespace qpp */

#endif /* SINGLETON_H_ */
