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
namespace internal
{
// macros for defining Singletons, just implement
// private ctors and dtors in defined classes

#define CLASS_SINGLETON(Foo)\
class Foo: public Singleton<Foo>\
{\
    friend class Singleton<Foo>;

#define CLASS_CONST_SINGLETON(Foo)\
class Foo: public Singleton<const Foo>\
{\
    friend class Singleton<const Foo>;

// Singleton policy class

template<typename T>
class Singleton
{
protected:
	Singleton() = default;
	virtual ~Singleton()
	{
	}
	;
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
};

} /* namespace internal */
} /* namespace qpp */

#endif /* SINGLETON_H_ */
