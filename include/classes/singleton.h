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

// Singleton policy class, use as
//
// class Foo: public Singleton</*const*/ Foo> // use const for const reference
// {
//	   friend class Singleton<Foo>;
//     Foo(/* parameters */){/* implement ctor here */}
// public:
//     // public members here
// };
//
// then get the instance as
//
// Foo& single_Foo = Foo::get_instance();
//
// can use pointers to the instance,
// but of course you cannot delete them
//
// Foo* psingle_Foo = & Foo::get_instance();
// delete psingle_Foo; // this doesn't compile, cannot delete

template<typename T>
class Singleton
{
protected:
	Singleton() = default;
	~Singleton() = default;
	Singleton(const Singleton&) = delete;
	Singleton& operator=(const Singleton&) = delete;
public:
	template<typename ... Args> // parameters forwarded to derived constructor
	static T& get_instance(Args ... args)
	{
		// Guaranteed to be destroyed.
		// Instantiated on first use.
		// Thread safe in C++11
		static T instance { std::forward<Args>(args)... };
		return instance;
	}
};

} /* namespace qpp */

#endif /* SINGLETON_H_ */
