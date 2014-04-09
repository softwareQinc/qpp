/*
 * measurement.h
 *
 *  Created on: Apr 9, 2014
 *      Author: vlad
 */

#ifndef MEASUREMENT_H_
#define MEASUREMENT_H_

#include <vector>
#include "types.h"
#include "qudit.h"

namespace qpp
{

class Measurement
{
public:
	Measurement(const Qudit& q, const const types::cmat& U){}; // measure in a basis specified by U
	Measurement(const Qudit& q, const std::vector<types::cmat>& Ks); // Kraus measurement
	virtual ~Measurement();
};

} /* namespace qpp */

#endif /* MEASUREMENT_H_ */
