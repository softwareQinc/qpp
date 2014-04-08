/* 
 * File:   qpp.h
 * Author: vlad
 *
 * Created on December 12, 2013, 10:42 PM
 */

// IMPORTANT: instantiation of global singletons Gates and RandomDevices
// Any additional singletons should be instantiated here
// ALWAYS include it in main.cpp
#ifndef QPP_H_
#define	QPP_H_

#include "types.h"
#include "constants.h"
#include "gates.h"
#include "stat.h"
#include "functions.h"
#include "random.h"
#include "entropies.h"
#include "io.h"
#include "timer.h"
#include "exception.h"
#include "channels.h"
#include "randevs.h"

namespace qpp
{

// initialize the random devices
RandomDevices *rdevs = RandomDevices::getInstance();

// initialize the gates
const Gates *gt = Gates::getInstance();

} /* namespace qpp */

#endif	/* QPP_H_ */

