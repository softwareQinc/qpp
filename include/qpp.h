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

#include "channels.h"
#include "constants.h"
#include "entanglement.h"
#include "entropies.h"
#include "functions.h"
#include "io.h"
#include "random.h"
#include "types.h"
#include "classes/exception.h"
#include "classes/gates.h"
#include "classes/qudit.h"
#include "classes/randevs.h"
#include "classes/stat.h"
#include "classes/states.h"
#include "classes/timer.h"

namespace qpp
{

// initialize the random devices singleton
RandomDevices& rdevs = RandomDevices::getInstance();

// initialize the gates singleton
const Gates& gt = Gates::getInstance();

// initialize the states singleton
const States& st = States::getInstance();

} /* namespace qpp */

#endif	/* QPP_H_ */

