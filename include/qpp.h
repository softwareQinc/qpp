/* 
 * File:   qpp.h
 * Author: vlad
 *
 * Created on December 12, 2013, 10:42 PM
 */

#ifndef QPP_H_
#define	QPP_H_

#include "types.h"
#include "util.h"
#include "constants.h"
#include "gates.h"
#include "stat.h"
#include "functional.h"
#include "random.h"
#include "entropy.h"
#include "io.h"
#include "timer.h"

// These are internal functions the user should not have access to
//#include "internal.h"

// MATLAB interface should not be included by default
//#include "matlab.h"

namespace qpp
{

int _init(); // initialization

}

#endif	/* QPP_H_ */

