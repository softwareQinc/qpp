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
//#include "internal.h" // DO NOT include this explicitly in qpp.h. These are internal functions the user should not have access to.

namespace qpp
{

int _init(); // initialization

}

#endif	/* QPP_H_ */

