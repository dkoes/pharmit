/*
 * MinimizationSupport.h
 *
 *  Created on: Oct 30, 2018
 *      Author: dkoes
 */

#ifndef MINIMIZATIONSUPPORT_H_
#define MINIMIZATIONSUPPORT_H_


//pull in the desired library for supporting minimization of compounds
#if USE_GNINA
#include "GninaConverter.h"
namespace MinimizeConverter = GninaConverter;
#else
#include "SminaConverter.h"
namespace MinimizeConverter = SminaConverter;
#endif


#endif /* MINIMIZATIONSUPPORT_H_ */
