/*
Pharmit
Copyright (c) David Ryan Koes, University of Pittsburgh and contributors.
All rights reserved.

Pharmit is licensed under both the BSD 3-clause license and the GNU
Public License version 2. Any use of the code that retains its reliance
on the GPL-licensed OpenBabel library is subject to the terms of the GPL2.

Use of the Pharmit code independently of OpenBabel (or any other
GPL2 licensed software) may choose between the BSD or GPL licenses.

See the LICENSE file provided with the distribution for more information.

*/
/*
 * ShapeDistance.h
 *
 *  Created on: Oct 17, 2011
 *      Author: dkoes
 */

#ifndef SHAPEDISTANCE_H_
#define SHAPEDISTANCE_H_
#include "MappableOctTree.h"

typedef
float (*shapeMetricFn)(const MappableOctTree* leftMIV, const MappableOctTree* leftMSV,
		const MappableOctTree* rightMIV, const MappableOctTree* rightMSV);

//a global function for performing shape comparisons between MIV/MSV pairs
extern shapeMetricFn shapeDistance;

enum DistanceFunction {RelativeVolume, AbsVolume, Hausdorff,RelativeTriple, AbsoluteTriple,IncludeExclude,RelVolExclude};

//set shapedistance to the requested function, hausdoff needs the dim
extern void setDistance(DistanceFunction df);

float searchVolumeDist(const MappableOctTree* obj, const MappableOctTree* MIV,
		const MappableOctTree* MSV, float& min, float& max);
float volumeDist(const MappableOctTree* x, const MappableOctTree* y);

#endif /* SHAPEDISTANCE_H_ */
