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
