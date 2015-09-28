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
 * ShapeObj.h
 *
 *  A shape object is initialized by a conformer of a molecule.  It provides the
 *  appropriate interface to GSSTreeCreator to construct a voxelized representation
 *  of a molecule, but stores as its object he position of the molecule and its
 *  molecular weight/rotatable bonds for faster filtering.
 *  Created on: Jul 28, 2015
 *      Author: dkoes
 */

#ifndef SHAPEOBJ_H_
#define SHAPEOBJ_H_

#include "OBMoleculeAnalytic.h"
#include "ThreePointData.h"

#define PHARMIT_RESOLUTION (0.5)
#define PHARMIT_DIMENSION (32)


class ShapeObj: public OBAMolecule
{
	static pair<double,double> calcSplitMoments(const vector<Eigen::Vector3d>& coords, unsigned C);

	public:
	struct MolInfo
	{

		unsigned long molPos : TPD_MOLDATA_BITS; //40 bits
		unsigned long pharmPos : TPD_MOLDATA_BITS; //40 bits - pharmacophore information
		unsigned weight : WEIGHT_BITS; //10 bits
		unsigned nrot : ROTATABLE_BITS; //4 bits
		unsigned extra :96 - (2*TPD_MOLDATA_BITS + WEIGHT_BITS + ROTATABLE_BITS);

		MolInfo() :
				molPos(0), pharmPos(0), weight(0), nrot(0), extra(0)
		{
		}
		MolInfo(OpenBabel::OBMol& m, unsigned long mPos) :
				molPos(mPos), pharmPos(0), extra(0)
		{
			weight = ThreePointData::reduceWeight(m.GetMolWt());
			nrot = ThreePointData::reduceRotatable(countRotatableBonds(m));
		}
	}__attribute__((__packed__));

	MolInfo minfo;

	//initializes shape to mol without any transformation
	ShapeObj(OBMol& mol, const MolInfo& info, float dimension,
			float resolution);

	ShapeObj(OBMol& mol, const Eigen::Vector3d translate, const Eigen::Matrix3d& rotate, const MolInfo& info, float dimension,
			float resolution);

	virtual ~ShapeObj()
	{
	}

	virtual void write(ostream& out) const; //writes molinfo

	static void normalizeMol(OBMol& mol);

	static void computeAndApplyNormalization(vector<Eigen::Vector3d>& coords, Eigen::Vector3d& translate, Eigen::Matrix3d& rotate);
};

#endif /* SHAPEOBJ_H_ */
