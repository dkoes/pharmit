/*
 * ShapeObj.h
 * GPLv2 License.
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

class ShapeObj: public OBAMolecule
{
	static pair<double,double> calcSplitMoments(const vector<Eigen::Vector3d>& coords, unsigned C);

	void normalizeMol(OBMol& mol);

	public:
	struct MolInfo
	{

		unsigned long molPos : TPD_MOLDATA_BITS; //40 bits
		unsigned weight : WEIGHT_BITS; //10 bits
		unsigned nrot : ROTATABLE_BITS; //4 bits
		unsigned extra :64 - (TPD_MOLDATA_BITS + WEIGHT_BITS + ROTATABLE_BITS);

		MolInfo() :
				molPos(0), weight(0), nrot(0), extra(0)
		{
		}
		MolInfo(OpenBabel::OBMol& m, unsigned long mPos) :
				molPos(mPos), extra(0)
		{
			weight = ThreePointData::reduceWeight(m.GetMolWt());
			nrot = ThreePointData::reduceRotatable(countRotatableBonds(m));
		}
	};

	MolInfo minfo;

	ShapeObj(OBMol& mol, const MolInfo& info, float dimension,
			float resolution);

	ShapeObj(OBMol& mol, const Eigen::Vector3d translate, const Eigen::Matrix3d& rotate, const MolInfo& info, float dimension,
			float resolution);

	virtual ~ShapeObj()
	{
	}

	virtual void write(ostream& out) const; //writes molinfo

	static void computeAndApplyNormalization(vector<Eigen::Vector3d>& coords, Eigen::Vector3d& translate, Eigen::Matrix3d& rotate);
};

#endif /* SHAPEOBJ_H_ */
