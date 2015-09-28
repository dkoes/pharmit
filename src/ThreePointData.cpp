/*
 Pharmit
 Copyright (C) 2011  David Ryan Koes and the University of Pittsburgh

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

/*
 * ThreePointData.cpp
 *
 *  Created on: Aug 6, 2010
 *      Author: dkoes
 */

#include "ThreePointData.h"
#include "Triplet.h"
#include "basis.h"
#include "SphereGrid.h"

using namespace boost;
using namespace OpenBabel;

//multplier for float values in reduction
#define REDUCED_FLOAT_SIG (1000.0)
//reduce down to accuracy of a short
short ThreePointData::reduceFloat(double val)
{
	int i = round(val * REDUCED_FLOAT_SIG);
	if (i > SHRT_MAX)
		return SHRT_MAX;
	else if (i < SHRT_MIN)
		return SHRT_MIN;
	return i;
}

//val is assumed to be between -pi and pi
unsigned ThreePointData::reduceAngle(double val)
{
	double v = (val + M_PI) / (2 * M_PI); //between 0 and 1
	return round(v * (1 << TPD_ANGLE_BITS));
}

double ThreePointData::unreduceAngle(int val)
{
	double a = val / (double) (1 << TPD_ANGLE_BITS); //between 0 and 1
	return a * (2 * M_PI) - M_PI;
}

unsigned short ThreePointData::reduceFloatUnsigned(double val)
{
	unsigned i = round(val * REDUCED_FLOAT_SIG);
	if (i > USHRT_MAX)
		return USHRT_MAX;
	return i;
}

unsigned ThreePointData::reduceFloatBigUnsigned(double val)
{
	if (val < 0)
		val = 0;
	unsigned i = round(val * REDUCED_FLOAT_SIG);
	return i;
}

double ThreePointData::unreduceFloat(int val)
{
	double ret = val;
	ret /= REDUCED_FLOAT_SIG;
	return ret;
}

unsigned ThreePointData::reduceWeight(double d)
{
	if (d < 128)
		return 0;
	d -= 128;
	//range up to 1024+128
	unsigned ret = round(d);
	if (ret >= 1024)
		return 1023;
	return ret;
}

double ThreePointData::unreduceWeight(unsigned v)
{
	return v + 128;
}

unsigned ThreePointData::reduceRotatable(unsigned nr)
{
	if (nr >= (1 << ROTATABLE_BITS))
		return (1 << ROTATABLE_BITS) - 1;
	return nr;

}

//compute value representing orientation of average vector relative to coordinate
//system of i,j,k triangle (if possible)
unsigned ThreePointData::vecValue(const PharmaPoint *pt, const PharmaPoint *I,
		const PharmaPoint *J, const PharmaPoint *K)
{
	if (pt->vecs.size() == 0)
		return 0;

	CoordinateBasis basis(*I, *J, *K, true);
	if (!basis.hasValidBasis())
		return 0;

	//take average of vectors, but not for aromatic (since they're opposite)
	vector3 v(0, 0, 0);
	if (pt->pharma->name == "Aromatic" || pt->vecs.size() == 1)
	{
		v = pt->vecs[0];
	}
	else
	{
		for (unsigned i = 0, n = pt->vecs.size(); i < n; i++)
			v += pt->vecs[i];
		v /= pt->vecs.size();
	}

	float x = 0;
	float y = 0;
	float z = 0;
	vector3 zero(0, 0, 0);
	basis.setTranslate(zero);
	basis.replot(v.x(), v.y(), v.z(), x, y, z);
	unsigned ret = sphereGrid.pointToGrid(x, y, z);
	ret++; //zero stands for no vector

	return ret;
}

//empty constructor
ThreePointData::ThreePointData() :
		molPos(0), i1(0), i2(0), i3(0), l1(0), l2(0), l3(0), x(0), y(0), z(0), theta2(
				0), phi2(0), theta3(0), phi3(0), extra1(0), extra2(0), extra3(
				0), weight(0), nrot(0)
{

}

//initialize x,y,z - used for creating min/max object with external sorting
ThreePointData::ThreePointData(signed val) :
		molPos(0), i1(0), i2(0), i3(0), l1(val), l2(val), l3(val), x(0), y(0), z(
				0), theta2(0), phi2(0), theta3(0), phi3(0), extra1(0), extra2(
				0), extra3(0), weight(0), nrot(0)
{

}

//create the canonical three point from the points and indexes i,j,k
//return simple finger class
ThreePointData::ThreePointData(unsigned offset, double w, unsigned nr,
		const vector<PharmaPoint>& points, unsigned i, unsigned j, unsigned k)
{
	//points are always ordered by pharma kind first
	Triplet tri(points, i, j, k);
	array<PharmaIndex, 3> PIs = tri.getPoints();
	array<TripletRange, 3> ranges = tri.getRanges();

	fingerprint.set(PIs[0].index, PIs[1].index, PIs[2].index, points);

	weight = reduceWeight(w);
	nrot = reduceRotatable(nr);

	//points and lines are in canonical order, initialize struct
	molPos = offset;
	l1 = ranges[0].length;
	l2 = ranges[1].length;
	l3 = ranges[2].length;

	i1 = PIs[0].index;
	i2 = PIs[1].index;
	i3 = PIs[2].index;

	//printf("P(%d,%d,%d)  I(%d,%d,%d)  L(%d,%d,%d)\n",PIs[0].pharmaIndex,
	//		PIs[1].pharmaIndex,PIs[2].pharmaIndex, (int)i1, (int)i2, (int)i3, (int)l1,(int)l2,(int)l3);

	double X = PIs[0].point->x;
	double Y = PIs[0].point->y;
	double Z = PIs[0].point->z;

	x = reduceFloat(X);
	y = reduceFloat(Y);
	z = reduceFloat(Z);

	double d2 = PharmaPoint::pharmaDist(*PIs[1].point, *PIs[0].point);
	theta2 = reduceAngle(atan2(PIs[1].point->y - Y, PIs[1].point->x - X));
	phi2 = reduceAngle(acos((PIs[1].point->z - Z) / d2));

	double d3 = PharmaPoint::pharmaDist(*PIs[2].point, *PIs[0].point);
	theta3 = reduceAngle(atan2(PIs[2].point->y - Y, PIs[2].point->x - X));
	phi3 = reduceAngle(acos((PIs[2].point->z - Z) / d3));

	extra1 = PIs[0].point->size;
	extra2 = PIs[1].point->size;
	extra3 = PIs[2].point->size;

	//alternatively, vector info
	if (PIs[0].point->pharma->getVectors != NULL)
		extra1 = vecValue(PIs[0].point, PIs[0].point, PIs[1].point,
				PIs[2].point);
	if (PIs[1].point->pharma->getVectors != NULL)
		extra2 = vecValue(PIs[1].point, PIs[0].point, PIs[1].point,
				PIs[2].point);
	if (PIs[2].point->pharma->getVectors != NULL)
		extra3 = vecValue(PIs[2].point, PIs[0].point, PIs[1].point,
				PIs[2].point);

}

//return number of rotatable bonds according to smarts expression above
unsigned countRotatableBonds(OBMol& mol)
{
	//count number of rotatable bounds
	return mol.NumRotors();
}
