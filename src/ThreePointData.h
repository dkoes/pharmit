/*
Pharmer: Efficient and Exact 3D Pharmacophore Search
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
 * ThreePointData.h
 *
 *  Created on: Aug 6, 2010
 *      Author: dkoes
 */

#include "pharmarec.h"
#include "TripletFingerprint.h"
#include "SimpleFingers.h"

#ifndef THREEPOINTDATA_H_
#define THREEPOINTDATA_H_

typedef unsigned long uint128_t __attribute__((mode(TI)));


//at most 2^24 (16.7 million) mols per a database
#define TPD_MOLID_BITS 24
#define TPD_MOLDATA_BITS 40
#define TPD_LENGTH_BITS 16
#define TPD_COORD_BITS 16
#define TPD_ANGLE_BITS 16
#define TPD_INDEX_BITS 8
#define EXTRA_BITS 6
#define WEIGHT_BITS 10
#define ROTATABLE_BITS 4
// a struct representing a single 3 point pharmacophore
struct ThreePointData
{
	unsigned long molPos: TPD_MOLDATA_BITS; //location within moldata, each unique mol starts such that upper MOLID bits are the mol's id

	//index of each pharmapoints, for correspondence
	unsigned i1: TPD_INDEX_BITS;
	unsigned i2: TPD_INDEX_BITS;
	unsigned i3: TPD_INDEX_BITS;
//64 bits
	//lengths of the pharmacophore triangle sides
	unsigned l1: TPD_LENGTH_BITS;
	unsigned l2: TPD_LENGTH_BITS;
	unsigned l3: TPD_LENGTH_BITS;
//112 bits
	//actual spatial coordinates of mol, reduced
	//we store the cartesian coordinates of the first point, and then
	//spherical coordinates for the other two points relative to the first
	//we could get away with one less angle, but I am concerned about lose of precision
	//and computational overhead
	signed x: TPD_COORD_BITS;
	signed y: TPD_COORD_BITS;
	signed z: TPD_COORD_BITS;
//160 bits
	unsigned theta2: TPD_ANGLE_BITS;
	unsigned phi2: TPD_ANGLE_BITS;
	unsigned theta3: TPD_ANGLE_BITS;
	unsigned phi3: TPD_ANGLE_BITS;
//224 bits
	unsigned extra1: EXTRA_BITS;
	unsigned extra2: EXTRA_BITS;
	unsigned extra3: EXTRA_BITS;
//242 bits
	unsigned weight: WEIGHT_BITS;
//252 bits
	unsigned nrot: ROTATABLE_BITS; //number of rotatable bounds
//256 bits
	TripletFingerprint fingerprint;
//512 bits - 64 bytes per triplet
	ThreePointData(unsigned offset, double w, unsigned numrot, const vector<PharmaPoint>& points, unsigned i, unsigned j, unsigned k);
	ThreePointData();
  ThreePointData(signed val);

	static unsigned reduceAngle(double d);
	static double unreduceAngle(int val);
	static short reduceFloat(double d);
	static unsigned reduceFloatBigUnsigned(double d); //don't truncate to short
	static unsigned short reduceFloatUnsigned(double d);
	static double unreduceFloat(int val);
	static unsigned vecValue(const PharmaPoint *pt, const PharmaPoint *I, const PharmaPoint *J, const PharmaPoint *K);

	static unsigned reduceWeight(double d);
	static double unreduceWeight(unsigned v);

	static unsigned reduceRotatable(unsigned nr);

	signed x1() const { return x; };
	signed y1() const { return y; };
	signed z1() const { return z; };

	signed x2() const { return round((x) + (l1) * cos(unreduceAngle(theta2)) * sin(unreduceAngle(phi2))); }
	signed y2() const { return round((y) + (l1)* sin(unreduceAngle(theta2)) * sin(unreduceAngle(phi2))); }
	signed z2() const { return round((z) + (l1) * cos(unreduceAngle(phi2))); }

	signed x3() const { return round((x) + (l3) * cos(unreduceAngle(theta3)) * sin(unreduceAngle(phi3))); }
	signed y3() const { return round((y) + (l3)* sin(unreduceAngle(theta3)) * sin(unreduceAngle(phi3))); }
	signed z3() const { return round((z) + (l3) * cos(unreduceAngle(phi3))); }

	unsigned molID() const { return molPos >> (TPD_MOLDATA_BITS-TPD_MOLID_BITS); }
};



#endif /* THREEPOINTDATA_H_ */
