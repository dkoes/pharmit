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
 * BoundingBox.h
 *
 *  Created on: Aug 2, 2010
 *      Author: dkoes
 */

#ifndef BOUNDINGBOX_H_
#define BOUNDINGBOX_H_

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <cstdlib>
#include <climits>
#include <iostream>

using namespace std;
using namespace boost;

enum SplitType
{
	NoSplit,
	SplitXFixed, //split by a fixed amount along x axis
	SplitYFixed, //y axis
	SplitZFixed,
//z axis
};

typedef tuple<int, int, int> TripleInt;

//maintain a bounding box
struct BoundingBox
{
	short minx;
	short maxx;
	short miny;
	short maxy;
	short minz;
	short maxz;

	BoundingBox() :
		minx(SHRT_MAX), maxx(SHRT_MIN), miny(SHRT_MAX), maxy(SHRT_MIN), minz(
				SHRT_MAX), maxz(SHRT_MIN)
	{
	}

	void update(const TripleInt& ijk)
	{
		update(ijk.get<0>(), ijk.get<1>(), ijk.get<2>());
	}

	void update(short x, short y, short z)
	{
		if (x < minx)
			minx = x;
		if (x > maxx)
			maxx = x;
		if (y < miny)
			miny = y;
		if (y > maxy)
			maxy = y;
		if (z < minz)
			minz = z;
		if (z > maxz)
			maxz = z;
	}

	double volume() const
	{
		return (double)(maxx-minx)*(maxy-miny)*(maxz-minz);
	}

	//make largest possible
	void maximize()
	{
		minx = miny = minz = SHRT_MIN;
		maxx = maxy = maxz = SHRT_MAX;
	}

	void setX(short min, short max)
	{
		minx = min;
		maxx = max;
	}

	void setY(short min, short max)
	{
		miny = min;
		maxy = max;
	}

	void setZ(short min, short max)
	{
		minz = min;
		maxz = max;
	}

	//return true if the passed ranges overlap
	static bool hasOverlap(double min1, double max1, double min2, double max2)
	{
		if (min1 <= max2 && min1 >= min2)
			return true;
		if (max1 <= max2 && max1 >= min2)
			return true;
		if (min2 <= max1 && min2 >= min1)
			return true; //fully contained
		return false;
	}

	//return true if this overlaps with x
	bool hasOverlap(const BoundingBox& x) const
	{
		if (hasOverlap(minx, maxx, x.minx, x.maxx)
				&& hasOverlap(miny, maxy, x.miny, x.maxy)
				&& hasOverlap(minz, maxz, x.minz, x.maxz))
			return true;
		return false;
	}

	//return true if this is fully contained within x
	bool containedIn(const BoundingBox& x) const
	{
		if(minx < x.minx) return false;
		if(maxx > x.maxx) return false;
		if(miny < x.miny) return false;
		if(maxy > x.maxy) return false;
		if(minz < x.minz) return false;
		if(maxz > x.maxz) return false;

		return true;
	}

	//return true if passed split bisects box
	bool planeInBox(SplitType split, short splitVal) const
	{
		switch (split)
		{
		case SplitXFixed:
			return splitVal > minx && splitVal < maxx;
		case SplitYFixed:
			return splitVal > miny && splitVal < maxy;
		case SplitZFixed:
			return splitVal > minz && splitVal < maxz;
		default:
			abort();
		}
	}

	//narrow along a single dimension
	void narrow(SplitType split, short min, short max)
	{
		switch(split)
		{
		case SplitXFixed:
			minx = min;
			maxx = max;
			break;
		case SplitYFixed:
			miny = min;
			maxy = max;
			break;
		case SplitZFixed:
			minz = min;
			maxz = max;
			break;
		default:
			break;
		}
	}

	friend ostream& operator<<(ostream& out, const BoundingBox& box);
};

#endif /* BOUNDINGBOX_H_ */
