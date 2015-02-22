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
 * FloatCoord
 *
 *  Created on: Aug 2, 2010
 *      Author: dkoes
 *
 *      A simple coordinate class
 */

#ifndef FLOATCOORD_H_
#define FLOATCOORD_H_

struct FloatCoord
{
	float x;
	float y;
	float z;

	FloatCoord(float a, float b, float c) :
		x(a), y(b), z(c)
	{
	}
	FloatCoord() :
		x(0), y(0), z(0)
	{
	}

	FloatCoord(float *xyz): x(xyz[0]), y(xyz[1]), z(xyz[2])
	{

	}

	FloatCoord& operator+=(const FloatCoord& rhs)
	{
		x += rhs.x;
		y += rhs.y;
		z += rhs.z;
		return *this;
	}

	FloatCoord& operator-=(const FloatCoord& rhs)
	{
		x -= rhs.x;
		y -= rhs.y;
		z -= rhs.z;
		return *this;
	}

	FloatCoord& operator/=(float r)
	{
		x /= r;
		y /= r;
		z /= r;
		return *this;
	}

	bool operator==(const FloatCoord& rhs) const
		{
		return x == rhs.x && y == rhs.y && z == rhs.z;
		}

	//square of distance to rhs
	double sqrDist(const FloatCoord& rhs) const
	{
		float d1 = x-rhs.x;
		float d2 = y-rhs.y;
		float d3 = z-rhs.z;
		return d1*d1+d2*d2+d3*d3;
	}
};


#endif /* FLOATCOORD_H_ */
