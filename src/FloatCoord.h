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
 * FloatCoord
 *
 *  Created on: Aug 2, 2010
 *      Author: dkoes
 *
 *      A simple coordinate class
 */

#ifndef PHARMITSERVER_FLOATCOORD_H_
#define PHARMITSERVER_FLOATCOORD_H_

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


#endif /* PHARMITSERVER_FLOATCOORD_H_ */
