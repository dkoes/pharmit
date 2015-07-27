/*
 * Cube.h
 *
 * A simple class for representing a cube.
 *  Created on: Sep 28, 2011
 *      Author: dkoes
 */

#ifndef CUBE_H_
#define CUBE_H_

#include <cmath>
#include <cfloat>
#include <algorithm>
using namespace std;

/* A cube. Functions for consistently dividing into octants */
class Cube
{
public:
	float x, y, z; //bottom corner
	float dim;

	inline float squared(float v) const
	{
		return v * v;
	}

	//find the minimum distance between line segements a1-a2 and b1-b2
	inline float min1dist(float a1, float a2, float b1, float b2) const
	{
		//eh, just try all combos
		return min(min(fabs(a1-b2),fabs(a1-b1)),min(fabs(a2-b2),fabs(a2-b1)));
	}
public:

	Cube(float d) :
		x(0), y(0), z(0), dim(d)
	{

	}

	Cube(float d, float X, float Y, float Z): x(X), y(Y), z(Z), dim(d) {}

	float getDimension() const
	{
		return dim;
	}

	float volume() const
	{
		return dim*dim*dim;
	}

	void getCenter(float& cx, float& cy, float& cz) const
	{
		cx = x+dim/2;
		cy = y+dim/2;
		cz = z+dim/2;
	}

	void getBottomCorner(float& bx, float& by, float& bz) const
	{
		bx = x;
		by = y;
		bz = z;
	}

	//return i'th octant
	Cube getOctant(unsigned i) const
	{
		Cube res = *this;
		res.dim /= 2.0;

		switch (i)
		{
		case 0:
			break;
		case 1:
			res.x += res.dim;
			break;
		case 2:
			res.y += res.dim;
			break;
		case 3:
			res.x += res.dim;
			res.y += res.dim;
			break;
		case 4:
			res.z += res.dim;
			break;
		case 5:
			res.x += res.dim;
			res.z += res.dim;
			break;
		case 6:
			res.y += res.dim;
			res.z += res.dim;
			break;
		case 7:
			res.x += res.dim;
			res.y += res.dim;
			res.z += res.dim;
			break;
		default:
			abort();
			break;
		}
		return res;
	}

	//shortest distance from this cube to rhs
	float minDist(const Cube& rhs) const
	{
		float minx = min1dist(x, x+dim, rhs.x, rhs.x+rhs.dim);
		float miny = min1dist(y, y+dim, rhs.y, rhs.y+rhs.dim);
		float minz = min1dist(z, z+dim, rhs.z, rhs.z+rhs.dim);

		return sqrt(minx*minx+miny*miny+minz*minz);
	}
};

#endif /* CUBE_H_ */
