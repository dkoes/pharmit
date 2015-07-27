/*
 * MGrid.cpp
 *
 *  Created on: Nov 21, 2011
 *      Author: dkoes
 */

#include "MGrid.h"
#include "molecules/Molecule.h"

void MGrid::gridToPoint(unsigned long g, double& x, double& y, double& z) const
{
	double res = resolution;
	unsigned length = dimension / res;
	z = g % length;
	g /= length;
	y = g % length;
	g /= length;
	x = g;

	x *= res;
	y *= res;
	z *= res;

	x -= dimension / 2;
	y -= dimension / 2;
	z -= dimension / 2;

	x += res / 2.0;
	y += res / 2.0;
	z += res / 2.0;
}

int MGrid::pointToGrid(double x, double y, double z) const
{
	double dim = dimension;
	double res = resolution;
	unsigned len = dim / res;
	x += dim / 2;
	y += dim / 2;
	z += dim / 2;

	if (x < 0)
		return -1;
	if (x >= dim)
		return -1;
	if (y < 0)
		return -1;
	if (y >= dim)
		return -1;
	if (z < 0)
		return -1;
	if (z >= dim)
		return -1;

	x /= res;
	y /= res;
	z /= res;

	unsigned ret = x;
	ret = ret * len + y;
	ret = ret * len + z;
	return ret;
}

//return the "radius" (half length) of a chord a distance d from the
//center of a circle of radius r
static double chordRadius(double r, double d)
{
	if (d > r)
		return 0;
	return sqrt(r * r - d * d);
}

//mark a line length 2*r in the grid that's a single chord of a sphere at x,y, z
void MGrid::markZChord(double x, double y, double z, double r)
{
	int start = pointToGrid(x, y, z - r);
	int end = pointToGrid(x, y, z + r);
	if (start < 0 || end < 0)
		return;
	grid.set_range(start, end, true);
}

//makr a circle in grid centered at x,y,z with radius r
void MGrid::markYZCircle(double x, double y, double z, double r)
{
	for (double d = 0; d <= r; d += resolution)
	{
		double cr = chordRadius(r, d);
		markZChord(x, y + d, z, cr);
		if (d != 0)
			markZChord(x, y - d, z, cr);
	}
}

void MGrid::markXYZSphere(double x, double y, double z, double r)
{
	if (sphereInGrid(x, y, z, r))
	{
		//mark all yz circles
		for (double d = 0; d <= r; d += resolution)
		{
			double cr = chordRadius(r, d);
			markYZCircle(x + d, y, z, cr);
			if (d != 0)
				markYZCircle(x - d, y, z, cr);
		}
	}
}

bool MGrid::sphereInGrid(float x, float y, float z, float r) const
{
	float mind = -dimension / 2 - r;
	float maxd = dimension / 2 + r;

	if (x < mind || y < mind || z < mind)
		return false;
	if (x > maxd || y > maxd || z > maxd)
		return false;
	return true;
}

//this downsamples from from
void MGrid::copyFrom(const MGrid& from)
{
	grid.clear();
	bvect::enumerator en = from.grid.first();
	bvect::enumerator en_end = from.grid.end();

	while (en < en_end)
	{
		unsigned g = *en;
		double x, y, z;
		from.gridToPoint(g, x, y, z);
		grid.set(pointToGrid(x, y, z));
		++en;
	}
}

void MGrid::setPoint(float x, float y, float z)
{
	unsigned p  = pointToGrid(x,y,z);
	grid.set(p, true);
}

//given the solvent accessible surface, make an approximation to the
//solvent excluded molecular surface by drawing spheres from the boundary
//of the sa and noting where they don't intersect
void MGrid::makeSurface(const MGrid& sagrid, const MGrid& lesssagrid,
		double probe)
{
	bvect shell = sagrid.grid - lesssagrid.grid;
	bvect::enumerator en = shell.first();
	bvect::enumerator en_end = shell.end();

	MGrid reachable(dimension, resolution);
	while (en < en_end)
	{
		unsigned g = *en;
		double x, y, z;
		sagrid.gridToPoint(g, x, y, z);
		if (!sagrid.isInteriorPoint(x, y, z))
		{
			reachable.markXYZSphere(x, y, z, probe);
		}
		++en;
	}

	grid |= sagrid.grid - reachable.grid;
}

//return true if all neighbors are set
bool MGrid::isInteriorPoint(float x, float y, float z) const
{
	//all neighbors must be set
	for (int xm = -1; xm <= 1; xm++)
	{
		float xc = x + resolution * xm;
		for (int ym = -1; ym <= 1; ym++)
		{
			float yc = y + resolution * ym;
			for (int zm = -1; zm <= 1; zm++)
			{
				float zc = z + resolution * zm;
				if (!test(xc, yc, zc))
					return false;
			}
		}
	}
	return true;
}

//a point that doesn't have neighbors on all its faces
//(ignoring edge/corner neighbors)
bool MGrid::isExposedPoint(float x, float y, float z) const
{
	if (test(x + resolution, y, z) && test(x - resolution, y, z)
			&& test(x, y + resolution, z) && test(x, y - resolution, z)
			&& test(x, y, z + resolution) && test(x, y, z - resolution))
		return false;
	return true;
}

//a point that is all by itself
bool MGrid::isSolitaryPoint(float x, float y, float z) const
{
	if (!test(x + resolution, y, z) && !test(x - resolution, y, z)
			&& !test(x, y + resolution, z) && !test(x, y - resolution, z)
			&& !test(x, y, z + resolution) && !test(x, y, z - resolution))
		return true;
	return false;
}


//remove all gridpoints that are exposed points
void MGrid::shrinkByOne()
{
	bvect::enumerator en = grid.first();
	bvect::enumerator en_end = grid.end();

	bvect toRemove;
	while (en < en_end)
	{
		unsigned g = *en;
		double x, y, z;
		gridToPoint(g, x, y, z);
		if (isExposedPoint(x, y, z))
		{
			toRemove.set(g);
		}
		++en;
	}

	grid -= toRemove;
}

void MGrid::growByOne()
{
	bvect::enumerator en = grid.first();
	bvect::enumerator en_end = grid.end();

	bvect toAdd;
	while (en < en_end)
	{
		unsigned g = *en;
		double x, y, z;
		gridToPoint(g, x, y, z);

		if (!test(x+resolution,y,z))
			toAdd.set(pointToGrid(x+resolution,y,z));
		if (!test(x-resolution,y,z))
			toAdd.set(pointToGrid(x-resolution,y,z));

		if (!test(x,y+resolution,z))
			toAdd.set(pointToGrid(x,y+resolution,z));
		if (!test(x,y-resolution,z))
			toAdd.set(pointToGrid(x,y-resolution,z));

		if (!test(x,y,z+resolution))
			toAdd.set(pointToGrid(x,y,z+resolution));
		if (!test(x,y,z-resolution))
			toAdd.set(pointToGrid(x,y,z-resolution));

		++en;
	}
	grid |= toAdd;
}


//reduce the size of the object by the specified amount
void MGrid::shrink(double amount)
{
	unsigned num = ceil(amount / resolution);
	for (unsigned i = 0; i < num; i++)
	{
		shrinkByOne();
		if(grid.count() == 0)
			break;
	}

}

//grow the size of the object by the specified amount
void MGrid::grow(double amount)
{
	unsigned num = ceil(amount / resolution);
	for (unsigned i = 0; i < num; i++)
	{
		growByOne();
	}
}

//return all set points in vector
void MGrid::getSetPoints(vector<MGrid::Point>& points) const
{
	points.clear();
	bvect::enumerator en = grid.first();
	bvect::enumerator en_end = grid.end();
	while (en < en_end)
	{
		unsigned g = *en;
		double x, y, z;
		gridToPoint(g, x, y, z);
		points.push_back(Point(x,y,z));
		++en;
	}

}
