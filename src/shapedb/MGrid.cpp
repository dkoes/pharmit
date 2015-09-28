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
void MGrid::makeSurface(const MGrid& sagrid, double probe)
{
	MGrid lesssagrid = sagrid;
	lesssagrid.shrinkByOne();

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

//make a face a distance r from pt along which axis
void MGrid::makeFace(Eigen::Vector3f pt, int which, double r,vector<Eigen::Vector3f>& vertices, vector<Eigen::Vector3f>& normals, vector<int>& faces)
{
	using namespace Eigen;
	Vector3f neigh = pt;
	neigh[which] += 2*r;
	if(!test(neigh.x(), neigh.y(), neigh.z()))
	{
		//exposed, make verts
		pt[which] += r;

		//one dimension is only in r, the other are +/- both ways to make corners
		Vector3f v1 = pt, v2 = pt, v3 = pt, v4 = pt;
		v1[(which+1)%3] += r;
		v1[(which+2)%3] += r;

		v2[(which+1)%3] += r;
		v2[(which+2)%3] -= r;

		v3[(which+1)%3] -= r;
		v3[(which+2)%3] -= r;

		v4[(which+1)%3] -= r;
		v4[(which+2)%3] += r;

		int i1,i2,i3,i4;
		i1 = vertices.size();
		vertices.push_back(v1);

		i2 = vertices.size();
		vertices.push_back(v2);

		i3 = vertices.size();
		vertices.push_back(v3);

		i4 = vertices.size();
		vertices.push_back(v4);

		//the normals are all the same - r
		Vector3f norm(0,0,0);
		norm[which] = -1.0; //becaue of counterclockwise winding??
		normals.push_back(norm);
		normals.push_back(norm);
		normals.push_back(norm);
		normals.push_back(norm);

		//two faces [v1,v2,v3] and [v1,v3,v4]
		faces.push_back(i1);
		faces.push_back(i2);
		faces.push_back(i3);

		faces.push_back(i1);
		faces.push_back(i3);
		faces.push_back(i4);
	}
}

//create a voxel mesh from a grid
//does not clear vectors, instead appends info
void MGrid::makeMesh(vector<Eigen::Vector3f>& vertices, vector<Eigen::Vector3f>& normals, vector<int>& faces)
{
	using namespace Eigen;
	bvect::enumerator en = grid.first();
	bvect::enumerator en_end = grid.end();
	double x = 0, y = 0, z = 0;

	while (en < en_end)
	{
		unsigned g = *en;
		++en;
		gridToPoint(g, x, y , z);

		if(isExposedPoint(x,y,z))
		{
			double r = resolution/2.0;
			Vector3f pt(x,y,z);
			for(unsigned i = 0; i < 3; i++)
			{
				//doesn't make a face unless actually exposed
				makeFace(pt, i, r, vertices, normals, faces);
				makeFace(pt, i, -r, vertices, normals, faces);
			}
		}
	}
}

