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
 * SphereGrid.cpp
 *
 *  Created on: May 6, 2010
 *      Author: dkoes
 */


#include "SphereGrid.h"
#include <queue>

#define ADDPT(x,y,z) pts.push_back(vector3(x,y,z))

//twenty uniformly spaced points
void SphereGrid::makeDodecahedron()
{
	vector<vector3> pts;
	pts.reserve(20);

	ADDPT(-1,-1,-1);
	ADDPT(-1,-1,1);
	ADDPT(-1, 1,-1);
	ADDPT(-1, 1, 1);

	ADDPT(1,-1,-1);
	ADDPT(1,-1,1);
	ADDPT(1, 1,-1);
	ADDPT(1, 1, 1);

	double golden = (1+sqrt(5))/2;

	ADDPT(0,-1.0/golden, -golden);
	ADDPT(0,-1.0/golden, golden);
	ADDPT(0,1.0/golden, -golden);
	ADDPT(0,1.0/golden, golden);

	ADDPT(-1.0/golden, -golden,0);
	ADDPT(-1.0/golden, golden,0);
	ADDPT(1.0/golden, -golden,0);
	ADDPT(1.0/golden, golden,0);

	ADDPT(-golden, 0, -1.0/golden);
	ADDPT(-golden, 0, 1.0/golden);
	ADDPT(golden, 0, -1.0/golden);
	ADDPT(golden, 0, 1.0/golden);

	addCoords(pts);
}

//12 points, dual of dodecahedron
void SphereGrid::makeIcosahedron()
{
	vector<vector3> pts;
	pts.reserve(12);

	double a = sqrt((5+sqrt(5))/10);
	double b = sqrt((5-sqrt(5))/10);

	ADDPT(a,b,0);
	ADDPT(a,-b,0);
	ADDPT(-a,b,0);
	ADDPT(-a,-b,0);

	ADDPT(0,a,b);
	ADDPT(0,a,-b);
	ADDPT(0,-a,b);
	ADDPT(0,-a,-b);

	ADDPT(b,0,a);
	ADDPT(b,0,-a);
	ADDPT(-b,0,a);
	ADDPT(-b,0,-a);

	addCoords(pts);
}


//connect nearest neighbors to points
void SphereGrid::connectPoints()
{
	//just use brute force approach
	for(unsigned i = 0, n = points.size(); i < n; i++)
	{
		double mindist = HUGE_VAL;
		//find min dist
		for(unsigned j = 0; j < n; j++)
		{
			if(j != i)
			{
				double d = points[i].sqdistance(points[j]);
				if(d < mindist)
					mindist = d;
			}
		}
		//find all neighbors with mindist
		for(unsigned j = 0; j < n; j++)
		{
			if(j != i)
			{
				double d = points[i].sqdistance(points[j]);
				if(d < mindist+.001)
				{
					points[i].neighbors.push_back(j);
				}
			}
		}
	}
}

//normalize and create sphere grid points
void SphereGrid::addCoords(const vector<vector3>& coords)
{
	for(unsigned i = 0, n = coords.size(); i < n; i++)
	{
		vector3 c = coords[i];
		c.normalize();
		points.push_back(SphereGridPoint(points.size(), c.x(), c.y(), c.z()));
	}
}


//return a grid point (0..32)
unsigned SphereGrid::pointToGrid(double x, double y, double z)
{
	//inefficient, lazy implementation; if performance ends up mattering,
	//do analytically
	vector3 pt(x,y,z);
	pt.normalize();

	unsigned ret = 0;
	double mindist = HUGE_VAL;
	for(unsigned i = 0, n = points.size(); i < n; i++)
	{
		double d = pt.distSq(points[i].pnt);
		if(d < mindist)
		{
			mindist = d;
			ret = i;
		}
	}
	return ret;
}

//return a mask of all grid points within angle (degrees) of pt
//angle should be <= pi
unsigned long SphereGrid::searchMask(unsigned pt, double angle)
{
	unsigned long ret = 0;
	angle = DEG_TO_RAD*angle;
	//compute max distance for angle
	double maxdist = 2*sin(angle/2);
	maxdist *= maxdist;
	//do a bfs of points from pt
	vector<bool> seen(points.size(), false);
	queue<unsigned> Q;
	Q.push(pt);
	seen[pt] = true;
	while(!Q.empty())
	{
		unsigned p = Q.front();
		Q.pop();
		ret |= (1UL << p);

		for(unsigned i = 0, n = points[p].neighbors.size(); i < n; i++)
		{
			unsigned next = points[p].neighbors[i];
			if(seen[next]) continue;
			seen[next] = true;
			//is this close enough?
			double dist = points[next].sqdistance(points[pt]);
			if(dist <= maxdist)
			{
				ret |= (1UL<<next);
				Q.push(next);
			}
		}
	}

	return ret;
}

SphereGrid sphereGrid; //global

