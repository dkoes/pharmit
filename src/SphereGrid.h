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
 * SphereGrid.h
 *
 * Representation of gridded points on a unit sphere.
 *  Created on: May 6, 2010
 *      Author: dkoes
 */

#ifndef PHARMITSERVER_SPHEREGRID_H_
#define PHARMITSERVER_SPHEREGRID_H_

#include <vector>
#include <openbabel/math/vector3.h>

using namespace std;

#define SPHEREGRID_BITS (5)
class SphereGrid
{
	struct SphereGridPoint
	{
		unsigned index;
		vector<unsigned> neighbors;
		OpenBabel::vector3 pnt;

		SphereGridPoint(): index(0) {}
		SphereGridPoint(unsigned i, double X, double Y, double Z): index(i), pnt(X,Y,Z) {}

		double sqdistance(const SphereGridPoint& rhs)
		{
			return pnt.distSq(rhs.pnt);
		}
	};


	vector<SphereGridPoint> points; //indexed by index

	void addCoords(const vector<OpenBabel::vector3>& coords);
	void makeDodecahedron(); //20 points
	void makeIcosahedron(); //12 points, dual of dodecahedron
	void connectPoints();
public:

	//return a grid point (0..32)
	unsigned pointToGrid(double x, double y, double z);
	//return a mask of all grid points within angle of pt
	unsigned long searchMask(unsigned pt, double angle);

	SphereGrid()
	{
		//by default combine two platonic solids and get 32 points
		//not quite evenly spaced, but good enough
		makeDodecahedron();
		makeIcosahedron();
		connectPoints();
	}
	virtual ~SphereGrid() {}
};

extern SphereGrid sphereGrid;

#endif /* PHARMITSERVER_SPHEREGRID_H_ */
