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
