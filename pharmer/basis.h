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
 * basis.h
 *
 *  Created on: Dec 14, 2009
 *      Author: dkoes
 */

#ifndef BASIS_H_
#define BASIS_H_

#include "pharmarec.h"

#include <openbabel/math/vector3.h>
#include <openbabel/math/matrix3x3.h>


//given a trip of points, constructs a new basis and then can
//reposition passed coordinates
class CoordinateBasis
{
	OpenBabel::vector3 origin;
	OpenBabel::matrix3x3 transform;
	int basisOrder[3]; //canonicalized order of basis points
	float basisDistances[3];
	bool valid;
public:
	CoordinateBasis():valid(false)
	{

	}

	CoordinateBasis(const PharmaPoint&i, const PharmaPoint& j,
			const PharmaPoint& k, bool fixedOrder) :
		valid(false)
	{
		memset(basisOrder, -1, sizeof(basisOrder));
		memset(basisDistances, 0, sizeof(basisDistances));
		constructBasis(i, j, k, fixedOrder);
	}

	void makeIdentity()
	{
		valid = true;
		origin = OpenBabel::vector3(0,0,0);
		transform = OpenBabel::matrix3x3(1.0);
	}

	bool constructBasis(const PharmaPoint&i, const PharmaPoint& j,
			const PharmaPoint& k, bool fixedOrder);

	bool hasValidBasis()
	{
		return valid;
	}

	OpenBabel::vector3 getTranslate() const { return origin; }
	void setTranslate(OpenBabel::vector3 o) { origin = o; }

	int basisPoint(int i)
	{
		return basisOrder[i];
	} //return ith basis point

	void replot(double x, double y, double z, float &nx, float& ny, float& nz) const;
	void replot(double x, double y, double z, double &nx, double& ny, double& nz) const
	{
		float fx, fy, fz;
		replot(x,y,z,fx,fy,fz);
		nx = fx;
		ny = fy;
		nz = fz;
	}

	//keep rooted at origin
	void replot_vector(double x, double y, double z, float &nx, float& ny, float& nz) const;

	friend ostream& operator<<(ostream&, const CoordinateBasis& x);
};


#endif /* BASIS_H_ */
