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
 * basis.h
 *
 *  Created on: Dec 14, 2009
 *      Author: dkoes
 */

#ifndef PHARMITSERVER_BASIS_H_
#define PHARMITSERVER_BASIS_H_

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


#endif /* PHARMITSERVER_BASIS_H_ */
