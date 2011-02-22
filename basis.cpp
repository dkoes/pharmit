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
 * basis.cpp
 *
 *  Created on: Dec 14, 2009
 *      Author: dkoes
 */

#include "basis.h"

const double minDist = .2; //minimum distance between basis points


ostream& operator<<(ostream& out, const CoordinateBasis& x)
{
	out << x.origin << "\n" << x.transform;
	return out;
}

//given the passed points, generate a new canonical basis that can be
//used to replot the other points
//use the order specified if fixedOrder is true
//return true if a valid basis is found
bool CoordinateBasis::constructBasis(const PharmaPoint&i, const PharmaPoint& j,
		const PharmaPoint& k, bool fixedOrder)
{
	unsigned iindex = 0;
	unsigned jindex = 1;
	unsigned kindex = 2;
	valid = false;
	//first check distances
	double d_ij = PharmaPoint::pharmaDist(i, j);
	if (d_ij < minDist)
		return false;
	double d_jk = PharmaPoint::pharmaDist(j, k);
	if (d_jk < minDist)
		return false;
	double d_ki = PharmaPoint::pharmaDist(k, i);
	if (d_ki < minDist)
		return false;

	//find the largest and second largest vectors with a joined tail
	vector3 first;
	vector3 second;

	if (fixedOrder)
	{
		//center on k
		origin.Set(k.x, k.y, k.z);
		first.Set(i.x - k.x, i.y - k.y, i.z - k.z);
		second.Set(j.x - k.x, j.y - k.y, j.z - k.z);
		basisOrder[0] = kindex;
		basisOrder[1] = iindex;
		basisOrder[2] = jindex;
		basisDistances[0] = d_ki;
		basisDistances[1] = d_jk;
		basisDistances[2] = d_ij;
	}
	else
	{
		if (d_ij > d_jk && d_ij > d_ki)
		{ //ij largest
			if (d_jk > d_ki)
			{
				//centered at j
				origin.Set(j.x, j.y, j.z);
				first.Set(i.x - j.x, i.y - j.y, i.z - j.z);
				second.Set(k.x - j.x, k.y - j.y, k.z - j.z);
				basisOrder[0] = jindex;
				basisOrder[1] = iindex;
				basisOrder[2] = kindex;
				basisDistances[0] = d_ij;
				basisDistances[1] = d_jk;
				basisDistances[2] = d_ki;
			}
			else
			{ //d_ki second largest, centered at i
				origin.Set(i.x, i.y, i.z);
				first.Set(j.x - i.x, j.y - i.y, j.z - i.z);
				second.Set(k.x - i.x, k.y - i.y, k.z - i.z);
				basisOrder[0] = iindex;
				basisOrder[1] = jindex;
				basisOrder[2] = kindex;
				basisDistances[0] = d_ij;
				basisDistances[1] = d_ki;
				basisDistances[2] = d_jk;
			}
		}
		else if (d_jk > d_ij && d_jk > d_ki)
		{ //jk largest
			if (d_ij > d_ki)
			{ //centered at j
				origin.Set(j.x, j.y, j.z);
				first.Set(k.x - j.x, k.y - j.y, k.z - j.z);
				second.Set(i.x - j.x, i.y - j.y, i.z - j.z);
				basisOrder[0] = jindex;
				basisOrder[1] = kindex;
				basisOrder[2] = iindex;
				basisDistances[0] = d_jk;
				basisDistances[1] = d_ij;
				basisDistances[2] = d_ki;
			}
			else
			{ //d_ki largest, centered at k
				origin.Set(k.x, k.y, k.z);
				first.Set(j.x - k.x, j.y - k.y, j.z - k.z);
				second.Set(i.x - k.x, i.y - k.y, i.z - k.z);
				basisOrder[0] = kindex;
				basisOrder[1] = jindex;
				basisOrder[2] = iindex;
				basisDistances[0] = d_jk;
				basisDistances[1] = d_ki;
				basisDistances[2] = d_ij;
			}
		}
		else
		{//ki largest (or equilateral, which we don't handle properly)
			if (d_ij > d_jk)
			{ //centered on i
				origin.Set(i.x, i.y, i.z);
				first.Set(k.x - i.x, k.y - i.y, k.z - i.z);
				second.Set(j.x - i.x, j.y - i.y, j.z - i.z);
				basisOrder[0] = iindex;
				basisOrder[1] = kindex;
				basisOrder[2] = jindex;
				basisDistances[0] = d_ki;
				basisDistances[1] = d_ij;
				basisDistances[2] = d_jk;
			}
			else
			{ //jk larger, centered on k
				origin.Set(k.x, k.y, k.z);
				first.Set(i.x - k.x, i.y - k.y, i.z - k.z);
				second.Set(j.x - k.x, j.y - k.y, j.z - k.z);
				basisOrder[0] = kindex;
				basisOrder[1] = iindex;
				basisOrder[2] = jindex;
				basisDistances[0] = d_ki;
				basisDistances[1] = d_jk;
				basisDistances[2] = d_ij;
			}
		}
	}

	//now take the cross product
	vector3 product = cross(first, second);
	if (product.length() < minDist)
		return false; //too close to colinear, a smaller cutoff coudl probably be used

	//first, second, and product define a new basis
	//orthonormalize it to eliminate skew (Gram-Schmidt)
	first.normalize();
	//u2 = v2 - proj_u1(v2)
	second = second - dot(second, first) * first;
	second.normalize();
	//u3 = v3 - proj_u1(v3) - proj_u2(v3)
	product = product - dot(product, first) * first - dot(product, second)
			* second;
	product.normalize();

	//now have an orthonormal basis
	//compute the transformation matrix
	transform = matrix3x3(first, second, product).transpose().inverse();

	valid = true;
	return true;
}


//generate new coordinates for x,y, and z using current basis
void CoordinateBasis::replot(double x, double y, double z, float &nx,
		float& ny, float& nz) const
{
	vector3 orig(x, y, z);
	orig -= origin;
	vector3 res = transform * orig;
	nx = res.x();
	ny = res.y();
	nz = res.z();
	cout.precision(8);
}

void CoordinateBasis::replot_vector(double x, double y, double z, float &nx,
		float& ny, float& nz) const
{
	vector3 orig(x, y, z);
	vector3 res = transform * orig;
	nx = res.x();
	ny = res.y();
	nz = res.z();
}

