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
 * RMSD.h
 *
 *  Created on: Feb 3, 2010
 *      Author: dkoes
 */

#ifndef RMSD_H_
#define RMSD_H_

#include <openbabel/mol.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include <iostream>
using namespace Eigen;
using namespace OpenBabel;
using namespace std;

typedef double FloatType;
typedef Matrix4d Mat4x4;
typedef Vector4d Vec4;
typedef Vector3d Vec3;
typedef Matrix3d Mat3x3;

class RMSDResult
{
	//floats to save space, actual computation must be done in double precision though
	float val;

	Matrix3f rotation;
	Vector3f translation;

	friend class RMSDCalculator;
public:
	RMSDResult(): val(0), rotation(Matrix3f::Identity()),translation(Vector3f::Zero())
	{

	}

	RMSDResult(double v, const Vec3& t, const Mat3x3& r): val(v), rotation(r.cast<float>()), translation(t.cast<float>()) {}

	void clear()
	{
		val = 0;
		rotation = Matrix3f::Identity();
		translation = Vector3f::Zero();
	}

	double value() const { return val; }

	//modify points by rot/trans
	void reorient(vector<Vector3f>& pnts) const
	{
		for(unsigned i = 0, n = pnts.size(); i < n; i++)
		{
			pnts[i] = rotation*pnts[i]+translation;
		}
	}

	void reorient(unsigned n, double *coords) const
	{
		vector<Vector3f> pnts(n);
		for(unsigned i = 0; i < n; i++)
		{
			pnts[i] = Vector3f(coords[3*i],coords[3*i+1],coords[3*i+2]);
		}

		reorient(pnts);

		for(unsigned i = 0; i < n; i++)
		{
			coords[3*i] = pnts[i].coeff(0);
			coords[3*i+1] = pnts[i].coeff(1);
			coords[3*i+2] = pnts[i].coeff(2);
		}
	}

	void reorient(unsigned n, float *coords) const
	{
		vector<Vector3f> pnts(n);
		for(unsigned i = 0; i < n; i++)
		{
			pnts[i] = Vector3f(coords[3*i],coords[3*i+1],coords[3*i+2]);
		}

		reorient(pnts);

		for(unsigned i = 0; i < n; i++)
		{
			coords[3*i] = pnts[i].coeff(0);
			coords[3*i+1] = pnts[i].coeff(1);
			coords[3*i+2] = pnts[i].coeff(2);
		}
	}

	void reorient(OBMol &mol) const
	{
		unsigned n = mol.NumAtoms();
		double *coords = mol.GetCoordinates();
		reorient(n, coords);
	}

	const Matrix3f& rotationMatrix() const { return rotation; }
	const Vector3f& translationVector() const { return translation; }

	friend ostream& operator<<(ostream& out, const RMSDResult& r);
};

RMSDResult calculateRMSD(const double *ref, const double *fit, unsigned n);
RMSDResult calculateRMSD(const double *ref, const double *fit, const double *weights, unsigned n);




#endif /* RMSD_H_ */
