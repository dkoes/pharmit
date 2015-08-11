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
 * RMSD.h
 *
 *  Created on: Feb 3, 2010
 *      Author: dkoes
 */

#ifndef PHARMITSERVER_RMSD_H_
#define PHARMITSERVER_RMSD_H_

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include <iostream>
#include <openbabel/mol.h>

using namespace std;

typedef double FloatType;
typedef Eigen::Matrix4d Mat4x4;
typedef Eigen::Vector4d Vec4;
typedef Eigen::Vector3d Vec3;
typedef Eigen::Matrix3d Mat3x3;

/*
 * This class stores the transformation of a molecules with its query rmsd.
 * The notion of rmsd is overloaded for shape searches to be the similarity
 */
class RMSDResult
{
	//floats to save space, actual computation must be done in double precision though
	float val;

	Eigen::Matrix3f rotation;
	Eigen::Vector3f translation;

	friend class RMSDCalculator;
public:
	RMSDResult(): val(0), rotation(Eigen::Matrix3f::Identity()),translation(Eigen::Vector3f::Zero())
	{

	}

	RMSDResult(double v, const Vec3& t, const Mat3x3& r): val(v), rotation(r.cast<float>()), translation(t.cast<float>()) {}

	void clear()
	{
		using namespace Eigen;
		val = 0;
		rotation = Matrix3f::Identity();
		translation = Vector3f::Zero();
	}

	double value() const { return val; }
	void setValue(double v) { val = v; }

	//modify points by rot/trans
	void reorient(vector<Eigen::Vector3f>& pnts) const
	{
		for(unsigned i = 0, n = pnts.size(); i < n; i++)
		{
			pnts[i] = rotation*pnts[i]+translation;
		}
	}

	void reorient(unsigned n, double *coords) const
	{
		using namespace Eigen;
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
		using namespace Eigen;
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

	void reorient(OpenBabel::OBMol &mol) const
	{
		unsigned n = mol.NumAtoms();
		double *coords = mol.GetCoordinates();
		reorient(n, coords);
	}

	const Eigen::Matrix3f& rotationMatrix() const { return rotation; }
	const Eigen::Vector3f& translationVector() const { return translation; }

	friend ostream& operator<<(ostream& out, const RMSDResult& r);
};

RMSDResult calculateRMSD(const double *ref, const double *fit, unsigned n);
RMSDResult calculateRMSD(const double *ref, const double *fit, const double *weights, unsigned n);




#endif /* PHARMITSERVER_RMSD_H_ */
