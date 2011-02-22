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
 * RMSD.cpp
 *
 *  Created on: Feb 3, 2010
 *      Author: dkoes
 */

#include "RMSD.h"

ostream& operator<<(ostream& out, const RMSDResult& r)
{
	out << r.rotation << "\n";
	out << r.translation << "\n";
	return out;
}


//calculate  rmsd of passed n points
//use dual number quaternions (Walker, Shao, and Volz, Estimating 3-D Location Parameters Using Dual Number Quaternions, 1991
RMSDResult calculateRMSD(const double *ref, const double *fit, unsigned n)
{
	Mat4x4 C1 = Mat4x4::Zero();
	Mat4x4 C3 = Mat4x4::Zero();

	//calculate C matrices
	//C2 first; this is the full inverse(c_2+c_2^t)
	double sum = n;
	double c2mult = 1.0/(2.0*sum);

	double constant = 0;
	//C1 and C3, require W and Q
	for(unsigned index = 0; index < n; index++)
	{
		double beta = 1.0;
		double rx = ref[3*index]/2.0;
		double ry = ref[3*index+1]/2.0;
		double rz = ref[3*index+2]/2.0;

		double fx = fit[3*index]/2.0;
		double fy = fit[3*index+1]/2.0;
		double fz = fit[3*index+2]/2.0;

		FloatType qdata[] __attribute__ ((aligned (16))) =
		{0, rz, -ry, -rx,
		-rz, 0, rx, -ry,
		ry, -rx, 0, -rz,
		rx, ry, rz, 0};
		Map<Mat4x4> q(qdata);

		FloatType wdata[]  __attribute__ ((aligned (16))) =
		{0, -fz, fy, -fx,
		fz, 0, -fx, -fy,
		-fy, fx, 0, -fz,
		fx, fy, fz, 0};
		Map<Mat4x4> w(wdata);

		constant += beta*(rx*rx+ry*ry+rz*rz + fx*fx+fy*fy+fz*fz);

		C1 += beta*q.transpose() * w;
		C3 += beta*(w - q);
	}

	C3 *= 2;
	//now compute A
	Mat4x4 A = 0.5*(C3.transpose()*c2mult*C3) + C1 + C1.transpose();
	SelfAdjointEigenSolver<Mat4x4> esolver(A);

	const Vec4& evals = esolver.eigenvalues();
	const Mat4x4& evecs = esolver.eigenvectors();

	//find max, is it always the last one?
	unsigned maxIndex = 0;
	double maxValue = -HUGE_VAL;
	for(unsigned i = 0; i < 4; i++)
	{
		if(evals[i] > maxValue)
		{
			maxIndex = i;
			maxValue = evals[i];
		}
	}

	Vec4 r = evecs.col(maxIndex);
	Vec4 s = -1.0*c2mult*C3*r;

	FloatType wrtdata[] __attribute__ ((aligned (16))) = {
	r.coeff(3), r.coeff(2), -r.coeff(1), r.coeff(0),
	-r.coeff(2), r.coeff(3), r.coeff(0), r.coeff(1),
	r.coeff(1), -r.coeff(0), r.coeff(3), r.coeff(2),
	-r.coeff(0), -r.coeff(1), -r.coeff(2), r.coeff(3)};

	Map<Mat4x4> wt(wrtdata);

	FloatType qrdata[] __attribute__ ((aligned (16))) = {
	r.coeff(3), r.coeff(2),-r.coeff(1), -r.coeff(0),
	-r.coeff(2), r.coeff(3), r.coeff(0), -r.coeff(1),
	r.coeff(1), -r.coeff(0), r.coeff(3), -r.coeff(2),
	r.coeff(0), r.coeff(1), r.coeff(2), r.coeff(3)};

	Map<Mat4x4> q(qrdata);

	Vec4 t = 2*wt*s;
	Mat4x4 R = wt*q;
	double val = 0;

	if(maxValue < constant) //for zero rmsd can be off by a little bit and end up with sqrt of negative
		val = sqrt((4*(constant - maxValue))/n);


	RMSDResult res(val, t.block<3,1>(0,0), R.block<3,3>(0,0));
	return res;
}



//calculate weighted rmsd of passed n points
//use dual number quaternions (Walker, Shao, and Volz, Estimating 3-D Location Parameters Using Dual Number Quaternions, 1991
RMSDResult calculateRMSD(const double *ref, const double *fit, const double *weights, unsigned n)
{
	Mat4x4 C1 = Mat4x4::Zero();
	Mat4x4 C3 = Mat4x4::Zero();

	//calculate C matrices
	//C2 first; this is the full inverse(c_2+c_2^t)
	double sum = 0;
	for(unsigned i = 0; i < n; i++)
		sum+= weights[i];
	double c2mult = 1.0/(2.0*sum);

	double constant = 0;
	//C1 and C3, require W and Q
	for(unsigned index = 0; index < n; index++)
	{
		double beta = weights[index];
		double rx = ref[3*index]/2.0;
		double ry = ref[3*index+1]/2.0;
		double rz = ref[3*index+2]/2.0;

		double fx = fit[3*index]/2.0;
		double fy = fit[3*index+1]/2.0;
		double fz = fit[3*index+2]/2.0;

		FloatType qdata[] __attribute__ ((aligned (16))) =
		{0, rz, -ry, -rx,
		-rz, 0, rx, -ry,
		ry, -rx, 0, -rz,
		rx, ry, rz, 0};
		Map<Mat4x4> q(qdata);

		FloatType wdata[]  __attribute__ ((aligned (16))) =
		{0, -fz, fy, -fx,
		fz, 0, -fx, -fy,
		-fy, fx, 0, -fz,
		fx, fy, fz, 0};
		Map<Mat4x4> w(wdata);

		constant += beta*(rx*rx+ry*ry+rz*rz + fx*fx+fy*fy+fz*fz);

		C1 += beta*q.transpose() * w;
		C3 += beta*(w - q);
	}

	C3 *= 2;
	//now compute A
	Mat4x4 A = 0.5*(C3.transpose()*c2mult*C3) + C1 + C1.transpose();
	SelfAdjointEigenSolver<Mat4x4> esolver(A);

	const Vec4& evals = esolver.eigenvalues();
	const Mat4x4& evecs = esolver.eigenvectors();

	//find max, is it always the last one?
	unsigned maxIndex = 0;
	double maxValue = -HUGE_VAL;
	for(unsigned i = 0; i < 4; i++)
	{
		if(evals[i] > maxValue)
		{
			maxIndex = i;
			maxValue = evals[i];
		}
	}

	Vec4 r = evecs.col(maxIndex);
	Vec4 s = -1.0*c2mult*C3*r;

	FloatType wrtdata[] __attribute__ ((aligned (16))) = {
	r.coeff(3), r.coeff(2), -r.coeff(1), r.coeff(0),
	-r.coeff(2), r.coeff(3), r.coeff(0), r.coeff(1),
	r.coeff(1), -r.coeff(0), r.coeff(3), r.coeff(2),
	-r.coeff(0), -r.coeff(1), -r.coeff(2), r.coeff(3)};

	Map<Mat4x4> wt(wrtdata);

	FloatType qrdata[] __attribute__ ((aligned (16))) = {
	r.coeff(3), r.coeff(2),-r.coeff(1), -r.coeff(0),
	-r.coeff(2), r.coeff(3), r.coeff(0), -r.coeff(1),
	r.coeff(1), -r.coeff(0), r.coeff(3), -r.coeff(2),
	r.coeff(0), r.coeff(1), r.coeff(2), r.coeff(3)};

	Map<Mat4x4> q(qrdata);

	Vec4 t = 2*wt*s;
	Mat4x4 R = wt*q;
	double val = 0;

	if(maxValue < constant) //for zero rmsd can be off by a little bit and end up with sqrt of negative
		val = sqrt((4*(constant - maxValue))/n);

	RMSDResult res(val, t.block<3,1>(0,0), R.block<3,3>(0,0));
	return res;
}
