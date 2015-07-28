/*
 * ShapeObj.cpp
 *
 *  Created on: Jul 28, 2015
 *      Author: dkoes
 */

#include <ShapeObj.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>


//calculate the moment of intertia for axis C on both sides
pair<double,double> ShapeObj::calcSplitMoments(OBMol& mol, unsigned C)
{
	double above = 0;
	double below = 0;
	unsigned A = (C+1)%3;
	unsigned B = (C+2)%3;

	OBAtom *atom;
	vector<OBAtom*>::iterator j;
	//find center of mass
	for (atom = mol.BeginAtom(j); atom; atom = mol.NextAtom(j))
	{
		vector3 pt = atom->GetVector();
		double val = pt[A]*pt[A]+pt[B]*pt[B];
		if(pt[C] < 0)
		{
			below += val;
		}
		else
		{
			above += val;
		}
	}

	return pair<double,double>(below,above);
}


void ShapeObj::normalizeMol(OBMol& mol)
{
	using namespace Eigen;
	//we need to normalize the molecular pose
	double center[3] = {0,0,0};
	Matrix3d I = Matrix3d::Zero();

	OBAtom *atom;
	double cnt = 0;
	vector<OBAtom*>::iterator j;
	//find center of mass
	for (atom = mol.BeginAtom(j); atom; atom = mol.NextAtom(j))
	{
		center[0] += atom->x();
		center[1] += atom->y();
		center[2] += atom->z();
		cnt += 1;
	}

	center[0] /= cnt;
	center[1] /= cnt;
	center[2] /= cnt;

	mol.Translate(-vector3(center));

	double Ixx = 0, Iyy = 0, Izz = 0, Ixy = 0, Ixz = 0, Iyz = 0;
	//calculate inertial tensor
	for (atom = mol.BeginAtom(j); atom; atom = mol.NextAtom(j))
	{
		double x = atom->x() - center[0];
		double y = atom->y() - center[1];
		double z = atom->z() - center[2];

		Ixx += (y * y + z * z);
		Iyy += (x * x + z * z);
		Izz += (x * x + y * y);
		Ixy += -x * y;
		Ixz += -x * z;
		Iyz += -y * z;
	}

	I << Ixx, Ixy, Ixz,
			Ixy, Iyy, Iyz,
			Ixz, Iyz, Izz;

	SelfAdjointEigenSolver<Matrix3d> esolver(I);
	const Matrix3d& evecs = esolver.eigenvectors();

	Matrix3d principalAxes = evecs; //don't transpose since eigen is column major
	if (principalAxes.determinant() < 0)
		principalAxes = -principalAxes;

	double *rotMat = principalAxes.data();
	mol.Rotate(rotMat);

	pair<double,double> split;
	split = calcSplitMoments(mol, 1);
	if(split.first < split.second)
	{
		//rotate around x
		double xrot[] = {1,0,0, 0,-1,0, 0,0,-1};
		mol.Rotate(xrot);
	}
	split = calcSplitMoments(mol, 0);
	if(split.first < split.second)
	{
		//rotate around y
		double yrot[] = {-1,0,0, 0,1,0, 0,0,-1};
		mol.Rotate(yrot);
	}
}

ShapeObj::ShapeObj(OBMol& mol, const MolInfo& info, float dimension,
		float resolution) :
		minfo(info)
{
	OBMol m = mol;
	m.DeleteHydrogens();
	normalizeMol(m);
	//set in parent
	set(m, dimension, resolution);
}

void ShapeObj::write(ostream& out) const
{
	out.write((char*) &minfo, sizeof(minfo));
}
