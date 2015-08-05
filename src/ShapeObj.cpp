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
pair<double,double> ShapeObj::calcSplitMoments(const vector<Eigen::Vector3d>& coords, unsigned C)
{
	using namespace Eigen;

	double above = 0;
	double below = 0;
	unsigned A = (C+1)%3;
	unsigned B = (C+2)%3;

	//find center of mass
	for(unsigned i = 0, n = coords.size(); i < n; i++)
	{
		const Vector3d& pt = coords[i];
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


//align coords to their moment of inertia in a canonical way
//mutates coords and sets translate and rotate
void ShapeObj::computeAndApplyNormalization(vector<Eigen::Vector3d>& coords, Eigen::Vector3d& translate, Eigen::Matrix3d& rotate)
{
	using namespace Eigen;
	//we need to normalize the molecular pose
	Vector3d center = Vector3d::Zero();
	Matrix3d I = Matrix3d::Zero();

	const unsigned n = coords.size();
	for(unsigned i = 0; i < n; i++)
	{
		center += coords[i];
	}
	center /= n;
	translate = center*-1.0;

	double Ixx = 0, Iyy = 0, Izz = 0, Ixy = 0, Ixz = 0, Iyz = 0;
	//calculate inertial tensor
	for (unsigned i = 0; i < n; i++)
	{
		coords[i] -= center;
		double x = coords[i][0];
		double y = coords[i][1];
		double z = coords[i][2];

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

	Matrix3d principalAxes = evecs.transpose();
	if (principalAxes.determinant() < 0)
		principalAxes = -principalAxes;

	rotate = principalAxes;
	for (unsigned i = 0; i < n; i++)
	{
		coords[i] = rotate*coords[i];
	}

	pair<double,double> split;
	split = calcSplitMoments(coords, 1);
	if(split.first < split.second)
	{
		//rotate around x
		Matrix3d xrot; xrot  << 1,0,0, 0,-1,0, 0,0,-1;
		rotate = xrot*rotate;
		for (unsigned i = 0; i < n; i++)
		{
			coords[i] = xrot*coords[i];
		}

	}
	split = calcSplitMoments(coords, 0);
	if(split.first < split.second)
	{
		//rotate around y
		Matrix3d yrot; yrot << -1,0,0, 0,1,0, 0,0,-1;
		rotate = yrot*rotate;
		for (unsigned i = 0; i < n; i++)
		{
			coords[i] = yrot*coords[i];
		}

	}
}

void ShapeObj::normalizeMol(OBMol& mol)
{
	using namespace Eigen;

	vector<Vector3d> coords;
	OBAtom *atom;
	vector<OBAtom*>::iterator j;
	//find center of mass
	for (atom = mol.BeginAtom(j); atom; atom = mol.NextAtom(j))
	{
		Vector3d c(atom->x(), atom->y(), atom->z());
		coords.push_back(c);
	}

	Vector3d translate;
	Matrix3d rotate;

	computeAndApplyNormalization(coords, translate, rotate);

	unsigned i = 0;
	for (atom = mol.BeginAtom(j); atom; atom = mol.NextAtom(j), i++)
	{
		Vector3d& c = coords[i];
		atom->SetVector(c[0], c[1], c[2]);
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

//use provided transformation to move molecule
ShapeObj::ShapeObj(OBMol& mol, const Eigen::Vector3d translate, const Eigen::Matrix3d& rotate, const MolInfo& info, float dimension,
		float resolution)
{
	OBMol m = mol;
	m.DeleteHydrogens();
	vector3 t(translate[0],translate[1],translate[2]);
	m.Translate(t);

	Eigen::Matrix3d trans = rotate.transpose(); //get out of column major order
	double *rot = trans.data();
	m.Rotate(rot);

	//set in parent
	set(m, dimension, resolution);
}

void ShapeObj::write(ostream& out) const
{
	out.write((char*) &minfo, sizeof(minfo));
}
